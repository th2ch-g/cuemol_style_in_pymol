/* Render a manifest through an explicitly supplied CueMol native module. */
const fs = require('node:fs');
const path = require('node:path');
const crypto = require('node:crypto');

const [modulePath, configPath, manifestPath] = process.argv.slice(2);
if (!manifestPath) throw new Error('Usage: node appearance_reference.cjs MODULE CONFIG MANIFEST');
const core = require(path.resolve(modulePath));
core.initCueMol(path.resolve(configPath));
const manifest = JSON.parse(fs.readFileSync(manifestPath));
const directory = path.dirname(path.resolve(manifestPath));
const invoke = (object, method, ...args) => object.invokeMethod(method, ...args);
const sceneManager = core.getService('SceneManager');
const streams = core.getService('StreamManager');
const scene = invoke(sceneManager, 'createScene');
const selection = expression => {
  const value = core.createObj('SelCommand');
  if (!invoke(value, 'compile', expression, scene.getProp('uid')))
    throw new Error(`Invalid selection: ${expression}`);
  return value;
};
const color = rgb => {
  const value = core.createObj('Color');
  invoke(value, 'setRGB', ...rgb);
  return value;
};
const reader = invoke(streams, 'createHandler', 'pdb', 0);
reader.setProp('build2ndry', false);
const molecule = invoke(reader, 'createDefaultObj');
invoke(reader, 'setPath', path.join(directory, manifest.structure));
invoke(reader, 'attach', molecule);
invoke(reader, 'read');
invoke(reader, 'detach');
invoke(scene, 'addObject', molecule);
const colors = new Map();
for (const atom of manifest.atoms) {
  const target = invoke(molecule, 'getAtom', atom.chain || '_', atom.resi, atom.name);
  if (!target) throw new Error(`Missing reference atom: ${atom.chain}.${atom.resi}.${atom.name}`);
  const position = core.createObj('Vector');
  invoke(position, 'set3', ...atom.position);
  target.setProp('pos', position);
  const id = target.getProp('id');
  const key = JSON.stringify(atom.color);
  if (!colors.has(key)) colors.set(key, []);
  colors.get(key).push(id);
}
const secondary = {};
for (const type of ['helix', 'sheet']) {
  const selected = invoke(molecule, 'getSelArray', selection(`rprop secondary=${type} and name CA`));
  secondary[type] = [];
  for (let i = 0; i < selected.getProp('elemCount'); ++i) {
    const atom = invoke(molecule, 'getAtomByID', invoke(selected, 'getAt', i));
    secondary[type].push(`${atom.getProp('chainName')}.${atom.getProp('residIndex')}`);
  }
}
fs.writeFileSync(path.join(directory, manifest.output.replace('.png', '-secondary.json')), JSON.stringify(secondary, null, 2));
for (const [kind, code] of [['helix', 'H'], ['sheet', 'S']]) {
  const expected = manifest.atoms.filter(a => a.name === 'CA' && a.ss === code).map(a => `${a.chain || '_'}.${a.resi}`);
  if (JSON.stringify(secondary[kind]) !== JSON.stringify(expected))
    throw new Error(`Reference secondary structure differs: ${kind}: ${secondary[kind]} != ${expected}`);
}
invoke(molecule, 'fireAtomsMoved');
const paint = core.createObj('PaintColoring');
if (manifest.native_colors) {
  // CueMol GUI's initial molecular paint, independently of the plugin colors.
  const styleManager = core.getService('StyleManager');
  for (const [sel, name] of [['sheet', 'SteelBlue'], ['helix', 'khaki'],
                             ['nucleic', 'yellow'], ['*', 'FloralWhite']])
    invoke(paint, 'append', selection(sel), invoke(styleManager, 'compileColor', name, scene.getProp('uid')));
} else {
  for (const [rgb, ids] of colors)
    invoke(paint, 'append', selection(`aid ${ids.join(',')}`), color(JSON.parse(rgb)));
}
molecule.setProp('coloring', paint);
const renderer = invoke(molecule, 'createRenderer', manifest.renderer);
if (manifest.native_colors) {
  const coloringStyle = ['ribbon', 'cartoon', 'tube', 'nucl'].includes(manifest.renderer)
    ? 'DefaultHSCPaint' : 'DefaultCPKColoring';
  invoke(renderer, 'applyStyles', [manifest.styles, coloringStyle].filter(Boolean).join(','));
} else {
  if (manifest.styles) invoke(renderer, 'applyStyles', manifest.styles);
  renderer.setProp('coloring', paint);
}
for (const [key, value] of Object.entries(manifest.properties || {})) renderer.setProp(key, value);
scene.setProp('bgcolor', color(manifest.background));
const camera = core.createObj('Camera');
const rotation = core.createObj('Quat');
for (const [index, key] of ['x', 'y', 'z', 'a'].entries()) rotation.setProp(key, manifest.rotation[index]);
camera.setProp('rotation', rotation);
const center = core.createObj('Vector');
invoke(center, 'set3', ...manifest.center);
camera.setProp('center', center);
camera.setProp('zoom', manifest.zoom);
camera.setProp('distance', manifest.distance);
camera.setProp('slab', manifest.slab);
camera.setProp('perspec', !!manifest.perspective);
invoke(scene, 'setCamera', 'audit', camera);
const backend = manifest.backend || 'umbreon';
const exporter = invoke(streams, 'createHandler', backend, 2);
if (backend === 'umbreon') {
  const settings = invoke(scene, 'getCreateAppData', 'render', 'RenderSettings');
  invoke(exporter, 'applyRenderSettings', settings, manifest.hatching ? 'umbreon_npr' : 'umbreon');
  exporter.setProp('useGI', false);
  exporter.setProp('aoSamples', 0);
  exporter.setProp('shadows', false);
  exporter.setProp('lightIntensity', manifest.hatching ? 1.55 : 1.3);
  exporter.setProp('flashFraction', 0.6);
  exporter.setProp('ambientFraction', manifest.hatching ? 0.16 : 0.0);
  exporter.setProp('transparentBackground', manifest.transparent_background ?? !!manifest.hatching);
  exporter.setProp('supersample', 3);
}
for (const [key, value] of Object.entries({camera: 'audit', width: manifest.width,
  height: manifest.height, perspective: !!manifest.perspective})) exporter.setProp(key, value);
invoke(exporter, 'attach', scene);
invoke(exporter, 'setPath', path.join(directory, manifest.output));
invoke(exporter, 'write');
invoke(exporter, 'detach');
fs.writeFileSync(path.join(directory, 'reference-runtime.json'), JSON.stringify({
  version: sceneManager.getProp('version'), build: sceneManager.getProp('build'),
  module_sha256: crypto.createHash('sha256').update(fs.readFileSync(modulePath)).digest('hex'),
  source_reference_revision: manifest.reference_revision,
}, null, 2) + '\n');
