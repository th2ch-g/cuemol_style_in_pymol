"""Compare direct sample composition with native pixel triangles in real PyMOL."""

from unittest.mock import patch

import numpy as np
from PIL import Image

from cuemol_style_in_pymol import cuemol_style as style
from cuemol_style_in_pymol import ray
from cuemol_style_in_pymol.controller import manager_for
from cuemol_style_in_pymol.export import settings


def compare_ray_paths(cmd, output):
    output.mkdir(parents=True, exist_ok=True)
    report = {}
    with settings(
        cmd,
        {
            "gamma": 1,
            "orthoscopic": 1,
            "ray_orthoscopic": -1,
            "ray_opaque_background": 1,
            "bg_gradient": 0,
        },
    ):
        for name, opacity, ortho, transparent in (
            ("toon1", 0, 1, False),
            ("richardson", 0, 0, False),
            ("metallic_copper", 0, 1, True),
            ("toon1", 0.4, 0, True),
            ("richardson", 0.4, 1, False),
            ("metallic_chrome", 0.4, 0, False),
        ):
            style("reset", name="all", quiet=1, _self=cmd)
            cmd.delete("all")
            cmd.pseudoatom("front", pos=[-0.6, 0, 0.4], elem="O", vdw=1.5)
            cmd.pseudoatom("back", pos=[0.6, 0, -0.4], elem="N", vdw=1.5)
            cmd.hide("everything")
            cmd.reset()
            cmd.zoom(buffer=2)
            cmd.bg_color("white")
            cmd.set("orthoscopic", ortho)
            cmd.set("ray_opaque_background", int(not transparent))
            style(
                name,
                selection="front",
                representation="cpk",
                name="front_style",
                transparency=opacity,
                quiet=1,
                _self=cmd,
            )
            style(
                "matte",
                selection="back",
                representation="cpk",
                name="back_style",
                transparency=opacity / 2,
                quiet=1,
                _self=cmd,
            )
            manager = manager_for(cmd)
            assert ray.can_compose(manager)
            arrays = []
            for fast in (True, False):
                path = output / f"{name}-{opacity}-{'fast' if fast else 'native'}.png"
                with patch.object(ray, "can_compose", return_value=fast):
                    style("ray", filename=path, width=160, height=120, _self=cmd)
                arrays.append(np.asarray(Image.open(path).convert("RGBA"), dtype=int))
            delta = abs(arrays[0] - arrays[1])
            # Native float32 pixel triangles can lose boundary subsamples.
            # Interior colors and the complete image must remain equivalent.
            assert delta.mean() < 0.1, (name, opacity, ortho, delta.mean())
            assert np.quantile(delta.max(axis=2), 0.99) <= 1
            assert delta.max() <= 64
            report[f"{name}-{opacity}"] = float(delta.mean())
            cmd.pseudoatom("unmanaged", pos=[0, 0, 1], vdw=0.3)
            cmd.show("spheres", "unmanaged")
            assert not ray.can_compose(manager)
            cmd.delete("unmanaged")
            cmd.set("gamma", 1.5)
            assert not ray.can_compose(manager)
            cmd.set("gamma", 1)
        style("reset", name="all", quiet=1, _self=cmd)
        cmd.delete("all")
    return report
