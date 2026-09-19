# CueMol風のPyMOLスタイル

[Package README](../../README.md) | [English](../pymol_cuemol.md)

`cuemol_style_in_pymol` はPyMOLに分子形状と操作中のGPU材質描画を追加する独立パッケージです。Qt版PyMOL 3.1、NumPy、SciPy、PyOpenGL、Pillow、および互換OpenGL 2.1 / GLSL 1.20コンテキストが必要です。`pixi install` でPyMOLを含む開発環境を導入できます。CueMolとmdtbxは不要です。PyMOLのソースコードと標準コマンドは変更しません。

## 大規模構造

既定の `atomic_mode=auto` は、2,000 原子以上のレイヤーで既定材質かつ輪郭なしの場合、球と結合を PyMOL の球・分割円柱命令で描きます。ポリマーのリボンは独自形状を使います。水・水素・脂質・イオンを含む選択原子を保持し、高速化のために原子を除外しません。元座標・色・選択・全読み込み状態・refresh・reset・標準 ray に対応します。

原子部分の照明と透明度は PyMOL 標準処理になり、独自材質シェーダーやサンプリングによる色合成とは異なります。大規模系でも従来の描画を使う場合は `atomic_mode=mesh`、小規模系で軽量表示を指定する場合は `atomic_mode=native` を使います。独自材質・輪郭・表面はメッシュ処理を継続し、非対応の明示的 native 指定はエラーになります。軽量形状と CGO もキャッシュ予算に含みます。選択原子や球・結合の半径は変えません。

軽量な半透明原子では、三角形への展開を避けるため PyMOL のシーン全体の設定 `transparency_mode=3` を使います。CueMol と Mol* は該当ビューの表示中、この設定を共有します。最後の軽量半透明ビューを reset すると、途中で明示的に変更していない限り以前のモードへ戻します。その間、他の半透明な標準オブジェクトにも同じ描画方式が適用されます。

```text
cuemol_style
cuemol_style ballstick, atomic_mode=native
cuemol_style ribbon, atomic_mode=mesh
```

読み込み・描画と準備時間を分けて計測できます。

```sh
uv run --no-project --python .pixi/envs/default/bin/python python tests/benchmark_large.py system.gro --output .cache/large-structure.json
```

原子数・状態数・準備時間・形状と CGO のサイズ・プロセス最大メモリを記録します。入力と計測結果はバージョン管理しません。

## 使い始める

このリポジトリで `pixi install` を実行するか、既存のPyMOL用Python環境に `uv pip install "git+https://github.com/th2ch-g/cuemol_style_in_pymol.git"` で導入します。PyMOLのPythonコンソールで次の行を実行するか、`.pymolrc.py` に追記してコマンドを登録します。

```python
from cuemol_style_in_pymol import __init_plugin__

__init_plugin__()
```

パッケージのインポートだけではPyMOLの設定を変更しません。mdtbxの `pymol_plugins` と連携する場合は自動登録されます。構造を読み込んだ後、PyMOLのコマンド欄で次のコマンドを実行します。

```text
cuemol_style
cuemol_style richardson
cuemol_style toon1
cuemol_style toon2, selection=chain A
cuemol_style matte, representation=surface, transparency=0.4
cuemol_style list
help cuemol_style
```

引数なしの `cuemol_style` は `ribbon` プリセットを適用します。既定のビュー名は `cuemol` です。同じ名前で別のスタイルを適用すると置き換わります。`richardson` は薄いリボン、シート矢印、明るいヘリックス裏面、色鉛筆状のハッチング、黒い輪郭線を使います。回転とズームはGPU描画に即座に反映されます。

`richardson` は暖色の紙、55度・-35度・80度の鉛筆線、分子色の線、線を描かないハイライトを使います。ヘリックスの平らな裏面だけを明るくし、丸い縁は元の色を保ちます。不透明な対話表示では鉛筆線と深度・法線による輪郭をGPUだけで計算します。3倍サンプリングのタイル描画で細い鉛筆線と滑らかな外周を保ち、回転やstate再生中のCPU画像計算とGPU画像転送を省きます。輪郭のない通常材質は直接描画します。

専用PNG/rayは精密な処理を維持します。法線と色を3倍で描画し、参照と同じハッシュ・ノイズ・線幅変化・サンプルごとの層の乗算をネイティブ処理で評価します。階調には深度フォグを含めます。GPU領域はタイルに分割します。輪郭は可視面の深度と法線から抽出し、遮蔽境界で連結して平滑化し、丸い線帯として描きます。内部の三角分割辺は輪郭抽出に使いません。

分子選択には元の構造を使い、`chimerax_style` が管理する描画用コピーは除外します。
セッション再読み込み後も、適用や refresh で描画用コピーを入力構造として扱いません。
通常の元構造は、名前がアンダースコアで始まる場合も利用できます。
各コマンドは独立したビューを管理します。同じ原子の表示を別のレンダラーへ切り替える
場合は、先に元のコマンドで reset し、その表示と設定を復元してください。

## スタイルと設定

スタイル名

| 分類 | 名前 |
| --- | --- |
| タンパク質の形状 | `richardson`, `ribbon`, `round_ribbon`, `fancy_ribbon`, `cartoon`, `round_cartoon`, `tube` |
| その他の分子形状 | `nucleic`, `ballstick`, `cpk`, `surface` |
| 照明と陰影 | `default`, `shadow`, `nolighting`, `matte`, `toon1`, `toon2` |
| 装飾材質 | `diff_metal`, `spec_metal`, `metallic_chrome`, `metallic_copper`, `stone35`, `wood31`, `wood14scl2` |
| 輪郭線 | `outline`, `silhouette` |

形状プリセットは名前に対応する表示形式を選びます。材質・輪郭線プリセットで表示形式を省略した場合は `ribbon` を使い、元がPyMOLのcartoon表示でもヘリックスを螺旋状に描きます。例えば `cuemol_style toon1` はリボン、`cuemol_style cartoon` はヘリックスを円柱にする表示です。明示的な表示形式の指定は既定値を上書きします。

材質・輪郭線スタイルに `representation=auto` を明示すると、既存のsticks・spheres・surface・cartoon・ribbonを独自形状で引き継ぎます。linesまたはnonbondedだけの原子は、タンパク質リボン、核酸主鎖と塩基対ロッド、または球棒モデルになります。この自動モードでは非表示原子を表示しません。マップ・ラベルなど未対応の表示はPyMOL標準のままです。

表示形式は `ribbon`、`cartoon`、`tube`、`nucleic`、`ballstick`、`sticks`、`cpk`、`surface` で指定できます。`cartoon` はヘリックスを円柱にし、`ribbon` は骨格に沿います。品質は `low`、`medium`（既定）、`high` です。高品質ほど準備時間とメモリ使用量が増えます。

`edge` には `auto`、`none`、`edges`、`silhouette`、`thin`、`normal`、`thick` を指定します。`edges` は自己遮蔽の輪郭も含み、`silhouette` は外周を描きます。折れ目の線は現行参照の既定値に合わせて無効です。`edge_width` の単位はオングストロームで、GPU画像では最小1ピクセルです。`edge_color` はPyMOLの色名を受け取ります。`transparency=keep` は元の表示形式の不透明度を保持し、0から1の数値はそれを上書きします。

```text
cuemol_style spec_metal, representation=cpk, edge=silhouette
cuemol_style richardson, edge_width=0.12, edge_color=black
cuemol_style wood31, representation=surface, quality=high
```

## 既定の配色

`color=cuemol` は、CueMol GUIで新しく読み込んだ分子の初期配色に合わせます。タンパク質のリボン・cartoon・tubeはDefaultHSCPaint、核酸形状はDefaultNucl、原子表現と表面はDefaultCPKColoringを使います。これらは分子の配色を参照するため、原子表現や表面の炭素にも分子の配色が引き継がれます。

CueMolの配色

| 対象 | 色 |
| --- | --- |
| ヘリックス / シート / コイル | khaki `#F0E68C` / SteelBlue `#4682B4` / FloralWhite `#FFFAF0` |
| 核酸の骨格と塩基 | 黄色 `#FFFF00` |
| 炭素 / 窒素 / 酸素 | 分子の配色 / 青 / 赤 |
| 水素 / 硫黄 / リン | シアン / 緑 / 黄 |

他のモードは `keep`（既存の原子色）、`chain`、`ss`（WoodyHSCの二次構造色）、`rainbow`、`element`（一般的な元素色）です。適用中にPyMOLの原子色を変更した場合、その色を使う配色モードで `refresh` を実行すると反映されます。元の原子色や二次構造の割り当ては変更しません。

```text
cuemol_style richardson, color=keep
cuemol_style nucleic, color=chain
```

## 選択とstate管理

独自形状をクリックすると、対応する元の原子を `sele` に選択します。原子・残基・鎖への展開はPyMOLのマウス選択モードに従い、Shiftで追加選択します。選択された原子位置にはピンクのマーカーを表示します。ドラッグによる回転・移動は標準操作を維持します。手前の不透明な標準形状は、背後の独自形状へのクリックを遮ります。クリック選択にはQt GUIが必要です。

再生前に読み込み済み全stateの基本メッシュを準備します。透明材質のサンプルは、視点・描画サイズ・フォグ・表示stateが変わると描画コールバックの外で再構築します。全体のstate変更、動画フレームからstateへの対応、オブジェクトごとのstate指定、`all_states` に対応します。座標、結合、色、二次構造、透明度、state数を編集した後は `refresh` が必要です。元の表示形式を変える場合はreset後に表示を変更し、スタイルを再適用します。元オブジェクトの表示状態とstate設定は短いQtタイマーで同期します。

```text
cuemol_style refresh
cuemol_style reset
cuemol_style reset, name=all
```

`reset` は元の標準表現を復元し、生成オブジェクトとGUIフックを削除します。適用、refresh、reset、スタイルの復元では背景を変更しません。背景はPyMOLの `bg_color` で指定します。元の座標や結合は変更しません。互いに重ならない選択には異なる `name` を使えます。管理対象が重なる選択は拒否します。管理グループや元オブジェクトを削除すると後処理を実行します。PyMOLセッションには設定を保存し、プラグインが利用可能なら読み込み時に形状を再生成します。

`cache_mb=2048` は準備済みメッシュ、標準CGOコマンドのデータ量、入力スナップショットの推定量にそれぞれ適用します。各分類に独立した上限があり、標準描画側の追加メモリも使うため、プロセス全体のメモリ上限ではありません。GPUキャッシュは最大256 MiBの頂点バッファを保持し、古いstateを追い出します。追い出されたバッファは準備済みCPUメッシュから再転送します。大きな軌跡、特に表面や半透明の全原子表示は準備の負荷が大きくなります。`cuemol_style list` は各表示の準備時間とメッシュ・標準CGOの保存量を表示します。全stateのray形状を保持する分、準備時間、メモリ使用量、保存セッションの容量が増えます。

対話表示では小さな三角形群を囲む範囲をメッシュとともに準備し、視点ごとの画面上の範囲を求めます。
形状のないタイルと合成領域を省き、静止時の3倍サンプリングは維持します。この範囲データも `cache_mb` に含みます。
鉛筆線は最大幅と傾きを考慮しても届かないサンプルの計算を省き、参照値を保ちます。
微小原子・回転・平行/透視投影で、省略の有無による画素一致を回帰テストします。

マウスでの回転・平行移動・ホイール操作中は1ピクセル1サンプルに切り替えます。Richardsonの3層の鉛筆線、配色、材質定数は保ち、移動中の輪郭のアンチエイリアスを軽くします。操作が150 ms止まると元の3倍表示へ戻ります。透明表示も操作中は準備済みの標準メッシュを使い、停止後に精密な視点サンプルへ戻します。PNG/rayは常に精密な経路を使います。大画面の測定例では静止用54.8 MiBに加え、再利用する操作用タイル24.1 MiBを保持します。

クリック判定はボタンを離した時に行い、ドラッグ開始時のGPU深度読み出しと三角形の交差判定を省きます。実際のクリックは事前計算した範囲で候補を絞ります。標準オブジェクト、Shift選択、Molstarとの同時使用も検証します。

対話表示のフレームバッファは3倍タイルの1サンプルあたり24 byteです。510ピクセルのタイルと両側2ピクセルの余白で約54.4 MiBを再利用します。太い輪郭では余白が増えます。上限は256 MiBとドライバのテクスチャ寸法です。

描画用タイルのGPU領域は最大約126 MiBで、これとは別に一時CPU配列を使います。輪郭の接続がタイル境界に依存しないよう、画面全体の3倍の深度・法線画像を使います。透明表示・管理外の形状を含む専用rayは可視面の1/3ピクセルの各サンプルを2枚の三角形にするため、基本メッシュより多くのメモリが必要です。非表示stateの密なCGOは粗い代替形状へ戻します。`cache_mb` 超過時は品質を下げず、エラーを報告します。

## 画像出力

```text
ray 2400, 1800
png figure_ray.png
cuemol_style png, filename=figure.png, width=2400, height=1800
cuemol_style ray, filename=figure_ray.png, width=640, height=480
```

標準の `ray` と `png, ray=1` は準備済み全stateの標準CGOを使い、ヘッドレス環境でも動作します。不透明メッシュはOpenGLに表示されないray専用三角形として保持します。透明CGOは重複させません。標準rayは粗い分子座標の材質サンプル、Richardsonの平均階調、およびユーザーのPyMOL照明を使います。ヘッドレス環境やQtタイマーを待たない呼び出しでは、表示・オブジェクトstateの変更後にrefreshしてください。

専用の `png` はQtを使い、背景を維持し、選択マーカーを除き、不透明形状を3倍サンプリングのタイルで描画します。透明な表示グループは、本体と輪郭をまとめた色・深度テクスチャを使い、不透明なGPUレイヤーとして個別に描画します。標準CGOによる余分な照明・フォグを輪郭に適用しません。同じ3倍格子を縮小した後、表示色空間でグループを合成します。転送バッファ、画素格納設定、フレームバッファ、描画領域、行列、シェーダー、GL属性は復元します。

専用の `ray` は管理対象の形状だけが見え、PyMOLのgammaが1の場合、同じ3倍の色・深度サンプルを直接合成します。大量の画素三角形への変換を省き、鉛筆線・輪郭・クリッピング・透明グループ・背景は保持します。標準の画素CGOはfloat32の交差判定で境界サンプルがわずかに変わるため、重なり回帰テストでは画像全体の平均差が0.1/255未満であることを確認します。

管理外の形状がある場合やgammaが1以外の場合は、視点に応じた画素CGOを作り、個々の鉛筆線とGPUと同じ連結済み輪郭を含めます。各可視サンプルには、その深度に手前向きの四角形を1枚置きます。隠れた三角形や別描画の輪郭円柱による透明度の重複を防ぎます。rayも同じ3倍格子とグループ合成を使い、ヘッドレスでも動作します。二重照明を避けるため、**出力シーン全体**の照明を一時的に中立化し、rayの影とフォグを無効にします。その画像内の通常オブジェクトも中立照明になります。設定・表示・再生状態は失敗時も復元します。通常オブジェクトの既存の陰影を保つ場合はGPU PNGを使ってください。

対話中の透明表示は、標準オブジェクトと共存するためPyMOLのCGOを使います。同じ不透明度の部品を一つの可視面としてサンプル化し、内部の重なりと輪郭色による不透明度の重複を防ぎます。色とフォグ補正は視点に追従し、元の透明度設定を維持します。対話表示でのオブジェクト間のソートと合成にはPyMOLの規則を使います。

専用PNG/ray出力はUmbreonと同じ `(1 - sum(alpha_i)) * B + sum(alpha_i * L_i)` で合成します。`B` は全透明グループを除いた背景シーン、`L_i` は背景シーンとグループ `i` を不透明に描いた画像です。陰影・輪郭・3倍画像の縮小後、RGBをsRGB空間で加算し、被覆alphaを線形に加算します。グループの不透明度の合計が1を超えると、背景の重みは負になります。管理外の形状も各レイヤーに含め、透明面の背後のオブジェクトとの遮蔽を反映します。失敗時も一時オブジェクトと表示状態を復元します。

オブジェクト変換行列、ステレオ・VRのクリック選択、原子のドラッグ編集、ヘッドレスでの対話GPU描画は対象外です。座標変換は元の原子へ適用してrefreshしてください。

## 表示形式の照合結果

参照は [CueMol 3173d8a](https://github.com/CueMol/cuemol2/tree/3173d8af62e211dd37b943ee53b3d6a632e6b5d7) と [Umbreon bf75c8a](https://github.com/CueMol/umbreon/tree/bf75c8adc05ed70a1344afbd718bcaab651c1070) です。数値定義と実際の出力画像を比較します。対象は現行Umbreonのdirect描画で、GI・AO・投影影を無効にします。旧OpenGL、POV-Rayの手続き的材質、GI描画は別方式であり、方式をまたいだ画素の完全一致は保証しません。

| 形状 | 参照寸法・構築方法 |
| --- | --- |
| `ribbon`, `round_ribbon` | ヘリックス半幅1.2、シート半幅1.4、半厚0.2、コイル半径0.35 A。弦長を節点間隔とする自然スプライン、シート基準点50%平滑化、運搬フレーム、断面変化の法線補正、接続gamma 2.2。矢印倍率1.8、gamma 2.2/1.2。 |
| `fancy_ribbon`, `richardson` | ヘリックス半幅1.3、縁の半径0.2、鋭さ0.3、シート半幅1.2、コイル半径0.25 A。平らなヘリックス裏面とシート側面だけHSV彩度を0.4下げます。矢印倍率1.6、gamma 1.0。 |
| `cartoon`, `round_cartoon` | 前後の残基を含む曲率罰則付き自然スプライン。ヘリックスrho 3.0、半径は平均軸距離+0.2 A。シート半幅/半厚1.4/0.2 A、rho 3.0/1.0、向きのrho 5.0。コイル半径0.2 A、rho -1/-2、重み付き端点とシート端微分の拘束。 |
| `tube`, `nucleic` | Tube半径0.35 A。核酸はP原子を基準点とする半軸1.25/0.5 Aの主鎖、半径0.5 Aの塩基対ロッド。スプライン端面は5段の半球。面内の互換水素結合から対を推定。 |
| `ballstick`, `sticks`, `cpk` | 球棒半径0.3/0.2 A、sticksは0.2/0.2 A。CPKのH/C/N/O/S/P半径は1.2/1.7/1.55/1.52/1.8/1.8 A、その他1.7 A。結合の中点で色を切り替え、密なメッシュで解析的な球・円柱を近似。 |
| `surface` | 現行CueMolの距離場SES。プローブ1.4 A、参照半径、原子球のSASとプローブ球の2段階contour、外向き勾配法線、内側成分の選別。low/medium/highのdetailは3/6/10、格子間隔は `1.43 / (1 + 0.2 * (detail - 1))` A。 |

表面は既定アルゴリズムがdistfieldとなったCueMol `af9509e` を対象に更新しました。格子・出力メモリを制限し、球の近傍だけを計算して改善しない距離の平方根を省きます。6400万格子点またはworkspace上限を超えた場合だけ格子間隔を広げます。旧EDTSurfは内部に残し、marching-cubes表の利用許諾を保持しています。[由来と利用許諾](../../native/edtsurf/README.md) を同梱しています。実行時のCueMol依存はありません。ソースからのビルドにはC++17コンパイラとpybind11が必要で、wheelには拡張を含めます。

リボン接続部は残基ごとのパラメータ表、断面変化の解析的な微分、平らな矢印の肩、固定断面、参照の端面形状を使います。メッシュ分割、鎖切断・altlocの選択、修飾塩基トポロジーには差が残ります。リボン接続面は共通化し、cartoonはRibbon2Rendererと同様に要素ごとに軸を求め、ヘリックス直前でもシート先端を細く保ちます。マップとラベルは標準表示のまま、元の原子・色・結合・二次構造を維持します。

## 材質と輪郭線の照合

PBR値は `UmbreonDisplayContext.cpp` に合わせています。ambient光は1、主光源は正規化した `(1,1,1)` で強度0.52、視線軸の補助光は鏡面反射なしで0.78です。GPUと標準サンプルはGGX、相関Smithマスキング、Schlick Fresnelを共通に使います。

| 材質 | Ambient | Diffuse | Metallic | Roughness | Specular | Reflection |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `default` | .20 | .80 | 0 | .3742032 | .40 | 0 |
| `matte` | .30 | .80 | 0 | .5 | 0 | 0 |
| `diff_metal` | .35 | .30 | 1 | .5491005 | .80 | .10 |
| `spec_metal` | .15 | .60 | 1 | .3742032 | .80 | .65 |
| `metallic_chrome` | .20 | .80 | 1 | .05 | .50 | 0 |
| `metallic_copper` | .20 | .80 | 1 | .15 | .50 | 0 |
| `stone35` | .20 | .80 | 0 | .85 | .25 | 0 |
| `wood31`, `wood14scl2` | .20 | .80 | 0 | .45 | .50 | 0 |

金属も分子色を保ち、背景を反射します。明示reflectionが0ならPBRのF0環境項を使います。木と石は現行backendのPBR材質を使い、旧POV-Rayの模様は描きません。`toon1` はdiffuse .8・brilliance 0、`toon2` はambient .3・diffuse .5・Phong 10000・size 50です。`shadow` はambient .75、`nolighting` は1です。toonとmetalの各組は異なる見た目になります。GI、AO、投影影、シーン内反射は未実装です。

`thin/normal/thick` 幅は .03/.06/.15 A、最小全幅は出力画像の1ピクセルです。遮蔽検査と深度フォグを使います。`outline` は自己遮蔽の輪郭も描き、`silhouette` は内側の輪郭を抑えます。折れ目の線は無効です。深度差12ピクセル、強弱の輪郭接続、法線の連続性、メッシュ内の線分交差検査を使い、連続した面の折れ目を除外します。4出力ピクセル未満の線を除き、遮蔽境界の連続線をつないでからChaikin平滑化2回と .4出力ピクセルの単純化を行います。丸い線帯は外側配置で、内側に半出力ピクセルの余白を持ちます。専用GPU PNG、透明サンプル、専用rayはこの共通処理を使い、三角分割の内部対角線には依存しません。

Richardsonの紙色は `#F0ECDD`、角度55/-35/80、閾値 .92/.62/.34、インク倍率1/.74/.38、筆圧下限 .4、明度差の下限 .15です。3倍時の実効間隔は2/3出力ピクセル、全幅 .45、線長/間隙50/5、幅の揺れ .45、長さの揺れ .5、先細り .35、角度の揺れ5度、紙の粒度 .15・scale 3です。ハッシュ・ノイズ・線の被覆率には独立した参照値テストがあります。各サンプルで層を乗算してから平均します。階調はdiffuse .85、ambient .05、wrap .5、rim power 3.5、rim bias .35、white point 1.2、gamma 2.4、ハイライト遷移 .81〜.86です。`GL_EXT_gpu_shader4` に対応した環境では、対話用GLSLも参照の整数ハッシュ・ノイズ・線幅変化・層の乗算を同じ3倍グリッドで評価します。非対応環境では紙色・角度・閾値・インク倍率・階調を保ち、浮動小数のノイズで近似します。対話表示の輪郭は局所的な深度と接平面の連続性を使い、専用出力の線の追跡や平滑化を近似します。保持用の標準ray形状は鉛筆の平均近似を使います。

## 検証

```sh
pixi install --locked
pixi run check
uv run --no-project --python .pixi/envs/default/bin/python python -m pytest tests -q
uv run --no-project --python .pixi/envs/default/bin/python python tests/check_cuemol_style.py --gui --benchmark --output .cache/validation
uv build --python .pixi/envs/default/bin/python
```

画像比較ではタンパク質PDBと互換CueMol Nodeモジュール・設定ファイルの場所を指定します。これらは検証専用です。

```sh
uv run --no-project --python .pixi/envs/default/bin/python python tests/audit_appearance.py \
  --structure structure.pdb --output .cache/appearance \
  --reference-module "$CUEMOL_MODULE" --reference-config "$CUEMOL_CONFIG"
uv run --no-project --python .pixi/envs/default/bin/python python tests/compare_appearance.py .cache/appearance
```

本家側の色はプラグインのRGB配列で上書きせず、GUIの初期paintingとDefaultHSCPaint/DefaultCPKColoringから取得します。全26プロファイルと明示的sticksの計8形状で、座標、二次構造、平行/透視投影カメラ、画像寸法、renderer設定をそろえます。Richardsonは両側とも不透明の紙色背景、それ以外は白背景です。プラグイン自体は背景を変更しません。`--angle`、`--zoom`、`--perspective`、`--transparency`、`--live`、`--ray` で条件を変えられます。`--live` は専用の精密出力の代わりに対話表示を取得します。実行版・モジュールのdigest・manifest・ログ・画像・数値差・比較一覧をignore対象の出力先へ保存します。前景IoUは陰影も含む指標で、純粋な形状精度ではありません。画像の位置合わせ処理は行いません。

以下の旧版での精密出力比較にはCueMol 2.3.13.523（`aeacb41`）と上記の定義を使いました。このbuildから `3173d8a` までの対象rendererの変更は選択用IDの追加で、今回使うdirect exporter、材質表、EDTSurf、形状寸法、スプライン計算は同じです。実行モジュールのdigestと参照ソースのrevisionは区別して保存します。

代表的な実測値は次のとおりです。MAEは両画像の前景領域の和集合におけるRGB各色の平均絶対差（0〜255）、平滑化MAEは2ピクセルのGaussianを適用した値です。CueMol、GPU PNG、専用rayとも3倍サンプリングです。タンパク質は1CRN・mediumで、カメラ・色・二次構造をそろえています。

| 出力・条件 | プロファイル | 前景IoU | MAE | 平滑化MAE |
| --- | --- | ---: | ---: | ---: |
| GPU・不透明・640×480 | `cartoon` | 0.9987 | 0.56 | 0.30 |
| GPU・不透明・640×480 | `default` | 0.9998 | 0.25 | 0.16 |
| GPU・不透明・640×480 | `richardson` | 0.9967 | 0.72 | 0.18 |
| GPU・透視投影・回転-45度・zoom .8・640×480 | `richardson` | 0.9893 | 3.28 | 0.51 |
| GPU・透明度.25・320×240 | `richardson` | 0.9889 | 2.29 | 0.41 |
| GPU・透明度.65・回転37度・zoom 1.2・320×240 | `richardson` | 0.9918 | 1.30 | 0.33 |
| 専用ray・不透明・320×240 | `default` | 0.9963 | 0.72 | 0.47 |
| 専用ray・不透明・320×240 | `richardson` | 0.9898 | 2.47 | 0.52 |
| 専用ray・透明度.4・320×240 | `richardson` | 0.9913 | 2.18 | 0.42 |
| 専用ray・透視投影・回転-45度・zoom .8・320×240 | `richardson` | 0.9875 | 4.19 | 0.99 |

輪郭位置とrenderer単位の透明合成には、上記の共通画面処理を使います。表は指定したカメラとシーンでの画像誤差で、任意の入力に対する画素単位の一致を保証するものではありません。

対話表示も、同じ条件の入力を使ってCueMol 2.3.15.530（`af9509e`）と直接比較しました。輪郭を近似するため、専用の精密出力より差が大きく、特に遮蔽境界の接続に差が残ります。

| 対話表示の条件 | プロファイル | 前景IoU | MAE | 平滑化MAE |
| --- | --- | ---: | ---: | ---: |
| 不透明・640×480 | `toon1` | 0.9830 | 6.34 | 1.15 |
| 不透明・640×480 | `toon2` | 0.9562 | 8.16 | 2.15 |
| 不透明・640×480 | `richardson` | 0.9719 | 9.18 | 1.33 |
| 不透明・640×480 | `silhouette` | 0.9870 | 5.03 | 0.99 |
| 透視投影・回転-45度・zoom .8・640×480 | `toon1` | 0.9768 | 9.05 | 1.65 |
| 透視投影・回転-45度・zoom .8・640×480 | `richardson` | 0.9654 | 11.50 | 2.00 |

追加のGUI回帰検証は次のコマンドで実行します。

```sh
uv run --no-project --python .pixi/envs/default/bin/python python tests/check_appearance_gui.py --output .cache/appearance-regression
```

追加の外観検証では、平行・透視投影で対話表示のタイル分割有無を比較し、fallback shaderのコンパイルと、GPU鉛筆色の参照サンプラーとの数値比較を行います。鉛筆色の平均誤差はRGBA8への量子化を含め0.154/255でした。

GUI検証には2400×1800出力、標準/独自の透明表示、選択、state対応、セッション再読込、失敗時の復元を含めます。追加の回帰検証ではsheetの12方向回転、タイル間の輪郭連続性、転送バッファの復元、透明な複数グループの重なりも確認します。重なり画像のGPU/ray間の前景MAEは0.23/255でした。ベンチマークは500残基・100合成state・1280×720で、準備・warmup時間、CPU/GPU保存量、最大RSS、回転FPS、state切替FPS、実動画描画速度を記録します。回転とstate切替の計測では、Qtがまとめる可能性のある再描画要求ではなく、完了したOpenGL描画を数えます。不透明な対話表示はCPUの輪郭・鉛筆処理を使いません。密な透明サンプルと精密な書き出しは引き続き高コストです。再生成可能な画像・ビルド・レポートはignoreし、README用gallery PNGだけを追跡して `tests/render_gallery.py` で再生成します。

個別のスタイルだけを測定する場合は、次を実行します。

```sh
uv run --no-project --python .pixi/envs/default/bin/python python tests/check_cuemol_style.py --gui --benchmark-only --benchmark-style richardson --output .cache/benchmark-richardson
```

`--benchmark-style toon1` や `ribbon` でも同じデータを測定できます。レポートには準備時間、最大RSS、メッシュ・標準CGO・GPUの保存量、実際の回転・再生速度を記録します。性能はOpenGLドライバと画面内の形状面積に依存します。環境固有の結果は `.cache/` に保存します。

現在のNumPy 2.5.3検証は1CRNを500残基へ繰り返し、100合成stateを準備します。Richardsonは準備34.51秒、1280×720で回転72.0 FPS、state切替61.2 states/s、実動画27.6 states/sでした。メッシュ・標準CGO・GPUは660.7/1984.7/254.3 MiB、静止用framebufferは54.4 MiB、最大RSSは5204.0 MiBです。100 stateの事前読み込みには依然大きなメモリが必要です。この通常品質の回転計測と、以下のマウス操作中の計測は区別しています。

マウス操作の測定は `uv run --no-project --python .pixi/envs/default/bin/python python tests/benchmark_drag.py cuemol toon1 richardson --structure local_protein.pdb` で実行できます。

NumPy 2.5.3で全26プロファイルとsticksを本家CueMol 2.3.15.530（af9509e）の既定色と再比較しました。surfaceの前景IoUは0.8304から0.9919、RGB平均差は60.55から1.13/255へ改善しています。1CRN・mediumの形状生成は約12 ms（旧EDTSurfは約10 ms）で、現行方式に合わせた2段階contourを使います。

metallic_chrome/copperはUmbreonの本家設定と分子色を保持し、粗さは0.05/0.15です。このbackendのcopperは自動的に橙色になりません。対話表示のRGB平均差は0.82/0.85、条件をそろえた精密PNGでは0.030/0.031です。旧POV-Ray材質とは定義が異なります。

OpenGLの整数uniformへPython整数を渡し、NumPy 2.5.3のnumpy.boolエラーを修正しました。クリック判定もpixel中心と3倍描画の深度範囲に合わせ、見えている面の選択と管理外の前景による遮蔽を両立しています。uniform位置はGL contextごとにcacheし、resetやcontext交換で破棄します。1GGG・3200×1800のマウス回転ではtoon1が28.7→111.0 FPS、Richardsonが15.2→58.5 FPSでした。カメラが変化しGPU描画も完了したフレームを数え、Qtへの再描画要求数は使いません。押下処理は52〜55 msから0.1 ms未満へ短縮しました。準備時間0.89/0.63秒、メッシュ16.0/17.0 MiB、GPU頂点領域14.0/14.8 MiBです。専用rayは640×480・標準4 threadsで2.61/2.34秒から0.80/0.88秒へ短縮しました。速度はsceneとドライバに依存します。全26プロファイルで、ドラッグ後に元の静止画像へ画素単位で戻ることを確認しました。

ray高速化後もCueMol 2.3.15.530と同じカメラ・320×240で再比較しました。前景MAEはmetallic_copper/metallic_chromeが0.086/0.046、Richardsonが2.284、toon1が2.982でした。金属の材質定数と分子色は変更していません。
