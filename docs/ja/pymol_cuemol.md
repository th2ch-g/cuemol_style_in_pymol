# CueMol風のPyMOLスタイル

[Package README](../../README.md) | [English](../pymol_cuemol.md)

`cuemol_style_in_pymol` はPyMOLに分子形状と操作中のGPU材質描画を追加する独立パッケージです。Qt版PyMOL 3.1、NumPy、SciPy、PyOpenGL、および互換OpenGL 2.1 / GLSL 1.20コンテキストが必要です。`pixi install` でPyMOLを含む開発環境を導入できます。CueMolとmdtbxは不要です。PyMOLのソースコードと標準コマンドは変更しません。

## 使い始める

このリポジトリで `pixi install` を実行するか、既存のPyMOL用Python環境に `uv pip install "git+https://github.com/th2ch-g/cuemol_style_in_pymol.git"` で導入します。PyMOLのPythonコンソールで次の行を実行するか、`.pymolrc.py` に追記してコマンドを登録します。

```python
from cuemol_style_in_pymol import __init_plugin__

__init_plugin__()
```

パッケージのインポートだけではPyMOLの設定を変更しません。mdtbxの `pymol_plugins` と連携する場合は自動登録されます。構造を読み込んだ後、PyMOLのコマンド欄で次のコマンドを実行します。

```text
cuemol_style richardson
cuemol_style toon1
cuemol_style toon2, selection=chain A
cuemol_style matte, representation=surface, transparency=0.4
cuemol_style list
help cuemol_style
```

既定のビュー名は `cuemol` です。同じ名前で別のスタイルを適用すると置き換わります。`richardson` は薄いリボン、シート矢印、明るいヘリックス裏面、色鉛筆状のハッチング、黒い輪郭線と折れ目を使います。回転とズームはGPU描画に即座に反映されます。

`richardson` のGPUシェーダーはCueMol 3のRichardson階調設定を参照します。暖色の紙、55度・-35度・80度の不均一な鉛筆線、分子色を使う線、線を描かないハイライトを組み合わせます。階調は表面法線と視線方向で決まります。ヘリックス外側は元の色を保ち、内側は明るくします。操作中のシェーダーはCueMol既定の3倍サンプリングに対応する細かな鉛筆線をフィルターします。独自の乱数を使い、Umbreonの遮蔽・深度フォグ・画面上の輪郭線処理を省略するため、操作速度を重視した近似です。合わせた設定値と残る差は、材質と輪郭線の照合欄に記載しています。

## スタイルと設定

スタイル名

| 分類 | 名前 |
| --- | --- |
| タンパク質の形状 | `richardson`, `ribbon`, `round_ribbon`, `fancy_ribbon`, `cartoon`, `round_cartoon`, `tube` |
| その他の分子形状 | `nucleic`, `ballstick`, `cpk`, `surface` |
| 照明と陰影 | `default`, `shadow`, `nolighting`, `matte`, `toon1`, `toon2` |
| 装飾材質 | `diff_metal`, `spec_metal`, `metallic_chrome`, `metallic_copper`, `stone35`, `wood31`, `wood14scl2` |
| 輪郭線 | `outline`, `silhouette` |

形状プリセットは表示形式を指定します。材質プリセットは `representation=auto` を使い、既存のsticks・spheres・surface・cartoon・ribbonを独自形状で引き継ぎます。linesまたはnonbondedだけの原子は、タンパク質リボン、核酸主鎖と塩基対ロッド、または球棒モデルになります。自動モードでは非表示原子を表示しません。マップ・ラベルなど未対応の表示はPyMOL標準のままです。

表示形式は `ribbon`、`cartoon`、`tube`、`nucleic`、`ballstick`、`sticks`、`cpk`、`surface` で指定できます。`cartoon` はヘリックスを円柱にし、`ribbon` は骨格に沿います。品質は `low`、`medium`（既定）、`high` です。高品質ほど準備時間とメモリ使用量が増えます。

`edge` には `auto`、`none`、`edges`、`silhouette`、`thin`、`normal`、`thick` を指定します。`edges` は見えている鋭い折れ目も描きます。`edge_width` の単位はオングストロームで、GPU画像では最小1ピクセルです。`edge_color` はPyMOLの色名を受け取ります。`transparency=keep` は元の表示形式の不透明度を保持し、0から1の数値はそれを上書きします。

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

再生前に読み込み済みの全stateを準備します。全体のstate変更、動画フレームからstateへの対応、オブジェクトごとのstate指定、`all_states` に対応します。座標、結合、色、二次構造、透明度、state数を編集した後は `refresh` が必要です。元の表示形式を変える場合はreset後に表示を変更し、スタイルを再適用します。元オブジェクトの表示状態とstate設定は短いQtタイマーで同期します。

```text
cuemol_style refresh
cuemol_style reset
cuemol_style reset, name=all
```

`reset` は元の標準表現を復元し、生成オブジェクトとGUIフックを削除します。適用、refresh、reset、スタイルの復元では背景を変更しません。背景はPyMOLの `bg_color` で指定します。元の座標や結合は変更しません。互いに重ならない選択には異なる `name` を使えます。管理対象が重なる選択は拒否します。管理グループや元オブジェクトを削除すると後処理を実行します。PyMOLセッションには設定を保存し、プラグインが利用可能なら読み込み時に形状を再生成します。

`cache_mb=2048` は準備済みメッシュ、標準CGOコマンドのデータ量、入力スナップショットの推定量にそれぞれ適用します。各分類に独立した上限があり、標準描画側の追加メモリも使うため、プロセス全体のメモリ上限ではありません。GPUキャッシュは最大256 MiBの頂点バッファを保持し、古いstateを追い出します。追い出されたバッファは準備済みCPUメッシュから再転送します。大きな軌跡、特に表面や半透明の全原子表示は準備の負荷が大きくなります。`cuemol_style list` は各表示の準備時間とメッシュ・標準CGOの保存量を表示します。全stateのray形状を保持する分、準備時間、メモリ使用量、保存セッションの容量が増えます。

## 画像出力

```text
ray 2400, 1800
png figure_ray.png
png figure_ray.png, width=2400, height=1800, ray=1
cuemol_style png, filename=figure.png, width=2400, height=1800
cuemol_style ray, filename=figure_ray.png, width=2400, height=1800
```

標準の `ray` と `png, ray=1` をそのまま使えます。ヘッドレス環境にも対応します。不透明メッシュはPyMOLのray専用三角形命令を使うCGOとしても保持します。この三角形はOpenGLには表示されず、クリック選択にも干渉しません。半透明メッシュにはすでに標準CGOがあるため、重複させません。準備済みの全stateにray形状が対応し、全体のstateと動画フレームの変更は即座に反映されます。元オブジェクトの表示状態とオブジェクトごとのstate指定は、対話表示と同じQtタイマーで同期します。ヘッドレス環境でこれらを変更した場合や、Qtタイマーを待たずにrayを実行する場合は、先に `cuemol_style refresh` を実行します。

2つの `cuemol_style` 出力操作は表示中のシーン全体を保存し、ファイル名が必須です。`png` はGPUの見た目を取得するためQt GUIが必要です。`ray` は現在のstateと視点に対する一時CGO形状を作り、輪郭・折れ目の線と視点に応じた材質サンプルを加えます。ヘッドレスPyMOLでも使えます。描画失敗時も表示状態、再生、一時設定を復元します。選択マーカーは出力から除きます。

不透明な本体は独立したGPUシェーダーで描画します。透明な本体は法線による照明を頂点色に焼き込んだ標準CGOを使い、PyMOLの透明オブジェクトと合成します。この照明は回転中も分子座標に固定されます。標準rayは同じ分子座標での材質評価と、ユーザーのPyMOL照明・輪郭設定を使います。専用rayはカメラ座標で材質を評価し、円柱で輪郭線を加えます。どちらのrayもGPUの陰影や線の見え方とは異なります。Richardsonではrayと透明CGOに色鉛筆線の平均被覆率を頂点階調として適用します。個々のハッチング線は不透明GPU描画と `cuemol_style png` で表示されます。木・石・金属は手続き的な近似であり、CueMolのPOV-Rayテクスチャを完全には再現しません。`shadow` は一定の陰影を使う材質であり、シーンの投影影を生成する機能ではありません。

オブジェクト変換行列、ステレオ・VRのクリック選択、ドラッグによる原子編集、ヘッドレス環境での対話GPU描画は対応範囲外です。必要な変換は元の原子座標に適用してrefreshしてください。

## 表示形式の照合結果

形状の既定値は [CueMol 3173d8a](https://github.com/CueMol/cuemol2/tree/3173d8af62e211dd37b943ee53b3d6a632e6b5d7) の `default_style.xml`、`TubeSection`、`RibbonRenderer`、`Ribbon2Renderer`、`NARenderer`、原子描画実装と照合しました。Richardsonの階調パラメータは [Umbreon bf75c8a](https://github.com/CueMol/umbreon/tree/bf75c8adc05ed70a1344afbd718bcaab651c1070) と照合しています。ソース上の定義を確認したものであり、画像の完全一致を保証するものではありません。

- `ribbon` と `round_ribbon`：ヘリックス半幅1.2、シート半幅1.4、半厚0.2、コイル半径0.35オングストローム。軸は点間距離を節点間隔とする自然3次スプラインを使い、シートの基準点を50%平滑化します。シート矢印の拡大率は1.8、gammaはそれぞれ2.2と1.2です。
- `fancy_ribbon` と `richardson`：ヘリックス半幅1.3、円弧状の縁の半径0.2、鋭さ0.3、シート半幅1.2、コイル半径0.25。ヘリックス裏面とシート側面はHSV彩度を0.4下げます。シート矢印の拡大率は1.6、gammaは1.0です。
- `cartoon` と `round_cartoon`：曲率罰則付き自然スプラインでヘリックス軸を求め、rhoは3.0です。円柱半径は基準点と軸の平均距離に0.2を加えます。シート半幅・半厚は1.4・0.2、平滑化rhoはそれぞれ3.0と1.0、コイル半径は0.2です。シート矢印の拡大率は1.8、gammaは1.0です。接続部と端点の拘束条件は近似です。
- `tube`：半径0.35、自然3次スプライン軸。`nucleic`：P原子を基準点とし、半軸1.25・0.5の楕円主鎖と半径0.5の塩基対ロッドを使います。塩基対は互換性のある面内の水素結合接触から推定します。対の割り当てや修飾塩基への対応は、CueMolの残基トポロジーや塩基対情報と異なる場合があります。
- `ballstick`：全原子の半径0.3、結合半径0.2。`sticks`：原子・結合とも半径0.2。結合色は中点で明確に切り替わります。`cpk`：H/C/N/O/S/Pの半径は1.2/1.7/1.55/1.52/1.8/1.8、その他の元素は1.7です。メッシュの分割数は本プラグインのquality設定に依存します。
- `surface`：プローブ半径1.4オングストローム、上記と同じ元素半径を使う溶媒排除表面です。PyMOLの表面生成法はCueMolのEDTSurf/MeshMSと異なるため、三角形分割や細部は一致しません。

スプラインの向き、鎖切断判定、断面間の接続、端面、表面と原子の対応付けには本プラグイン独自の実装を使います。既定寸法を合わせても、すべての分子表現が同一になるわけではありません。マップとラベルは既存のPyMOL標準表示を維持します。

## 材質と輪郭線の照合

OpenGL材質の係数はCueMolの `default_style.xml` に合わせています。

| 材質 | Ambient | Diffuse | Specular | Shininess |
| --- | ---: | ---: | ---: | ---: |
| `default` | 0.2 | 0.8 | 0.0 | 32.0 |
| `shadow` | 0.75 | 0.0 | 0.0 | 0.0 |
| `nolighting` | 1.0 | 0.0 | 0.0 | 0.0 |
| `matte` | 0.3 | 0.6 | 0.0 | 32.0 |
| `toon1` | 0.0 | 0.85 | 0.0 | 0.0 |
| `toon2` | 0.0 | 0.85 | 0.0 | 32.0 |
| `diff_metal`, `spec_metal` | 0.2 | 0.5 | 0.7 | 76.8 |

GPUと標準CGOの材質サンプルは同じ係数を使います。視点座標での光源方向 `(1, 1, 1.5)` はCueMolのOpenGL照明ソースを参照しています。この設定では `toon1` と `toon2` の拡散陰影は同じで、`diff_metal` と `spec_metal` も同じです。CueMolのPOV-Ray定義では、それぞれ異なるfinishで区別します。このプラグインはPOV-Rayの `brilliance`、`phong`、`F_MetalA`、`F_MetalD` を評価しません。CueMol側の描画方式、遮蔽、影、表示変換による差も残ります。

`metallic_chrome`、`metallic_copper`、`stone35`、`wood31`、`wood14scl2` は、対応するPOV-Rayテクスチャ名を参考にした近似です。模様や反射帯は独自実装であり、CueMolのPOV-Rayテクスチャと設定・画素が一致するものではありません。

輪郭線の `thin`、`normal`、`thick` は0.03、0.06、0.15オングストロームで、CueMolのEgLineスタイルと一致します。`outline` は輪郭と鋭い折れ目、`silhouette` は輪郭だけを描きます。複合プリセットの `richardson`、`toon1`、`toon2` は黒いnormal幅の線を選びます。GPUの線抽出と最小1ピクセル制限、rayの円柱状輪郭線は、CueMolの画面上の線描画と異なります。CueMolのrenderer単体の既定幅0.01は、名前付きnormalスタイルとは別の設定です。

RichardsonのkhakiとSteelBlueの色、紙色 `(0.941, 0.925, 0.867)`、線の角度 `55/-35/80`、閾値 `0.92/0.62/0.34`、幅のフェード係数 `10`、インク倍率 `1/0.74/0.38`、暗部の筆圧下限 `0.4`、最小コントラスト `0.15` は参照設定と同じです。指定間隔と線幅は出力画像の `0.5/0.45` ピクセルです。CueMol既定の3倍サンプリングと描画先の最小2ピクセル制限により、実効間隔は出力画像の `2/3` ピクセルになります。

操作中のシェーダーもこの間隔と幅を使い、9サンプルで線をフィルターします。線長と間隙 `50/5`、幅の揺れ `0.45`、長さの揺れ `0.5`、先細り `0.35`、角度の揺れ `5` 度、紙の粒度係数 `0.15` とスケール `3` も参照値に合わせています。乱数は独自で、9サンプル間で緩やかな線幅変化を共有し、各層を平均してから乗算します。Umbreonは各サンプルで線と層の乗算を評価します。遮蔽と深度フォグも省略するため、ピクセル単位では一致しません。ハイライトは紙色のままです。標準CGOとrayは平均被覆率の近似を使い、鉛筆の各層が完全に有効な場合は約0.256です。

## 検証

独立した検証スクリプトで、実PyMOLの全スタイル、失敗後の復元、複数state、セッション再読込、Qtクリック選択、回転、標準・独自形状の半透明表示を確認します。GUIモードでは目視確認用のGPU画像とray画像も保存します。リポジトリのpixiインタープリターでPythonを実行します。

```console
$ uv run --no-project --python .pixi/envs/default/bin/python python \
    tests/check_cuemol_style.py --output .cache/cuemol-headless
$ uv run --no-project --python .pixi/envs/default/bin/python python \
    tests/check_cuemol_style.py --gui --benchmark \
    --output .cache/cuemol-gui
```

ベンチマークは500残基・100個の合成stateを使い、1280×720の描画領域で中品質のリボンを描画します。準備時間、CPU・GPUメッシュ保存量、回転速度、明示的state切替速度、実際の動画描画速度、1 stateの表面準備を報告します。準備後の目標は回転30 FPS、再生15 FPSです。結果はGPUと入力形状に依存します。`.cache` 内の生成画像とレポートはGitの対象外です。

形状と名前は [CueMol ribbon renderer](https://cuemol.github.io/cuemol2_docs/cuemol2/RibbonRenderer/) およびCueMolのスタイル定義を参考にしています。メッシュ生成と描画は独自実装で、実行時にCueMolのソースツリーへ依存しません。
