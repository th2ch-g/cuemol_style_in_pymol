"""Measure and assemble camera-aligned audit images without registering them."""

import argparse
import json
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFilter


def rgb(path, background):
    source = Image.open(path).convert("RGBA")
    canvas = Image.new("RGBA", source.size, (*background, 255))
    return Image.alpha_composite(canvas, source).convert("RGB")


def compare(directory):
    rows, tiles = [], []
    for path in sorted(directory.glob("*-pymol.png")):
        name = path.name.removesuffix("-pymol.png")
        reference = directory / f"{name}-cuemol.png"
        if not reference.is_file():
            continue
        background = json.loads(path.with_name(f"{name}.json").read_text())[
            "background"
        ]
        background = np.array(background) * 255
        color = tuple(np.rint(background).astype(int))
        actual, expected = rgb(path, color), rgb(reference, color)
        if actual.size != expected.size:
            raise ValueError(f"Image dimensions differ for {name}")
        a, b = np.asarray(actual, dtype=float), np.asarray(expected, dtype=float)
        ma, mb = (
            np.abs(a - background).max(axis=2) > 10,
            np.abs(b - background).max(axis=2) > 10,
        )
        union = ma | mb
        blur_a = np.asarray(actual.filter(ImageFilter.GaussianBlur(2)), dtype=float)
        blur_b = np.asarray(expected.filter(ImageFilter.GaussianBlur(2)), dtype=float)
        rows.append(
            {
                "profile": name,
                "foreground_iou": float((ma & mb).sum() / max(1, union.sum())),
                "foreground_mae_255": float(np.abs(a - b)[union].mean())
                if union.any()
                else 0,
                "blur2_foreground_mae_255": float(np.abs(blur_a - blur_b)[union].mean())
                if union.any()
                else 0,
                "full_image_mae_255": float(np.abs(a - b).mean()),
            }
        )
        tile = Image.new("RGB", (960, 274), "white")
        draw = ImageDraw.Draw(tile)
        difference = Image.fromarray(np.uint8(np.clip(np.abs(a - b) * 3, 0, 255)))
        for col, (label, value) in enumerate(
            (
                ("PyMOL", actual),
                ("CueMol", expected),
                ("absolute difference x3", difference),
            )
        ):
            value.thumbnail((320, 240))
            tile.paste(value, (320 * col, 24))
            draw.text((320 * col + 4, 4), f"{name}: {label}", fill="black")
        tiles.append(tile)
    if not rows:
        raise ValueError("No paired audit images found")
    for start in range(0, len(tiles), 6):
        page = Image.new("RGB", (960, 274 * len(tiles[start : start + 6])), "white")
        for offset, tile in enumerate(tiles[start : start + 6]):
            page.paste(tile, (0, offset * 274))
        page.save(directory / f"comparison-{start // 6 + 1}.png")
    report = {
        "method": "No image registration; foreground differs from the manifest background by >10 in any channel; MAE uses 0..255 channels; blur uses sigma=2 pixels. Foreground IoU includes shading and is not a pure geometry score. Richardson uses an opaque paper background in both engines; other cases use white.",
        "cases": rows,
    }
    (directory / "comparison.json").write_text(json.dumps(report, indent=2) + "\n")
    for row in rows:
        print(
            f"{row['profile']:20s} IoU {row['foreground_iou']:.4f}  MAE {row['foreground_mae_255']:.3f}  blur {row['blur2_foreground_mae_255']:.3f}"
        )
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    compare(parser.parse_args().directory)
