"""Portable CueMol-inspired profiles; no CueMol runtime is required."""

from dataclasses import dataclass, replace


@dataclass(frozen=True)
class Profile:
    representation: str = "auto"
    section: str = "rectangle"
    material: str = "default"
    edges: str = "none"
    edge_width: float = 0.06
    back: bool = False


MATERIALS = (
    "default",
    "shadow",
    "nolighting",
    "matte",
    "toon1",
    "toon2",
    "diff_metal",
    "spec_metal",
    "metallic_chrome",
    "metallic_copper",
    "stone35",
    "wood31",
    "wood14scl2",
)
REPRESENTATIONS = (
    "auto",
    "ribbon",
    "cartoon",
    "tube",
    "nucleic",
    "ballstick",
    "sticks",
    "cpk",
    "surface",
)
QUALITIES = {"low": (4, 8), "medium": (8, 16), "high": (12, 24)}
COLORS = ("cuemol", "keep", "chain", "ss", "rainbow", "element")

# Optional WoodyHSCPaint, retained as color=ss.
SECONDARY_COLORS = {"H": (1.0, 127 / 255, 127 / 255), "S": (127 / 255, 1.0, 127 / 255)}
COIL_COLOR = (1.0, 1.0, 191 / 255)

# CueMol GUI's createDefPaintColoring, referenced by DefaultHSCPaint and
# DefaultCPKColoring through $molcol for newly loaded molecules.
CUEMOL_SECONDARY_COLORS = {
    "H": (240 / 255, 230 / 255, 140 / 255),
    "S": (70 / 255, 130 / 255, 180 / 255),
}
CUEMOL_OTHER_COLOR = (1.0, 250 / 255, 240 / 255)
CUEMOL_NUCLEIC_COLOR = (1.0, 1.0, 0.0)
CUEMOL_ELEMENTS = {
    "N": (0.0, 0.0, 1.0),
    "O": (1.0, 0.0, 0.0),
    "H": (0.0, 1.0, 1.0),
    "S": (0.0, 1.0, 0.0),
    "P": (1.0, 1.0, 0.0),
}
PROFILES = {name: Profile(material=name) for name in MATERIALS}
PROFILES.update(
    {
        "richardson": Profile("ribbon", "fancy", "toon1", "edges", back=True),
        "ribbon": Profile("ribbon"),
        "round_ribbon": Profile("ribbon", "ellipse"),
        "fancy_ribbon": Profile("ribbon", "fancy", back=True),
        "cartoon": Profile("cartoon"),
        "round_cartoon": Profile("cartoon", "ellipse"),
        "tube": Profile("tube"),
        "nucleic": Profile("nucleic"),
        "ballstick": Profile("ballstick"),
        "cpk": Profile("cpk"),
        "surface": Profile("surface"),
        "outline": Profile(edges="edges"),
        "silhouette": Profile(edges="silhouette"),
    }
)
for _name in ("toon1", "toon2"):
    PROFILES[_name] = replace(PROFILES[_name], edges="edges", back=True)


def resolve(style, representation="auto", edge="auto", edge_width=None):
    if style not in PROFILES:
        raise ValueError(f"Unknown style {style!r}; use 'cuemol_style list'.")
    if representation not in REPRESENTATIONS:
        raise ValueError(f"representation must be one of {REPRESENTATIONS}")
    p = PROFILES[style]
    if representation != "auto":
        p = replace(p, representation=representation)
    if edge in ("thin", "normal", "thick"):
        p = replace(
            p,
            edges="edges",
            edge_width={"thin": 0.03, "normal": 0.06, "thick": 0.15}[edge],
        )
    elif edge in ("none", "edges", "silhouette"):
        p = replace(p, edges=edge)
    elif edge != "auto":
        raise ValueError(
            "edge must be auto, none, edges, silhouette, thin, normal, or thick"
        )
    if edge_width is not None:
        import math

        width = float(edge_width)
        if not math.isfinite(width) or width <= 0:
            raise ValueError("edge_width must be finite and positive (angstroms)")
        p = replace(p, edge_width=width)
    return p
