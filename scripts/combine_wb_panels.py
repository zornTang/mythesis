from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "figures"

INPUTS = [
    ("(a) HL-60", FIG_DIR / "HL-60.png"),
    ("(b) 293T", FIG_DIR / "293T.png"),
]
OUTPUT = FIG_DIR / "WB_combined.png"


def load_font(size: int) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    candidates = [
        "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Supplemental/Helvetica.ttc",
        "/Library/Fonts/Arial.ttf",
    ]
    for candidate in candidates:
        path = Path(candidate)
        if path.exists():
            return ImageFont.truetype(str(path), size=size)
    return ImageFont.load_default()


def main() -> None:
    images = [(label, Image.open(path).convert("RGB")) for label, path in INPUTS]

    panel_width = max(image.width for _, image in images)
    panel_height = max(image.height for _, image in images)
    margin = 48
    gap = 36
    label_height = 52

    canvas_width = margin * 2 + panel_width * len(images) + gap * (len(images) - 1)
    canvas_height = margin * 2 + label_height + panel_height

    canvas = Image.new("RGB", (canvas_width, canvas_height), "white")
    draw = ImageDraw.Draw(canvas)
    font = load_font(28)

    x = margin
    for label, image in images:
        y = margin + label_height
        offset_x = x + (panel_width - image.width) // 2
        offset_y = y + (panel_height - image.height) // 2
        canvas.paste(image, (offset_x, offset_y))
        draw.text((x, margin), label, fill="black", font=font)
        x += panel_width + gap

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(OUTPUT, quality=100)
    print(f"Saved {OUTPUT}")


if __name__ == "__main__":
    main()
