"""
gif_to_png.py — конвертация GIF в PNG

Использование:
  python scripts/gif_to_png.py results/lab01/solution.gif
  python scripts/gif_to_png.py results/lab01/solution.gif --all
  python scripts/gif_to_png.py results/lab01/*.gif --all -o results/frames

Зависимости:
  pip install Pillow
"""

import argparse
import sys
from pathlib import Path
from PIL import Image


def gif_to_png(gif_path: Path, output_dir: Path, all_frames: bool) -> list[Path]:
    """Конвертирует один GIF файл. Возвращает список созданных PNG."""
    output_dir.mkdir(parents=True, exist_ok=True)
    created = []

    with Image.open(gif_path) as img:
        n_frames = getattr(img, "n_frames", 1)

        if not all_frames or n_frames == 1:
            # Только первый кадр
            frame = img.convert("RGBA")
            out_path = output_dir / (gif_path.stem + ".png")
            frame.save(out_path, "PNG")
            created.append(out_path)
            print(f"  ✓ {out_path}")
        else:
            # Все кадры
            digits = len(str(n_frames - 1))
            for i in range(n_frames):
                img.seek(i)
                frame = img.convert("RGBA")
                out_name = f"{gif_path.stem}_frame{str(i).zfill(digits)}.png"
                out_path = output_dir / out_name
                frame.save(out_path, "PNG")
                created.append(out_path)
                print(f"  ✓ {out_path}  (кадр {i + 1}/{n_frames})")

    return created


def main():
    parser = argparse.ArgumentParser(
        description="Конвертация GIF → PNG (через Pillow)"
    )
    parser.add_argument("input", nargs="+", help="GIF файл(ы) или glob-паттерн")
    parser.add_argument(
        "--all", "-a",
        action="store_true",
        help="Сохранить все кадры (по умолчанию — только первый)"
    )
    parser.add_argument(
        "--output", "-o",
        default=None,
        help="Папка для PNG (по умолчанию — рядом с GIF)"
    )
    args = parser.parse_args()

    # Собираем список файлов
    gif_files = []
    for pattern in args.input:
        paths = list(Path(".").glob(pattern)) if "*" in pattern else [Path(pattern)]
        gif_files.extend(paths)

    gif_files = [p for p in gif_files if p.suffix.lower() == ".gif"]

    if not gif_files:
        print("Ошибка: GIF файлы не найдены.", file=sys.stderr)
        sys.exit(1)

    total_created = []

    for gif_path in gif_files:
        if not gif_path.exists():
            print(f"Файл не найден: {gif_path}", file=sys.stderr)
            continue

        output_dir = Path(args.output) if args.output else gif_path.parent
        print(f"\n{gif_path}  →  {output_dir}/")

        try:
            created = gif_to_png(gif_path, output_dir, all_frames=args.all)
            total_created.extend(created)
        except Exception as e:
            print(f"  Ошибка: {e}", file=sys.stderr)

    print(f"\nГотово: создано {len(total_created)} файл(ов).")


if __name__ == "__main__":
    main()
