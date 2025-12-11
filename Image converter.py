from PIL import Image


def resize_image(input_path, output_path, size):
    original_image = Image.open(input_path)
    old_width, old_height = original_image.size
    target_width, target_height = size

    # Рассчитываем пропорции
    aspect_ratio = old_width / old_height
    target_ratio = target_width / target_height

    if target_ratio > aspect_ratio:
        # Если целевая картинка шире, меняем размер по ширине
        new_width = target_width
        new_height = int(target_width / aspect_ratio)
    else:
        # Иначе меняем по высоте
        new_height = target_height
        new_width = int(target_height * aspect_ratio)

    # Ресайз с высоким качеством
    resized_img = original_image.resize((new_width, new_height), Image.Resampling.LANCZOS)

    # Обрезка по центру (Crop)
    left = (new_width - target_width) / 2
    top = (new_height - target_height) / 2
    right = (new_width + target_width) / 2
    bottom = (new_height + target_height) / 2

    final_img = resized_img.crop((left, top, right, bottom))
    final_img.save(output_path)
    print(f"Готово! Картинка сохранена как {output_path}")


# Использование
resize_image("chemistry.jpg", "chemistry_640x360.jpg", (640, 360))