# Матрица миграции

| Current path | Target path | Type | Action | Reason |
| --- | --- | --- | --- | --- |
| `laba1`–`laba8` | `labs/master/lab01`–`lab08` | master labs | move | Сохранить учебную последовательность курса и разделить работы по темам. |
| `обратные задачи бакалавриат/практика N` | `education/bachelor/practice/practiceNN` | bachelor practice | move | Сохранить материалы бакалавриата с едиными каталогами. |
| `обратные задачи бакалавриат/лекция*.pdf` | `education/bachelor/lectures` | bachelor lectures | move/rename | Отделить лекционный материал от исходного каталога. |
| `inverse_problems_labs.pdf` | `references/methodics/inverse_problems_labs.pdf` | methodics | move | Методичка является источником структуры текущих лабораторных. |
| `lab5*.py` | `labs/master/lab05/variants` | experiment variants | move | Сохранить альтернативные реализации рядом с лабораторией. |
| `laba6*.py`, `laba7.py` | `labs/master/lab06/variants`, `lab07/variants` | FEM variants | move | Отделить варианты МКЭ от общего кода. |
| root `*.png`, `*.gif`, `*.npy`, `result/*` | `results/labXX` | generated results | move | Не смешивать результаты с исходным кодом. |
| `temp/*` | `archive/legacy/temp` | historical material | move | Сохранить содержательные черновики и результаты. |
| empty `practika.ipynb` | `archive/legacy/invalid-notebooks` | invalid historical notebook | move | Сохранить пустой исторический файл отдельно от проверяемых notebooks. |
| `gif_to_png.py` | `scripts/gif_to_png.py` | utility | move/update | Утилита не относится к конкретной лаборатории. |
| `jitfailure-*`, `__pycache__`, `*.pyc` | — | generated cache | delete | Воспроизводимые служебные файлы не являются исходными материалами. |
| identical copies in `lab04`/`lab08` | — | duplicate | delete | Удалены только подтверждённые хешем дубликаты. |

Git-индекс в среде разработки доступен только для чтения, поэтому физические
перемещения были выполнены обычным перемещением файлов. После получения
доступа к индексу Git сможет распознать их как renames.
