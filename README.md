[RU](https://github.com/Illumaria/weight-spectrum/blob/master/README.md) | [EN](https://github.com/Illumaria/weight-spectrum/blob/master/README.en.md)

# Весовой спектр линейного подпространства
Данный скрипт предназначен для вычисления весового спектра линейного подпространства.

## [Постановка задачи](https://github.com/Illumaria/weight-spectrum/blob/master/docs/problem.md)

## [Решение](https://github.com/Illumaria/weight-spectrum/blob/master/docs/solution.md)

# Использование
Для использования скрипта откройте терминал и выполните следующие команды:
```
$ git clone https://github.com/Illumaria/weight-spectrum.git
$ cd weight-spectrum/
$ python3 ./weight-spectrum.py
```
Для работы скрипта также необходимо указать как минимум путь к файлу с входными данными.
Полный список параметров запуска:
* `<input-file>`: путь к файлу с исходными данными;
* `<output-file>`: (необязательно) путь к файлу, куда будут сохранены выходные данные; при отсутствии значения выходные данные будут сохранены в файл `output.txt` в текущем каталоге;
* `-j <processes>`: (необязательно) количество процессов (ядер процессора), которое будет задействовано при вычислениях; при отсутствии значения будет использовано максимально возможное количество ядер, полученное автоматически.

# Лицензия
Скрипт может использоваться на условиях лицензии [MIT license](https://opensource.org/licenses/MIT).
