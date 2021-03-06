[RU](https://github.com/Illumaria/weight-spectrum/blob/master/docs/problem.md) | [EN](https://github.com/Illumaria/weight-spectrum/blob/master/docs/problem.en.md)

# Вычисление весового спектра линейного подпространства

1. Назовём вектором строку битов (значения 0 или 1) фиксированной длины *N*. Для заданного значения *N* возможно максимум ![powN](https://latex.codecogs.com/png.latex?%5Cinline%20%5Clarge%202%5E%7BN%7D) различных векторов.
2. Введём операцию сложения по модулю 2 векторов (операция *xor*), которая для двух векторов *a* и *b* получает вектор *a + b* той же длины *N*.
3. Пусть задано множество ![A](https://latex.codecogs.com/png.latex?%5Cinline%20%5Clarge%20A%20%3D%20%5C%7B%20a_%7Bi%7D%20%7C%20i%20%5Cin%201..K%20%5C%7D) из ![K](https://latex.codecogs.com/png.latex?%5Cinline%20%5Clarge%200%5Cleqslant%20K%5Cleqslant2%5E%7BN%7D) векторов. Назовём его порождающим: при помощи сложения векторов ![a_i](https://latex.codecogs.com/png.latex?%5Cinline%20%5Clarge%20a_%7Bi%7D) множества *A* можно получить ![powK](https://latex.codecogs.com/png.latex?%5Cinline%20%5Clarge%202%5E%7BK%7D) векторов вида ![vectors](https://latex.codecogs.com/png.latex?%5Cinline%20%5Clarge%20%5Csum_%7Bi%3D1%7D%5E%7BK%7D%7B%5Cbeta_%7Bi%7Da_%7Bi%7D%7D), где коэффициент ![beta_i](https://latex.codecogs.com/png.latex?%5Cinline%20%5Clarge%20%5Cbeta_%7Bi%7D) равен либо 0, либо 1.
4. Весом вектора назовём количество единичных (ненулевых) битов в векторе. Таким образом, вес — это натуральное число от 0 до *N*.

Необходимо для заданных порождающего множества векторов и числа *N* построить гистограмму (спектр) количества различных векторов по их весу.

# Формат входных данных
Текстовый файл из набора строк одинаковой длины по одному вектору в строке (символы 0 или 1 без разделителей).
Например:
```
001
100
110
```

# Формат выходных данных
Текстовый файл строк с парой значений вес/количество, разделённых символом табуляции, по одной паре на строку, сортированный по числовому значению веса.
Для примера входных данных выше:
```
0   1
1   3
2   3
3   1
```

# Дополнительно
1. По возможности использовать параллельные вычисления.
2. Оценить ресурсы, указать ограничения реализации, возможные методы дальнейшей алгоритмической и программной оптимизации.

[На главную](https://github.com/Illumaria/weight-spectrum/blob/master/README.md) | [Решение](https://github.com/Illumaria/weight-spectrum/blob/master/docs/solution.md)
