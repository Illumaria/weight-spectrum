[RU](https://github.com/Illumaria/weight-spectrum/blob/master/docs/problem.md) | [EN](https://github.com/Illumaria/weight-spectrum/blob/master/docs/problem.en.md)

# Linear subspace weight spectrum calculation
1. Let a vector be a string of bits (where a bit is either 0 or 1) of a constant length *N*. We can say that for a given value of *N* there can be no more than ![powN](https://latex.codecogs.com/png.latex?%5Cinline%20%5Clarge%202%5E%7BN%7D) different vectors.
2. Let's introduce the *addition modulo 2* operation (*exclusive OR*, *XOR*), which for two given vectors *a* and *b* produces a vector *a + b* of the same length *N*.
3. Let there be a set ![A](https://latex.codecogs.com/png.latex?%5Cinline%20%5Clarge%20A%20%3D%20%5C%7B%20a_%7Bi%7D%20%7C%20i%20%5Cin%201..K%20%5C%7D) of ![K](https://latex.codecogs.com/png.latex?%5Cinline%20%5Clarge%200%5Cleqslant%20K%5Cleqslant2%5E%7BN%7D) vectors. Let's call it a generating set: by adding up vectors ![a_i](https://latex.codecogs.com/png.latex?%5Cinline%20%5Clarge%20a_%7Bi%7D) of the *A* set we can get ![powK](https://latex.codecogs.com/png.latex?%5Cinline%20%5Clarge%202%5E%7BK%7D) vectors of type ![vectors](https://latex.codecogs.com/png.latex?%5Cinline%20%5Clarge%20%5Csum_%7Bi%3D1%7D%5E%7BK%7D%7B%5Cbeta_%7Bi%7Da_%7Bi%7D%7D), where ![beta_i](https://latex.codecogs.com/png.latex?%5Cinline%20%5Clarge%20%5Cbeta_%7Bi%7D) is either 0 or 1.
4. A natural number between 0 and *N* that is the sum of all ones in a vector is called *vector weight*.

The task is to derive a histogram (spectrum) of the number of different vectors by their weight for the given generating set of vectors and the number *N*.

# Input data format
A text file consisting of a set of strings of the same length, one vector (string) per line, where each vector is a set of numbers 0 or 1 without separators.
For example:
```
001
100
110
```

# Output data format
A text file of lines with a pair of weight/quantity values separated by a tab character, one pair per line, sorted by weight value.
For the example of input data given above:
```
0   1
1   3
2   3
3   1
```

# Optional
1. Use parallel computations where applicable.
2. Evaluate resources, specify limitations and possible methods for further algorithmic and software optimization.
