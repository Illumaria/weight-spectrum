#! /usr/bin/env python3
# -*- coding: UTF-8 -*-

import argparse
import sys
from math import log2
from multiprocessing import Process, Manager, cpu_count
from pathlib import Path


def read(path):
    """Read a file

    Args:
        path (str): Path to the file to read from.

    Returns:
        :rtype: (:obj:`list` of :obj:`int`, int, int): The list of
            vectors read from file, their number and their length.

    """
    with open(path, "r", encoding='utf-8') as fin:
        vectors = [x.rstrip() for x in fin.readlines()]
        vector_len = len(vectors[0])
        vector_num = len(vectors)
        vectors = list(map(lambda x: int(x, 2), vectors))
        return (vectors, vector_num, vector_len)


def get_zeros_pos(vector):
    """Get positions of zeros in a vector

    Args:
        vector (:obj:`list` of :obj:`int`): The input vector.

    Returns:
        :rtype: (:obj:`list` of :obj:`int`, int): The list of
            zeros positions in the vector and the length of
            the vector with all zeros removed.

    """
    zeros_pos = []
    vector_str = bin(vector)[2:]
    vector_len = len(vector_str)
    zeros_pos = [i for i, j in enumerate(vector_str[::-1]) if j == '0']
    new_vector_len = vector_len - len(zeros_pos)
    return (zeros_pos, new_vector_len)


def delete(vector, zeros_pos):
    """Remove zeros at given indexes from a single vector

    Args:
        vector (int): The ingeger value corresponding to the vector.
        zeros_pos (:obj:`list` of :obj:`int`): The list of
            zeros indexes in the vector, sorted in descending order.

    Returns:
        int: The ingeger value corresponding to the vector
            arfer removing zeros at given positions
            from vector's binary representation.

    Example:
        >>> print(delete(37, [3, 1]))
        11

    """
    bin_vector = bin(vector)[2:]
    bin_len = len(bin_vector)
    for zero_pos in zeros_pos[::-1]:
        if zero_pos < bin_len:
            bin_vector = (bin_vector[:bin_len-zero_pos-1]
                          + bin_vector[bin_len-zero_pos:])
    vector = int(bin_vector, 2)
    return vector


def delete_zeros(vectors, vector_len):
    """Delete common zeros from vectors

    Args:
        vectors (:obj:`list` of :obj:`int`): The list of
            vectors.
        vector_len (int): The length of vectors.

    Returns:
        :rtype: (:obj:`list` of :obj:`int`, int): The list of
            vectors with all common zeros removed and the new length
            of vectors.

    Example:
        >>> print(delete_zeros([5, 4], 3))
        ([3, 2], 2)

    """
    vec_sum = 0
    for vec in vectors:
        vec_sum |= vec
    zeros_pos, new_vector_len = get_zeros_pos(vec_sum)
    vectors = list(map(lambda x: delete(x, zeros_pos), vectors))
    return (vectors, new_vector_len)


def get_basis(vectors, vector_num, vector_len):
    """Get vectors basis

    Args:
        vectors (:obj:`list` of :obj:`int`): The list of
            vectors.
        vector_num (int): The number of vectors in the list.
        vector_len (int): The length of vectors in the list.

    Returns:
        :rtype: (:obj:`list` of :obj:`int`, int): The list of
            basis vectors and the rank of the basis.

    """
    # Initial rank equals to the current full rank
    rank = min(vector_len, vector_num)

    for r in range(rank):
        vectors = sorted(vectors, reverse=True)
        index = len(bin(vectors[r])[2:]) - 1
        for i in range(vector_num):
            if (vectors[i] & 1 << index) and (i != r):
                vectors[i] ^= vectors[r]

    basis = [vectors[i] for i in range(rank) if vectors[i]]
    # The final rank equals to the number of rows in basis matrix
    rank = len(basis)
    return (basis, rank)


def partition(start, end, cores):
    """Split a range into (exactly or almost) parts

    Args:
        start (int): The first index of the range.
        end (int): The last index of the range.
        cores (int): The number of parts to split into.

    Returns:
        :obj:`list` of :obj:`list` of :obj:`int`: The list of
            lists, where each inner list contains starting and
            ending index of a single part.

    Examples:
        >>> print(0, 100, 3)
        [[0, 33], [34, 67], [68, 100]]
        >>> print(10, 49, 4)
        [[10, 19], [20, 29], [30, 39], [40, 49]]

    """
    dn = round((end - start + 1) / cores)
    parts = []
    parts += [[start, start + dn - 1]]
    parts += [[start + dn*(i-1), start + i*dn - 1] for i in range(2, cores)]
    parts += [[start + (cores-1)*dn, end]]
    return parts


def count_ones(int_num):
    """Count ones in a binary representation of an integer number

    Args:
        int_num (int): The integer number.

    Returns:
        int: The number of ones in the binary representation
            of the given number.

    Examples:
        >>> print(count_ones(5))
        2
        >>> print(count_ones(11))
        3

    """
    return bin(int_num).count('1')


def gray_code(index):
    """Get the Gray code equivalent of the given integer number

    Args:
        index (int): The integer number.

    Returns:
        int: The integer number whose binary representation
            corresponds to Gray code of the index.

    Examples:
        >>> print(gray_code(5))
        7
        >>> print(gray_code(11))
        14

    """
    return index ^ (index // 2)


def get_spectrum(basis, vector_len, bounds, total_spectrum):
    """Get the spectrum of a basis

    Args:
        basis (:obj:`list` of :obj:`int`): The list of
            basis vectors.
        vector_len (int): The length of basis vectors.
        bounds (:obj:`list` of :obj:`int`): The list of
            starting and ending positions of the range
            to calculate spectrum for.
        total_spectrum (:obj:`list` of :obj:`list` of :obj:`int`):
            The list of lists where a single inner list
            contains spectrum for a single specific range.

    """
    spectrum = [0] * (vector_len + 1)

    # Calculate the first vector
    # in the given range and its weight
    current_vector = 0
    if bounds[0] != 0:
        i = 0
        gray_encode = gray_code(bounds[0])
        while gray_encode:
            if gray_encode % 2 == 1:
                current_vector ^= basis[i]
            i += 1
            gray_encode //= 2

    spectrum[count_ones(current_vector)] += 1

    # Calculate weights of other vectors
    for i in range(bounds[0], bounds[1]):
        bit_change_pos = int(log2((-1-i) & (1+i)))
        current_vector ^= basis[bit_change_pos]
        weight = count_ones(current_vector)
        spectrum[weight] += 1

    # Append weights from each process
    # to the main list
    total_spectrum.append(spectrum)


def process(basis, rank, vector_len_wz, vector_len, vector_num, cores):
    """The main process of the script for spectrum calculation

    Args:
        basis (:obj:`list` of :obj:`int`): The list of basis vectors.
        rank (int): The rank of the basis.
        vector_len_wz (int): The length of basis vectors after removing
            redundant zeros.
        vector_len (int): The starting basis vectors length.
        vector_num (int): The number of vectors in the basis.
        cores (int): The number of parallel processes to run.

    Returns:
        spectrum (:obj:`list` of :obj:`int`): The list of
            weights of the basis vectors.

    """
    spectrum = []
    if rank == vector_len_wz:
        spectrum = [1]
        for i in range(1, rank+1):
            spectrum.append(int(spectrum[i-1] * (rank-i+1) / i))
    else:
        if cores > 1:
            print("Using", cores, "cores for parallel computing.")
        else:
            print("Using 1 core for parallel computing.")
        parts = partition(0, 2**len(basis) - 1, cores)
        with Manager() as manager:
            total_spectrum = manager.list()

            processes = [
                Process(
                    target=get_spectrum,
                    args=(basis, vector_len, part, total_spectrum)
                    ) for part in parts
                ]

            for p in processes:
                p.start()
            for p in processes:
                p.join()
            spectrum = [sum(x) for x in zip(*total_spectrum)]
    spectrum = [int(weight * 2**(vector_num - rank)) for weight in spectrum]
    return spectrum


def write(spectrum, path):
    """Write spectrum to a file

    Args:
        spectrum (:obj:`list` of :obj:`int`): The list of
            weights of the basis vectors.
        path (str): Path to the file to write to.

    Returns:
        :rtype: (:obj:`list` of :obj:`int`, int, int): The list of
            vectors read from file, their number and their length.

    """
    with open(path, "w", encoding='utf-8') as fout:
        for i, elem in enumerate(spectrum[:-1]):
            fout.write("{}\t{}\n".format(i, elem))
        else:
            fout.write("{}\t{}".format(i+1, spectrum[-1]))


def arg_parser():
    """The function to fill command line arguments and parse them"""
    parser = argparse.ArgumentParser(
        description="Calculate linear subspace weight spectrum."
        )

    parser.add_argument("input", help="a path to the input text file")
    parser.add_argument("-o", "--output", default="./output.txt",
                        help="a path to the output"
                        + " text file; if not specified,"
                        " \"./output.txt\" will be used")
    parser.add_argument("-j", "--parallel", type=int,
                        choices=[i+1 for i in range(cpu_count())],
                        default=cpu_count(),
                        help="number of cores to use for parallel"
                        + " computing; if not specified, maximum"
                        + " number of cores will be used")

    args = parser.parse_args()
    return args


def main(input_file, output_file, cores):
    """The main function of the script

    Args:
        input_file (str): The path to the input file.
        output_file (str): The path to the output file.
        cores (int): The number of parallel processes to run.

    """
    # Read vectors from the file and get their number and length
    vectors, vector_num, vector_len = read(input_file)

    # Delete zeros and get the new
    # vector length (without zeros)
    vectors, vector_len_wz = delete_zeros(vectors, vector_len)

    # Find basis and rank of vectors matrix
    basis, rank = get_basis(vectors, vector_num, vector_len_wz)

    # Calculate spectrum
    spectrum = process(basis, rank, vector_len_wz,
                       vector_len, vector_num, cores)

    # Write spectrum to the file
    print("Writing output to {}...".format(output_file))
    write(spectrum, output_file)
    print("Done.")


if __name__ == "__main__":
    args = arg_parser()
    main(Path(args.input), Path(args.output), args.parallel)
