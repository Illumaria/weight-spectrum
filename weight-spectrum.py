#!/usr/bin/env python3

import argparse
import sys
from math import log2
from multiprocessing import Process, Manager, cpu_count


# Read the file as a list of integers (vectors)
def read(path):
    with open(path, "r", encoding='utf-8') as fin:
        vectors = [x.rstrip() for x in fin.readlines()]
        vector_len = len(vectors[0])
        vector_num = len(vectors)
        vectors = list(map(lambda x: int(x, 2), vectors))
        return (vectors, vector_num, vector_len)


# Get positions of zeros in a single vector
def get_zeros_pos(vector):
    zeros_pos = []
    vector_str = bin(vector)[2:]
    vector_len = len(vector_str)
    zeros_pos = [i for i, j in enumerate(vector_str[::-1]) if j == '0']
    new_vector_len = vector_len - len(zeros_pos)
    return (zeros_pos, new_vector_len)


# Delete zeros at given indexes from a single vector
def delete(vector, zeros_pos):
    bin_vector = bin(vector)[2:]
    bin_len = len(bin_vector)
    for zero_pos in zeros_pos[::-1]:
        if zero_pos < bin_len:
            bin_vector = (bin_vector[:bin_len-zero_pos-1]
                         + bin_vector[bin_len-zero_pos:])
    vector = int(bin_vector, 2)
    return vector


# Delete common zeros from vectors
def delete_zeros(vectors, vector_len):
    vec_sum = 0
    for vec in vectors:
        vec_sum |= vec
    zeros_pos, new_vector_len = get_zeros_pos(vec_sum)
    vectors = list(map(lambda x: delete(x, zeros_pos), vectors))
    return (vectors, new_vector_len)


# Get basis of vectors
def get_basis(vectors, vector_num, vector_len):
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


# Split a range into (exactly or almost) equal parts
def partition(start, end, cores):
    dn = round((end - start + 1) / cores)
    parts = []
    parts += [[start, start + dn - 1]]
    parts += [[dn * (i - 1), i * dn - 1] for i in range(2, cores)]
    parts += [[(cores - 1) * dn, end]]
    return parts


# Count ones in a binary representation
# of an integer number
def count_ones(int_num):
    # The code provided below does not invoke
    # any type convertions and only performs
    # as many cycles as there are '1' digits
    # in a binary representation of a number

    # count = 0
    # while(int_num != 0):
    #     int_num &= int_num - 1
    #     count += 1
    # return count

    # However, this solution turns out to be faster
    return bin(int_num).count('1')


# Get the Gray code of a given index
def gray_code(index):
    return index ^ (index // 2)


# Get the spectrum of a basis
def get_spectrum(basis, vector_len, bounds, total_spectrum):
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


# The main process of the script
# for calculating spectrum
def process(basis, rank, vector_len_wz, vector_len, vector_num, cores):
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


# Write spectrum to the file
def write(spectrum, path):
    with open(path, "w", encoding='utf-8') as fout:
        for i, elem in enumerate(spectrum[:-1]):
            fout.write("{}\t{}\n".format(i, elem))
        else:
            fout.write("{}\t{}".format(i+1, spectrum[-1]))

            
def arg_parser():
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
    main(args.input, args.output, args.parallel)


# in_24_32.txt: 1 loop, best of 3: 20.2 s per loop (-18.5% vs Git/Habr)
# in_31_32.txt: 1233.2595570087433 seconds         (-15.7% vs Git/Habr)
