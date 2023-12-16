import numpy as np

import sys

pairings = {'A': 'U', 'U': 'A', 'G':'C', 'C':'G', 'N': '0', 'Y': '0'}

# compute the back trace of the optimal structure and return the representaiton of that structure
def compute_dot_paren(M, pointers, min_hairpin):
    result = "."*len(M)
    # find end solution to start backtrace
    solution_idx = np.array(M)[:,-1].argmax()

    #run backtrace from solution
    i = solution_idx
    j = len(result) - 1
    k = []
    while i + min_hairpin < j or len(k) > 0:
        # print(f"{i}, {j}")
        if i + min_hairpin >= j:
            # go to next bifurcation
            k.sort()
            last_k = k.pop()
            i = last_k[1]
            j = last_k[0]
        elif pointers[i][j][1] == -1:
            # bifurcation
            k.append((pointers[i][j][0],i))
            i = pointers[i][j][0]+1
        elif pointers[i][j][0] == i + 1 and pointers[i][j][1] == j - 1:
            if M[i][j] > M[i+1][j-1]:
                result = result[:i]+'('+result[i+1:j]+')'+result[j+1:]
            i = i + 1
            j = j - 1
        elif pointers[i][j][1] == j:
            i = i + 1
        elif pointers[i][j][0] == i:
            j = j - 1
        else:
            return result
    return result

# return the dot-parenthese structure corresponding to the sequence
def nussinov_secondary_structure(sequence, min_hairpin_length=3):
    M = [[0 for _ in range(len(sequence))] for _ in range(len(sequence))]
    pointers = [[(0,0) for _ in range(len(sequence))] for _ in range(len(sequence))]
    for j in range(len(sequence)):
        for i in range(len(sequence))[::-1]:
            if i + min_hairpin_length >= j:
                M[i][j] = 0
            else:
                dir = (i+1, j-1)
                match = False
                try:
                    match = (sequence[i] == pairings[sequence[j]])
                except:
                    match = False
                if match:
                    M[i][j] = M[i+1][j-1] + 1
                else:
                    M[i][j] = M[i+1][j-1]
                if M[i+1][j] > M[i][j]:
                    M[i][j] = M[i+1][j]
                    dir = (i+1, j)
                if M[i][j-1] > M[i][j]:
                    M[i][j] = M[i][j-1]
                    dir = (i, j-1)
                max = M[i][j]
                max_k = -1
                for k in range(i+1 + min_hairpin_length, j-min_hairpin_length):
                    sum = M[i][k] + M[k+1][j]
                    if sum > max:
                        max = sum
                        max_k = k
                if max_k > -1:
                    M[i][j] = max
                    dir = (max_k, -1)
                pointers[i][j] = dir
    return compute_dot_paren(M, pointers, min_hairpin_length)


def compute_for_file(file_name, file_dir='.', sequences_to_read=None, lines_to_skip=0):
    # file_dir = 'data'
    # file_name = 'microgreen_id_rna'
    # file_dir = 'test'
    # file_name = "test_seq"
    structure_file_name = f'{file_name}_structure'
    if sequences_to_read != None:
        if len(sequences_to_read) == 1:
            structure_file_name += f'_{sequences_to_read[0]}'
        else:
            structure_file_name += f'_{len(sequences_to_read)}_{sequences_to_read[0]}'
    structure_file = open(f'output/{structure_file_name}.fasta', 'w')
    min_len = 3

    with open(f'{file_dir}/{file_name}.fasta', 'r') as f:
        line = f.readline()
        line_count = 0
        read_seq = False
        while line:
            if ((line[0] == '>' or len(line) == 0) and (sequences_to_read==None or line[1:-1] in sequences_to_read) and line_count>=lines_to_skip):
                structure_file.write(line)
                read_seq = True
            elif read_seq:
                structure_file.write(line)
                structure_file.write(nussinov_secondary_structure(line[:-1], min_len)+'\n')
                read_seq = False
            line = f.readline()
            line_count += 1
    structure_file.close()

compute_for_file('microgreen_id_rna', 'data', lines_to_skip=56)