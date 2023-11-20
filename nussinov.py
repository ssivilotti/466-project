import numpy as np

pairings = {'A': 'U', 'U': 'A', 'G':'C', 'C':'G'}

# compute the back trace of the optimal structure and return the representaiton of that structure
def compute_dot_paren(M, pointers):
    result = "."*len(M)
    # implement
    return result

# return the dot-parenthese structure corresponding to the sequence
def secondary_structure(sequence):
    M = [[0 for _ in range(len(sequence))] for _ in range(len(sequence))]
    pointers = [[(0,0) for _ in range(len(sequence))] for _ in range(len(sequence))]
    for j in range(len(sequence)):
        for i in range(sequence)[::-1]:
            if i >= j:
                M[i][j] = 0
            else:
                if sequence[i] == pairings[sequence[j]]:
                    M[i][j] = M[i+1][j-1] + 1
                else:
                    M[i][j] = M[i+1][j-1] + 1
                pointers[i][j] = (i+1, j-1)
                if M[i+1][j] > M[i][j]:
                    M[i][j] = M[i+1][j]
                    pointers = (i+1, j)
                if M[i][j-1] > M[i][j]:
                    M[i][j] = M[i][j-1]
                    pointers = (i, j-1)
                max = M[i][j]
                max_k = -1
                for k in range(i+1, j):
                    sum = M[i][k] + M[k+1][j]
                    if sum > max:
                        max = sum
                        max_k = k
                if max_k > -1:
                    M[i][j] = max
                    pointers = (k)
    return compute_dot_paren(M,pointers)

file_name = 'microgreen_id_rna'
structure_file = open(f'output/{file_name}_structure', 'w')

with open(f'{file_name}.fasta', 'r') as f:
    line = f.readline()
    while line:
        if (line[0] == '>'):
            structure_file.write(line)
        else:
            structure_file.write(secondary_structure(line))
        line = f.readline()
structure_file.close()