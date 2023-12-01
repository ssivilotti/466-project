import numpy as np

pairings = {'A': 'U', 'U': 'A', 'G':'C', 'C':'G'}

# compute the back trace of the optimal structure and return the representaiton of that structure
def compute_dot_paren(M, pointers):
    result = "."*len(M)
    # find end solution to start backtrace
    solution_idx = np.array(M)[:,-1].argmax()

    #run backtrace from solution
    i = solution_idx
    j = len(result) - 1
    k = []
    while i < j or len(k) > 0:
        # print(f"{i}, {j}")
        if i >= j:
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
def secondary_structure(sequence, min_hairpin_length=4):
    M = [[0 for _ in range(len(sequence))] for _ in range(len(sequence))]
    pointers = [[(0,0) for _ in range(len(sequence))] for _ in range(len(sequence))]
    for j in range(len(sequence)):
        for i in range(len(sequence))[::-1]:
            if i + min_hairpin_length >= j:
                M[i][j] = 0
            else:
                dir = (i+1, j-1)
                if sequence[i] == pairings[sequence[j]]:
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
                for k in range(i+1 + min_hairpin_length, j):
                    sum = M[i][k] + M[k+1][j]
                    if sum > max:
                        max = sum
                        max_k = k
                if max_k > -1:
                    M[i][j] = max
                    dir = (max_k, -1)
                pointers[i][j] = dir
    return compute_dot_paren(M,pointers)

# file_dir = 'data'
# file_name = 'microgreen_id_rna'
file_dir = 'test'
file_name = "test_seq"
structure_file = open(f'output/{file_name}_structure.fasta', 'w')

with open(f'{file_dir}/{file_name}.fasta', 'r') as f:
    line = f.readline()
    line_count = 0
    while line and line_count < 4:
        if (line[0] == '>' or len(line) == 0):
            structure_file.write(line)
        else:
            structure_file.write(line)
            structure_file.write(secondary_structure(line[:-1])+'\n')
        line = f.readline()
        line_count += 1
structure_file.close()