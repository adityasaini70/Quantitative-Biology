from copy import deepcopy
import pandas as pd
import sys
import matplotlib.pyplot as plt

def create_2darray(l_1,l_2):
    array = [[0 for i in range(len(l_2))] for j in range(len(l_1))] 
    for i in range(len(l_1)):#row
        for j in range(len(l_2)):#col
            if(l_1[i]==l_2[j]):
                array[i][j] = 1
    return array

def find_maxval_col_row(array, col_idx, row_idx, flag):
    #flag = 0 : find max value in given row, by changing the column value
    #flag = 1 : find max value in given col, by changing the row value
    
    ans = -1 
    ans_1 = -1
    idx = []
    
    if flag == 0:
        for col in range(col_idx, len(array[0])):
            ans = max(ans, array[row_idx][col])
            if ans_1 != ans:
                idx = [row_idx, col]
            ans_1 = ans
    else:
        #flag = 1
        for row in range(row_idx, len(array)):
            ans = max(ans, array[row][col_idx])
            if ans_1 != ans:
                idx = [row, col_idx]
            ans_1 = ans
    return [ans, idx]

def nw_algo(array):
    ans = deepcopy(array)
    row_max = len(ans) - 1
    col_max = len(ans[0]) - 1
    for row_idx in range(row_max, -1, -1):
        for col_idx in range(col_max, -1, -1):
            max_val = 0
            if(col_idx + 1 < len(ans[0]) and row_idx + 1 < len(ans)):
                #Diagonal case
                max_val = max(max_val, ans[row_idx+1][col_idx+1])
            if(col_idx + 2 < len(ans[0]) and row_idx + 1 < len(ans)):
                #Down a row making column gap
                max_val = max(max_val, find_maxval_col_row(ans, col_idx + 2, row_idx + 1, 0)[0])
            if(col_idx + 1 < len(ans[0]) and row_idx + 2 < len(ans)):
                #Down a column making row gap
                max_val = max(max_val, find_maxval_col_row(ans, col_idx + 1, row_idx + 2, 1)[0])           
            ans[row_idx][col_idx] += max_val
    return ans

def align(row, col, array, ans):
    if (row > len(array) or col > len(array[0])):
        return
    
    max_val = max(find_maxval_col_row(array, col + 1, row, 0)[0], find_maxval_col_row(array, col, row + 1, 1)[0], array[row][col])
    
    if(max_val == array[row][col]):
        aligned_s.append(s[col])
        aligned_t.append(t[row])
        if(row + 1 < len(array) and col + 1 < len(array[0])):
            align(row+1, col+1, array, ans)
    
    elif(max_val == find_maxval_col_row(array, col + 1, row, 0)[0]):
        #Down a row, skipping one column at a time
        col_range = find_maxval_col_row(array, col + 1, row, 0)[1][1]
        # print(col_range)
        for i in range(col, col_range + 1):
            aligned_s.append(s[i])
            aligned_t.append('-')
            # print(aligned_s, aligned_t)
        aligned_t[-1] = t[row]
        if(row + 1 < len(array) and col_range + 1 < len(array[0])):
            align(row + 1, col_range + 1, array,ans)
        
        
    elif(max_val == find_maxval_col_row(array, col, row + 1, 1)[0]):
        row_range = find_maxval_col_row(array, col, row + 1, 1)[1][0]
        for i in range(row, row_range + 1):
            aligned_t.append(t[i])
            aligned_s.append('-')
        aligned_s[-1] = s[col]
        if(row_range + 1 < len(array) and col + 1 < len(array[0])):
            align(row_range + 1, col + 1, array , ans)
    




#Extracting data from the given FASTA file
input = sys.argv[2]
file = open(input,'r')
file_content = file.readlines()
s = list(file_content[1])
t = list(file_content[4])
# s =['A','B','C','N','Y','R','Q','C','L','C','R','P','M'] #For testing purposes
# t =['A','Y','C','Y','N','R','C','K','C','R','B','P'] #For testing purposes
s.pop()
t.pop()
# print(len(s), len(t))
# print(s)
# print(t)


#Creating the similarity/identity matrix (zero gap penalty)
sim_matrix = create_2darray(t,s)

#Applying the NW algo to the similarity matrix
nw_matrix = nw_algo(sim_matrix)

aligned_s = [] #Will contain the aligned version of s
aligned_t = [] #Will contain the aligned version of t
align(0, 0, nw_matrix, sim_matrix)
print('Input : ')
print('1st Protein sequence =>',''.join(s))
print('2nd Protein sequence =>',''.join(t))
out_s = ''.join(aligned_s)
out_t = ''.join(aligned_t)

print('Aligned output : ')

b = ''
for i in range(len(out_s)):
	if(out_s[i]==out_t[i]):
		b += '|'
	else :
		b += ' '

print(out_s)
print(b)
print(out_t)
# print(len(s), len(t))
# print(aligned_s)
# print(aligned_t)
# print(nw_matrix[0])

print('Sum matrix : ')
print(nw_matrix)

print('Dotplot')
print(sim_matrix)
