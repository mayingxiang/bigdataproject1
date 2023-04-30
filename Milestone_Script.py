#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 29 19:36:01 2023

@author: petertian
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import time
import multiprocessing as mp

pd.set_option('display.max_columns', None)
# Cut Reads and Reference Bins into K-mers

k = 15
read_kmers = []  
reads = []  

with open('reads.fq', 'r') as f:
    for i, line in enumerate(f):
        if i % 4 == 1:  
            read = line.strip()  
            reads.append(read)  
            for j in range(len(read) - k + 1):
                kmer = read[j:j+k]  
                read_kmers.append(kmer)  

#print(len(read_kmers))

bin_kmers = []  
refs = []  

with open('reference_chr21_20000000_20050000.fa', 'r') as f:
    for line in f.readlines()[1:]:
        sequence = line.strip().split(',')[2]
        refs.append(sequence)  
        for j in range(len(sequence) - k + 1):
            kmer = sequence[j:j+k]  
            bin_kmers.append(kmer)  


#print(len(bin_kmers))


#Build a Distinct K-mer Set

all_kmers = read_kmers + bin_kmers
kmer_set = set(all_kmers)
kmers_distinct = len(kmer_set)

#print(len(all_kmers))
#print("Number of distinct kmers: ", kmers_distinct)


#Encode Reads and Reference Bins using One-Hot Encoding
# Reading Encoding
read_result = []
for kmer in kmer_set:
    current_result = []
    current_result.append(kmer)
    for read in reads:
        if kmer in read:
            current_result.append(1)
        else:
            current_result.append(0)
    read_result.append(current_result)   

columnnames = []
for i in range(0,2000):
    columnnames.append('Read'+str(i))
columnnames = ["kmer"] + columnnames
read_df = pd.DataFrame(read_result, columns=columnnames)

#print(read_df.head(`0))

# Reference Encoding
bin_result = []
for kmer in kmer_set:
    current_result = []
    current_result.append(kmer)
    for read in refs:
        if kmer in read:
            current_result.append(1)
        else:
            current_result.append(0)
    bin_result.append(current_result)  


columnnames = []
for i in range(0,500):
    columnnames.append('bin'+str(i))
columnnames = ["kmer"] + columnnames
bin_df = pd.DataFrame(bin_result, columns=columnnames)
#print(bin_df.head(10))


#MinHash
# For Sequential Purpose
def minhash(encode_result):
    encode_minhash = []
    for index in range(0,1000):
        np.random.seed(index)
        permuation = np.random.permutation(72530)
        curpermutation = []
        for bin in range(1,len(encode_result[0])):
            exist = 0
            for curindex in permuation:
                if encode_result[curindex][bin] == 1:
                    exist +=1
                    curpermutation.append(exist)
                    break
                else:
                    exist +=1
        encode_minhash.append(curpermutation)
    return encode_minhash

# Read sequential minhash
start_time = time.time()
readminhash_single = minhash(read_result)
end_time = time.time()
print("Sequential method time:", end_time - start_time, "seconds")

# Bin sequential minhash
start_time = time.time()
binminhash_single = minhash(bin_result)
end_time = time.time()
print("Sequential method time:", end_time - start_time, "seconds")


#For Parallel Purpose
def minhash_single(encode_result, loopnumber):
    encode_minhash = []
    for index in range(0,loopnumber):
     
        permuation = np.random.permutation(72530)
        curpermutation = []
        for bin in range(1,len(encode_result[0])):
            exist = 0
            for curindex in permuation:
                if encode_result[curindex][bin] == 1:
                    exist +=1
                    curpermutation.append(exist)
                    break
                else:
                    exist +=1
        encode_minhash.append(curpermutation)
    return encode_minhash

def minhash_parallel(encode_result):
    def collect_result(result):
        global encode_minhash
        encode_minhash.append(result)
    with mp.Pool(20) as pool:
        for i in range(0,20):
            pool.apply_async(encode_minhash, args=(encode_result, 50), callback=collect_result)
        pool.close()
        pool.join()
    return encode_minhash

# Read parallel minhash
encode_minhash = []
start_time = time.time()
readminhash_parallel = minhash_parallel(read_result)
end_time = time.time()
print("Parallel method time:", end_time - start_time, "seconds")
result = []
for i in readminhash_parallel:
  result = result + i
readminhash_parallel = result

# Bin Parallel minhash
encode_minhash = []
start_time = time.time()
binminhash_parallel = minhash_parallel(bin_result)
end_time = time.time()
print("Parallel method time:", end_time - start_time, "seconds")
result = []
for i in binminhash_parallel:
  result = result + i
binminhash_parallel = result

#Calculate Jaccard Similarity Based on Parallel
jaccard = []
for readnum in range(0,2000):
    maxbinsim = 0
    for binnum in range(0,500):
        num = 0
        for rownum in range(0,1000):
            if readminhash_parallel[rownum][readnum] == binminhash_parallel[rownum][binnum]:
                num += 1
        cursim = num / 1000
        if cursim > maxbinsim:
            maxbinsim = cursim
            maxbinname = 'bin' + str(20000000 + 100 * binnum) + "_" + str(20000100 + 100 * binnum)
    jaccard.append(["read" + str(readnum), maxbinname, maxbinsim])
#print(jaccard)


#Results and Evaluation Based on Parallel
df = pd.read_csv('read_position_benchmark.csv')
reference_end = df['reference_end'].tolist()
numbers_parallel = [int(x[1].split('_')[1]) for x in jaccard]

# Calculate the Pearson correlation between the two lists Based on Parallel
correlation_parallel, p_value_parallel = pearsonr(reference_end, numbers_parallel)
print("Pearson correlation: ", correlation_parallel)

# Calculate the mean square error Based on Parallel
mse_parallel = np.square(np.subtract(reference_end, numbers_parallel)).mean()
print("Mean square error: ", mse_parallel)


# Plot the scatter plot
#plt.scatter(reference_end, numbers_parallel,alpha=0.5)
#plt.xlabel("Benchmark Positions")
#plt.ylabel("Predicted Positions")
#plt.title(f"Scatter plot of positions (Pearson correlation: {correlation_parallel:.2f})")
#plt.show()



