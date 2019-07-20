import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def convert_to_Ts(seq):
    new_seq = ""
    for e in seq:
        if e == "U":
            new_seq += "T"
        else:
            new_seq += e
    return new_seq

def convert_seq_to_int_array(seq):
    ints = []
    for e in seq:
        if   e == 'A':
            ints.append(0)
        elif e == 'C':
            ints.append(1)
        elif e == 'G':
            ints.append(2)
        elif e == 'T':
            ints.append(3)
        else:
            ints.append(4)
    return ints

def get_best_alignment(seq_ints_1, seq_ints_2, threshold=20):
    best_score = 10000
    best_i = 0
    for i in range(len(seq_ints_1)):
        score = 0
        if len(seq_ints_2) >= len(seq_ints_1) - i:
            break
        for j in range(len(seq_ints_2)):
            if seq_ints_1[i+j] != seq_ints_2[j]:
                score += 1
        if score < best_score:
            best_i = i
            best_score = score

    return [best_i, best_score]

def get_mutation_profile(org_seq, org_seq_ints, reads):
    dist = [0 for x in range(len(org_seq_ints))]
    for r in reads:
        r_ints = convert_seq_to_int_array(r)
        a = get_best_alignment(r_ints, org_seq_ints)
        offset = a[0]
        for i in range(len(org_seq_ints)):
            if r_ints[i + offset] != org_seq_ints[i]:
                dist[i] += 1
    return dist


def get_included_seqs(filename):
    f = open(filename)
    seqs = [convert_to_Ts(l.rstrip()) for l in f.readlines()]
    f.close()
    return seqs

def get_seq_reads(filename, incuded_seqs):
    included_seq_ints = []
    for seq in included_seqs:
        included_seq_ints.append(convert_seq_to_int_array(seq))

    f = open(filename)
    lines = f.readlines()
    f.close()

    reads = []
    for i in range(len(lines) - 1):
        if lines[i][0] == "@" and lines[i + 1][0] != "@":
            reads.append(lines[i + 1].rstrip())
    #reads = reads[500:100000]

    i = 1
    scores = [[] for x in range(len(included_seq_ints))]
    seq_reads = [[] for x in range(len(included_seq_ints))]
    used_reads = 0
    j = -1
    for r in reads:
        j += 1
        best_score = 1000
        best_i = 0
        r_ints = convert_seq_to_int_array(r)
        i = 0
        for seq_int in included_seq_ints:
            scores[i] = get_best_alignment(r_ints, seq_int)
            if scores[i][1] < best_score:
                best_score = scores[i][1]
                best_i = i
            i += 1

        if best_score > 3:
            continue
        #used_reads += 1
        seq_reads[best_i].append(r)

    return seq_reads

included_seqs = get_included_seqs("included_seqs.csv")
seq_reads_dms = get_seq_reads("data/TAGTACTGCCAG/Sample1_S1_L001_R2_001_test.fastq", included_seqs)
seq_reads_nomod = get_seq_reads("data/AGTCGTGATGTT/Sample1_S1_L001_R2_001_test.fastq", included_seqs)
included_seq_ints = []
for seq in included_seqs:
    included_seq_ints.append(convert_seq_to_int_array(seq))


for i, seq in enumerate(included_seqs):
    print seq, len(seq_reads_dms[i]), len(seq_reads_nomod[i])
    #print seq, len(seq_reads_dms[i])
    exit()
    if i != 0:
        continue
    continue
    dist_dms = get_mutation_profile(seq, included_seq_ints[i], seq_reads_dms[i])
    dist_nomod = get_mutation_profile(seq, included_seq_ints[i], seq_reads_nomod[i])
    dist_substract = list(dist_dms)
    for i in range(len(dist_dms)):
        dist_substract[i] = dist_dms[i] - dist_nomod[i]
        if dist_substract[i] < 0:
            dist_substract[i] = 0

    fig, axs = plt.subplots(3, 1)

    sns.barplot(range(0, len(seq)), dist_dms, ax=axs[0])
    axs[0].set_xticks(np.arange(len(seq)))
    axs[0].set_xticklabels(list(seq))

    sns.barplot(range(0, len(seq)), dist_nomod, ax=axs[1])
    axs[1].set_xticks(np.arange(len(seq)))
    axs[1].set_xticklabels(list(seq))

    sns.barplot(range(0, len(seq)), dist_substract, ax=axs[2])
    axs[2].set_xticks(np.arange(len(seq)))
    axs[2].set_xticklabels(list(seq))
    plt.show()


    exit()


#print used_reads, len(reads)

