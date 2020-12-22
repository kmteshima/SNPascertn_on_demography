import subprocess
import multiprocessing
import time
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import random
import re
import glob

msdatafile = "ms.txt"
eigenfile = "Is3_reseq.txt"

def parse_ms_data(filename):

    # regular expressions
    ms_match = re.compile('ms\s(\d+)\s(\d+)')
    seed_match = re.compile('seed|^\d+\s\d+\s\d+\s*$')
    blank_match = re.compile('^\s*$')
    graph_match = re.compile('^\(.+\);$')
    newdata_match = re.compile('//')
    segsites_match = re.compile('segsites:\s(\d+)')
    pos_match = re.compile('positions')

    #
    # init with fake values
    #
    samplesize = -1
    sn = 0

    pos = []
    graph = []

    # open data file
    with open(filename, 'r') as f:
        for line in f:

            # remove newline from the tail
            line = line.rstrip()

            # record sample_size and num_replication
            m = ms_match.search(line)
            if m:
                samplesize = int(m.group(1))
                # nrep = int(m.group(2))
                continue

            # skip random number seed
            if seed_match.search(line):
                continue

            # skip blank line
            if blank_match.search(line):
                continue

            # tree info
            if graph_match.search(line):
                graph.append(line)
                continue

            # pos info
            if pos_match.search(line):
                pos = line.split()[1:]
                continue

            # new data start
            if newdata_match.search(line):
                # initialize variables
                seq = []
                continue

            # record num of segsites
            m = segsites_match.search(line)
            if m:
                segsites = int(m.group(1))

                # when there is no variation,
                # returun array of '0'
                if segsites == 0:
                    # seq = ['0' for i in range(samplesize)]
                    seq = ['0'] * samplesize
                    sn = 0
                    # yield seq
                    yield {'seq': seq, 'pos': [], 'graph': None}

                continue

            # read and append a seq
            seq.append(line)
            sn += 1

            # when all samples are read
            if sn == samplesize:
                sn = 0
                # yield seq
                yield {'seq': seq, 'pos': pos, 'graph': graph}


def make_discovery_panel(msdata,
                         nsam_discovery,
                         nsam_type):
    # initialize
    discovery_panel = list()
    start = 0
    end = 0

    # extract seq of discovery panel
    for d, t in zip(nsam_discovery, nsam_type):

        end += d

        if d > 0:
            discovery_panel.extend(msdata[start: end])
        start += d + t
        end += t

        # break

    # return discovery_panel
    return discovery_panel

def candidate_snp_markers(discovery, MAF):

    nsam = len(discovery)

    # ms data -> list of list of int
    # calc sum of derived alleles at each site
    # snp_counts is a list of numbers of derived SNPs
    snp_counts = np.sum([[int(i) for i in list(m)]
                         for m in discovery],
                        axis=0)

    # return boolean list
    # True if 'polymorphic' AND '>MAF'
    return np.all([snp_counts > 0, snp_counts < nsam, snp_counts >= MAF*nsam, snp_counts<=(nsam-MAF*nsam)], axis=0)



def random_snp_markers(candidate_marker_bool_list, marker_size):

    # bool -> list of indexes where SNPs exist
    snp_index_list = np.where(candidate_marker_bool_list == True)[0]

    # list of indexes of randomly sampled SNP markers
    random_snp_indexes = np.sort(np.random.choice(snp_index_list,
                                                  marker_size,
                                                  replace=False))

    return random_snp_indexes


def typing(msdata, nsam_discovery, nsam_type, snp_marker):

    # initialize
    typed = list()
    start = 0
    end = 0

    # type each population
    for d, t in zip(nsam_discovery, nsam_type):

        # start and end indexes of sequences
        start += d
        end += d + t

        # type if typing sample size >0
        if t > 0:
            # read sample seq
            # take SNPs at marker site
            # char -> int
            typed.append([[int(m[i:i + 1])
                           for i in snp_marker]
                          for m in msdata[start:end]])
        start += +t

    # return list of list of int
    return typed


msdir = "../data/ms"
eigendir = "../data/eigen_input"

nsam_discovery = [100, 100, 100]
nsam_type = [100, 100, 100]
# maf_list = [0,0.01,0.05,0.1,0.2]
# marker_size = 50

for msfile in glob.glob("ms.txt"):
    # print(msfile)
    eigenfile = "Is3_reseq.txt"
    with open("{}".format(eigenfile), "w") as f:
        # print(eigenfile)
        print("ms 300 1 -t 20 -r 20 20000 -I 3 100 100 100 0.3", file=f)
        print("40972 1142 16277", file=f)

        for m in parse_ms_data(msfile):
            position = m["pos"]
            msdata = m["seq"]

            print("", file=f)
            print("//", file=f)
            print("segsites: {}".format(len(position)), file=f)
            position = " ".join(position)
            print("positions: ", position, file=f)
            for l in msdata[100:200]:
                print(l, file=f)
            for l in msdata[300:400]:
                print(l, file=f)
            for l in msdata[500:600]:
                print(l, file=f)

            # print(msdata[300:400], file=f)
            # print(msdata[500:600], file=f)

            # break
        # break

