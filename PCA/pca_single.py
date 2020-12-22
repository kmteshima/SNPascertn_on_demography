import subprocess
import re
import numpy as np
import pprint


def run_command(cmd):

    #print(cmd)

    try:
        res = subprocess.check_output(cmd,
                                      stderr=subprocess.PIPE,
                                      shell=True)

        print(res.decode())

    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.cmd)
        print(e.output.decode())


def parse_ms_data(filename):
    """
    generator function

    input:
        filename of ms-output

    yield:
        dict of 0-1 string sequences, pos, graph
    """

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
    # nrep = -1
    # segsites = -1
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


nrep = 0
n1 = 50
n2 = 50
n3 = 50
nind = 150

for rep in parse_ms_data("Is3_single.txt"):

    #
    # .geno
    #

    # list of str -> np.array of list of int
    lseq = np.array([list(map(int, i))
                     for i in list(map(list, rep["seq"]))])

    # haploid data -> individual genotype data
    # 0: alt homo, 1: hetero, 2: ref homo
    genotype = np.array([2 - (lseq[i * 2] + lseq[i * 2 + 1])
                         for i in range(nind)])

    with open(f"rep{nrep:02}.geno", "w") as f:
        for i in genotype.T:
            f.write("".join(list(map(str, i))) + "\n")

    print()

    #
    # .snp
    #
    with open(f"rep{nrep:02}.snp", "w") as f:
        nsnp = 0
        for p in rep["pos"]:
            f.write(f"prs{nsnp:05}\t1\t{float(p) / 100:.6f}\t{int(float(p) * 1000000)}\n")
            nsnp += 1

    break

with open("mssamples.ind", "w") as f:
    for i in range(nind):

        if i < n1:
            f.write(f"sample{i:02}\tU\tPop1\n")
        elif i < n1 + n2:
            f.write(f"sample{i:02}\tU\tPop2\n")
        else:
            f.write(f"sample{i:02}\tU\tPop3\n")

command = "./run_smartpca.perl 00"
run_command(command)