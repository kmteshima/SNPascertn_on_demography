# import module
import re
import numpy as np
import pandas as pd
import random

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


def make_typing_panel(msdata, n1, n2, lsam):

    # large sampling
    pop1_lsam = msdata[n1 - lsam:n1]
    pop2_lsam = msdata[n1 + n2 - lsam:n1 + n2]

    return pop1_lsam, pop2_lsam


def calc_within_pi(pop_typing):
    """calculation within pi

    pop1_typing : popA typing set
    pop2_typing : popB typing set

    -----------------------------------------------------------------------

    within pi was calculation from typing data
    """
    pop_sam = len(pop_typing[0])
    pop_combi = pop_sam * (pop_sam - 1) / 2

    return sum([(pop_sam - sum(pop_typing[i])) * sum(pop_typing[i])
                / pop_combi for i in range(len(pop_typing))])


def calc_between_pi(pop1_typing, pop2_typing):

    pop1_sam = len(pop1_typing[0])
    pop2_sam = len(pop2_typing[0])

    # pop1
    derived_pop1 = [sum(i) for i in pop1_typing]
    ans_pop1 = pop1_sam * np.ones(len(pop1_typing)) - derived_pop1

    # pop2
    derived_pop2 = [sum(i) for i in pop2_typing]
    ans_pop2 = pop2_sam * np.ones(len(pop2_typing)) - derived_pop2

    return sum(ans_pop1 * derived_pop2 + ans_pop2 * derived_pop1) / (pop1_sam * pop2_sam)


def main(lsam, n1, n2, msfile, save_dir, file_name):

    pop1_sfs = []
    pop2_sfs = []
    stats_data = []
    for msdata in [m["seq"] for m in parse_ms_data(msfile)]:
        # Sampling
        pop1_lsam, pop2_lsam = make_typing_panel(msdata, n1, n2, lsam)

        #make SNP list
        pop1_snp_list = [list(map(int, s))
                         for s in zip(*[list(m) for m in pop1_lsam])]
        pop2_snp_list = [list(map(int, s))
                         for s in zip(*[list(m) for m in pop2_lsam])]

        # calculation site frequency spectrum
        pop1_sfs.append([round(sum(i), 0) for i in pop1_snp_list])
        pop2_sfs.append([round(sum(i), 0) for i in pop2_snp_list])


        # calculation within pi
        wpi_1 = calc_within_pi(pop1_snp_list)
        wpi_2 = calc_within_pi(pop2_snp_list)

        # calculation between pi
        bpi_12 = calc_between_pi(pop1_snp_list, pop2_snp_list)

        stats_data.append([wpi_1, wpi_2, bpi_12])

        # continue

    pop1_df = pd.DataFrame(pop1_sfs)
    pop2_df = pd.DataFrame(pop2_sfs)

    df = pd.DataFrame(stats_data,
                      columns=["wpi1", "wpi2", "bpi12"])

    pop1_df.to_csv("{}/noasc_pop1_sfs_{}.csv".format(save_dir, file_name))
    pop2_df.to_csv("{}/noasc_pop2_sfs_{}.csv".format(save_dir, file_name))

    df.to_csv("{}/noasc_stats_{}.csv".format(save_dir, file_name))

if __name__ == "__main__":
    ### modify

    lsam = 100
    n1 = 100
    n2 = 100

    mrate = 0.1
    MAF = 0.1

    # input data
    msfile = "Is2_reseq_mrate01.txt"

    # save directory
    save_dir = ""
    file_name = "Is2_reseq_MAF{}".format(int(MAF * 100))

    main(lsam, n1, n2, msfile, save_dir, file_name)