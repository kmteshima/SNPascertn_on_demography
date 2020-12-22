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

    # return discovery_panel
    return discovery_panel


def candidate_snp_markers(discovery_panel, MAF):

    # ms data -> list of list of int
    # calc sum of derived alleles at each site
    # snp_counts is a list of numbers of derived SNPs
    nsam = len(discovery_panel[0:100])

    snp_counts1 = np.sum([[int(i) for i in list(m)] for m in discovery_panel[0:100]], axis=0)
    snp_counts2 = np.sum([[int(i) for i in list(m)] for m in discovery_panel[100:200]], axis=0)
    snp_counts3 = np.sum([[int(i) for i in list(m)] for m in discovery_panel[200:300]], axis=0)

    snp1 = np.all([snp_counts1 > 0, snp_counts1 < nsam, snp_counts1 >= MAF * nsam, snp_counts1 <= (nsam - MAF * nsam)],
                  axis=0)
    snp2 = np.all([snp_counts2 > 0, snp_counts2 < nsam, snp_counts2 >= MAF * nsam, snp_counts2 <= (nsam - MAF * nsam)],
                  axis=0)
    snp3 = np.all([snp_counts3 > 0, snp_counts3 < nsam, snp_counts3 >= MAF * nsam, snp_counts3 <= (nsam - MAF * nsam)],
                  axis=0)

    candidate = np.any([snp1 == True, snp2 == True, snp3 == True], axis=0)

    # return boolean list
    # True if 'polymorphic' AND '>MAF'
    return candidate


def random_snp_markers(candidate_marker_bool_list, marker_size):

    # bool -> list of indexes where SNPs exist
    snp_index_list = np.where(candidate_marker_bool_list == True)[0]

    # list of indexes of randomly sampled SNP markers
    random_snp_indexes = np.sort(np.random.choice(snp_index_list,
                                                  marker_size,
                                                  replace=False))

    return random_snp_indexes

    # deplicated #
    #
    # initialize boolean list with all False
    # snp_bool_list = np.zeros(nsnp, dtype=bool)

    # positions of SNP markers are set True
    # snp_bool_list[random_snp_indexes] = True

    # return snp_bool_list


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


def calc_SFS(typing_data):

    return [round(sum(i), 0) for i in np.array(typing_data).T]


def calc_within_pi(typing_data):

    nsam = len(typing_data)
    comb = nsam * (nsam - 1) / 2.
    lseq = len(typing_data[0])

    derived = np.sum(typing_data, 0)
    ancestral = nsam * np.ones(lseq) - derived

    return np.sum(derived * ancestral) / comb


def calc_between_pi(pop1_typing, pop2_typing):

    pop1_sam = len(pop1_typing)
    pop2_sam = len(pop2_typing)

    # pop1
    derived_pop1 = np.sum(pop1_typing, axis=0)
    ans_pop1 = pop1_sam * np.ones(len(pop1_typing[0])) - derived_pop1

    # pop2
    derived_pop2 = np.sum(pop2_typing, axis=0)
    ans_pop2 = pop2_sam * np.ones(len(pop2_typing[0])) - derived_pop2

    return sum(ans_pop1 * derived_pop2 + ans_pop2 * derived_pop1) / (pop1_sam * pop2_sam)


def main(nsam_discovery, nsam_type, MAF, marker_size, msfile, save_dir, file_name):

    # make lists
    pop1_sfs = []
    pop2_sfs = []
    pop3_sfs = []
    stats_data = []

    # for msdata in [m["seq"] for m in parse_ms_data(msfile)]:
    for m in parse_ms_data(msfile):

        msdata = m["seq"]

        # sample discovery panel to generate SNP markers
        discovery_panel = make_discovery_panel(msdata,
                                               nsam_discovery,
                                               nsam_type)
        # make a list of candidate SNPs
        # return a boolean list (True if candidate SNPs, False otherwise)
        marker_cand = candidate_snp_markers(discovery_panel, MAF)

        # random choice SNP marker
        if sum(marker_cand) >= marker_size:
            # randomly sample SNPs
            # return list of indexes of markers
            snp_marker = random_snp_markers(marker_cand, marker_size)

            # genotype typing panels
            # return list of list of int
            pop1_typing, pop2_typing, pop3_typing = typing(msdata,nsam_discovery,nsam_type, snp_marker)

            # calculation site frequency spectrum
            pop1_sfs.append(calc_SFS(pop1_typing))
            pop2_sfs.append(calc_SFS(pop2_typing))
            pop3_sfs.append(calc_SFS(pop3_typing))

            # calculate within-pi
            wpi_1 = calc_within_pi(pop1_typing)
            wpi_2 = calc_within_pi(pop2_typing)
            wpi_3 = calc_within_pi(pop3_typing)

            # calculate between-pi
            bpi_12 = calc_between_pi(pop1_typing, pop2_typing)
            bpi_13 = calc_between_pi(pop1_typing, pop3_typing)
            bpi_23 = calc_between_pi(pop2_typing, pop3_typing)

            stats_data.append([wpi_1, wpi_2, wpi_3, bpi_12, bpi_13, bpi_23])

            # continue

    pop1_df = pd.DataFrame(pop1_sfs)
    pop2_df = pd.DataFrame(pop2_sfs)
    pop3_df = pd.DataFrame(pop3_sfs)

    df = pd.DataFrame(stats_data,
                      columns=["wpi1", "wpi2", "wpi3", "bpi12", "bpi13", "bpi23"])

    pop1_df.to_csv("{}/asc3_pop1_sfs_{}.csv".format(save_dir, file_name))
    pop2_df.to_csv("{}/asc3_pop2_sfs_{}.csv".format(save_dir, file_name))
    pop3_df.to_csv("{}/asc3_pop3_sfs{}.csv".format(save_dir, file_name))

    df.to_csv("{}/asc3_stats_{}.csv".format(save_dir, file_name))

if __name__ == "__main__":
    pass