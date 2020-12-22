
# module
import subprocess
import re
import numpy as np
import pandas as pd
import random

def run_command(cmd):
    try:
        res = subprocess.check_output(cmd,stderr=subprocess.PIPE,shell=True)
        print(res.decode())
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.cmd)
        print(e.output.decode())


def parse_ms_data(filename):
    """ Parse ms-output file

    generator function

    Args:
        filename (str): path to the ms-output file

    yield:
        dict: 0-1 string sequences, pos, graph
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


def make_ms_command(pop1_small, pop1_large,
                    pop2_small, pop2_large,
                    replicate, theta, mrate, fNr, rsite,
                    save_dir, AP_num, file_name):
    # initial value
    n1 = pop1_small + pop1_large
    n2 = pop2_small + pop2_large

    # make command list
    commands = []

    # command
    cmd = "ms {} {} ".format(n1 + n2, replicate)
    cmd += "-t {} ".format(theta)
    cmd += "-r {} {} ".format(fNr, rsite)
    cmd += "-I 2 {} {} {} ".format(n1, n2, mrate)
    cmd += "> {}/AP{}BI_SFS_{}.txt".format(save_dir, AP_num, file_name)

    return cmd


def small_sampling(msdata, pop1_small, pop1_large, pop2_small, pop2_large):
    # small sampling
    discovery_panel = np.r_[msdata[pop1_small:pop1_small + pop1_large]]

    return discovery_panel


def large_sampling(msdata, pop1_small, pop1_large, pop2_small, pop2_large):
    # large sampling
    pop1_lsam = msdata[pop1_small:pop1_small + pop1_large]
    pop2_lsam = msdata[pop1_small + pop1_large + pop2_small:
                       pop1_small + pop1_large + pop2_small + pop2_large]

    return pop1_lsam, pop2_lsam


def make_snp_marker(discovery, MAF):
    # make SNP marker list at random
    site_allele_freq = [np.mean(list(map(int, s)))
                        for s in zip(*[list(m) for m in discovery])]
    site = np.array([s for s in range(len(site_allele_freq))])

    return [s for s in site if abs(0.5 - site_allele_freq[s]) <= 0.5 - MAF]


def choice_snp_marker(marker, marker_size):
    import random

    return random.sample(marker, marker_size)


def typing(typing_panel, marker):
    # typing
    return [[list(map(int, s)) for s in zip(*[list(m) for m in typing_panel])][int(i)]
            for i in marker]


def calc_within_pi(pop_typing):
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


def main(pop1_small, pop1_large, pop2_small, pop2_large,
         replicate, theta, fNr, mrate, rsite,
         MAF, marker_size,
         save_dir, AP_num, file_name):

    cmd = make_ms_command(pop1_small, pop1_large, pop2_small, pop2_large,
                          replicate, theta, mrate, fNr, rsite,
                          save_dir, AP_num, file_name)

    run_command(cmd)

    msfile = "{}/AP{}BI_SFS_{}.txt".format(save_dir, AP_num, file_name)

    # make lists
    pop1_sfs = []
    pop2_sfs = []
    stats_data = []

    for msdata in [m["seq"] for m in parse_ms_data(msfile)]:
        # Sampling
        discovery_panel = small_sampling(msdata,
                                         pop1_small, pop1_large, pop2_small, pop2_large)

        pop1_lsam, pop2_lsam = large_sampling(msdata,
                                              pop1_small, pop1_large, pop2_small, pop2_large)

        # Make SNP Marker
        marker_cand = make_snp_marker(discovery_panel, MAF)

        # random choice SNP marker
        if len(marker_cand) >= marker_size:
            snp_marker = choice_snp_marker(marker_cand, marker_size)

            # genotyping
            pop1_typing = typing(pop1_lsam, snp_marker)
            pop2_typing = typing(pop2_lsam, snp_marker)

            # change minor allele frequency
            pop1_sfs.append([round(sum(i),0) for i in pop1_typing])
            pop2_sfs.append([round(sum(i),0) for i in pop2_typing])

            # calculation within pi
            pop1_wpi = calc_within_pi(pop1_typing)
            pop2_wpi = calc_within_pi(pop2_typing)

            # calculation between pi
            bpi = calc_between_pi(pop1_typing, pop2_typing)

            # make statistic data
            stats_data.append([pop1_wpi, pop2_wpi, bpi])

    pop1_df = pd.DataFrame(pop1_sfs)
    pop2_df = pd.DataFrame(pop2_sfs)
    stats_df = pd.DataFrame(stats_data, columns=["pop1_wpi", "pop2_wpi", "bpi"])

    pop1_df.to_csv("{}/AP{}asc2_pop1_SFS_{}.csv".format(save_dir, AP_num, file_name))
    pop2_df.to_csv("{}/AP{}asc2_pop2_SFS_{}.csv".format(save_dir, AP_num, file_name))
    stats_df.to_csv("{}/AP{}asc2_stats_{}.csv".format(save_dir, AP_num, file_name))


if __name__ == "__main__":
    ### modify

    # initial value
    pop1_small = 0
    pop1_large = 100
    pop2_small = 0
    pop2_large = 100
    replicate = 10000
    theta = 20
    fNr = 20
    rsite = 200000
    marker_size = 50

    mrate = 0.1

    save_dir = ""

    MAF =0.1

    AP_num = "s"
    file_name = "S5_single_MAF{}".format(int(MAF * 100))

    main(pop1_small, pop1_large, pop2_small, pop2_large, replicate, theta, fNr, mrate, rsite, MAF, marker_size, save_dir, AP_num, file_name)


