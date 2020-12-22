# README: three-population model

## Text and scripts of the folder

0. 3pop_ms_command.txt
1. 3pop_single.py: single population scheme
2. 3pop_Independent.py: independent panel scheme
3. 3pop_merged.py: merged panel scheme
4. 3pop_reseq.py: re-sequence data


## How to use

1. Perform *ms* simulation as described in 3pop_ms_command.txt
1. Execute .py scripts
    - scripts 1-3
        ```
        main(nsam_discovery, nsam_type, MAF, marker_size, msfile, save_dir, file_name)
        ```
        - nsam_discovery: number of individuals in each population samples as a discovery panel. (eg. [100,0,0])
        - nsam_type: number of individuals in each population samples as a typing panel. (eg. [100,100,100])
        - MAF: threshold frequency of the marker selection
        - markersize: number of SNP markers to be selected
        - msfile: input ms data
        - save_dir: directory of output file
        - file_name: name of output file
    - script 4
        ```
        main(lsam, n1, n2, n3, msfile, save_dir, file_name)
        ```
        - lsam: the number of individuals sampled from each population as a typing sample in (1)-(3).
        - n1: the number of individuals sampled from population 1 in (1)-(3).
        - n2: the number of individuals sampled from population 2 in (1)-(3).
        - n3: the number of individuals sampled from population 3 in (1)-(3).
        - msfile: input ms data
        - save_dir: directory of output file
        - file_name: name of output file
1. Outputs
   - three csv files are created as output files
   - A file with "sfs" in the file name: site frequency spectrum of the population for each replicate
   - A file with "stats" in the file name: $\pi_w$ and $\pi_b$ for each replicate.
