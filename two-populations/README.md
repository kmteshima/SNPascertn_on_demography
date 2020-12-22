# README: two-population model


## Text and scripts in this folder

0. 2pop_ms_command.txt
1. 2pop_single.py
2. 2pop_Independent.py
3. 2pop_merged.py
4. 2pop_reseq.py
5. S5_single.py
6. S5_merged.py
7. S5_Independent.py

## How to use

1. *ms* simulation
   - Perform *ms* as described in 2pop_ms_command.txt.
2. ascertainment process
   - scripts  
      - 2pop_single.py reproduces the single population scheme  
      - 2pop_Independent.py reproduces the independent panel scheme  
      - 2pop_merged.py reproduces the merged panel scheme  
      - 2pop_reseq.py reproduces the re-sequence
   - usage
      - scripts 1-3: please modify arguments of the 'main' function in the scripts, and execute.
        ```
        main(nsam_discovery, nsam_type, MAF, marker_size, msfile, save_dir, file_name)
        ```
        - nsam_discovery: number of individuals in each population samples as a discovery panel. (eg. [100,0])
        - nsam_type: number of individuals in each population samples as a typing panel. (eg. [100,100])
        - MAF: threshold frequency of the marker selection
        - markersize: number of SNP markers to be selected
        - msfile: input ms data
        - save_dir: directory of output file
        - file_name: name of output file
      - scripts 4: set the following arguments and execute.
        ```
        main(lsam, n1, n2, msfile, save_dir, file_name)
        ```
        - lsam: the number of individuals sampled from each population as a typing sample in (1)-(3).
        - n1: the number of individuals sampled from population 1 in (1)-(3).
        - n2: the number of individuals sampled from population 2 in (1)-(3).
        - msfile: input ms data
        - save_dir: directory of output file
        - file_name: name of output file
      - outputs
        - three csv files will be created as outputs.
        - A file with "sfs" in the file name: site frequency spectrum of the population for each replicate.
        - A file with "stats" in the file name: $\pi_w$ and $pi_b for each replicate.
3. typing process
   - scripts
      - S5_single.py  
      - S5_merged.py  
      - S5_Independent.py  
   - usage
        ```
        main(pop1_small, pop1_large, pop2_small, pop2_large, replicate, theta, fNr, mrate, rsite, MAF, marker_size, save_dir, AP_num, file_name)
        ```
        - pop1_small: number of individuals in population 1 samples as a discovery panel.
        - pop1_large: number of individuals in population 1 samples as a typing panel.
        - pop2_small: number of individuals in population 2 samples as a discovery panel.
        - pop2_large: number of individuals in population 2 samples as a typing panel.
        - replicate: ms parameters, in this case replicate
        - theta: mutation parameter used in ms, $\Theta=4Nu$
        - fNr: recombination parameter used in ms, $\rho = 4Nr$
        - mrate: migration parameter in ms, $4Nm$
        - rsite: number of recombination site
        - MAF: threshold frequency of the marker selection
        - marker_size: number of SNP markers to be selected
        - save_dir: directory of output file
        - AP_num: Indicators that identify the ascertainment scheme
        - file_name: name of output file
   - output
     - three csv files will be created. Formats are the same as ascertainment process (outputs of scripts 1-4).
