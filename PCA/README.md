# README: PCA

Text and scripts in this folder explains how to perform *ms*, apply typing process and conduct PCA with  EIGENSTRAT.


## Text and scripts in the folder

1. pca_ms_commands.txt
   - ms.txt
2. pca_make_input_single.py
3. pca_make_input_merged.py
4. pca_make_input_independent.py
5. pca_make_input_reseq.py
6. pca_single.py
7. pca_merged.py
8. pca_independent.py
9. pca_reseq.py
10. run_smartpca.perl (called from scripts 6-9)


## How to used

1. Perform *ms* as described in pca_ms_commands.txt. When *ms* is performed, you will obtain an output file (ms.txt).
2. Apply ascertainment process with scripts 2-5 depending on the ascertainment scheme.
   1. Set filenames in the head of scripts
      1. msdatafile: name of *ms* output data
      2. eigenfile: name of the output file of this script. This file becomes the input of the next script.
   2. execute the script
3. Conduct PCA by running EIGENSTRAT with scripts 6-9 depending on the ascertainment scheme.
    1. Set the input filename in the head of the script.
       1. eigenfile: input file name. This name should be the same as 'eigenfile' in the previous script.
       2. indfile: individual file name. This name should be the same as 'eigenfile' in the previous script.
    2. execute the script.
       - Results of PCA will be saved as pdf file.
       - Note: Scripts 6-9 use EIGENSOFT. py scripts call run_smarpca.perl to perform PCA. run_smarpca.perl uses smartpca.perl script distributed in EIGENSOFT package. Users should appropriately set the path to 'smartpca.perl' in run_smarpca.perl.
