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
10. run_smartpca.perl


## How to used

1. Perform *ms* as described in pca_ms_commands.txt. When *ms* is performed, you will obtain an output file (ms.txt).
2. Apply ascertainment process with scripts 2-5.
   1. specify xxx in scripts
   2. execute the script
3. Conduct PCA by running EIGENSTRAT with scripts 6-9.
    1. specify xxx in scripts
    2. execute the script


## depricated

(9) is the EIGENSTRAT command used this research, which is called in (5)-(8) to execute the PCA.


(2)-(9) do not require input variables, but require a specific file in the directory where the script exists. (2)-(5) require ms.txt, the execution result of the ms command described in (1). Each of these will output one output file, which will be the data required for (6)-(9). (6)-(9) output the PCA results in pdf format based on the data in (2)-(5).

(11)は(1)の実行例であり、(2)-(5)に必要なファイルの例です。EIGENSTRATがコンパイルされた環境下で(2)-(5)を実行したのち、(6)-(9)を実行するとPCA結果であるPDFが出力されます。
