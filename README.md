# c\_elegans\_dRNAseq\_analysis
Github repository for scripts related to the analysis of _C elegans_ dRNA sequencing data.

Requirements:

Python:
- python2
- matplotlib
- seaborn
- pysam
- pybedtools
- numpy
- biopython
- scikit-learn
- rpy2
- pygenometracks


R:
- R
- ggplot2
- cowplot
- here
- scales
- eulerr
- gdata

Other:
- bedtools
- samtools
- minimap2
- cpat

If you want to replicate basecalling and poly(A) tail calling you will also need:
- poreplex
- albacore
- nanopolish

To replicate analysis downstream of basecalling and poly(A) tail length calling, run master\_script.bsh
To replicate poly(A) calling with nanopolish, and basecalling with poreplex and albacore, you will need to modify master\_script.bsh by commenting / uncommenting certain sections that are labeled in the file. Be warned that poly(A) calling and basecalling take a long time to run and require a very large amount of data storage, as the requisite fast5 files are quite large.

To regenerate the metagene data used to used to generate figure 1B, you will need to uncomment a line in scripts/08\_make\_figures/figure1.bsh, which is labeled. By default, the script uses precomputed metagene data included in this GitHub.