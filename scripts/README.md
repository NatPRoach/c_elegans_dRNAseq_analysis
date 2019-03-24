#Scripts directory README

Directions for replicating the analysis in "The full length transcriptome of C. elegans".
*At the moment the file paths in these scripts are explicit paths to files on my machine. This will change as I update the files for better & easier reproducibility. Until then you may have to change the file paths to reflect the directory structure of your machine.*

If you plan to replicate the poly(A) profiling estimation, start in the data directory and run download\_fast5.bsh to download the necessary data and set up the directory structure required. Then start in this directory by running the scripts in 00\_basecalling\_and\_adapter\_trimming/ and making your way through the directories in the numbered order they appear. 

If you plan to only replicate the clustering and UTR estimation, start in the data directory and run download\_fastq.bsh to download the necessary data and set up the directory structures required. Then start in this directory by running the scripts in 01\_alignment/ and making your way through the directories in the numbered order they appear. 

