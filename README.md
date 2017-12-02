# RT-Mispriming
Pipeline to remove mispriming artifacts from RNA-seq datasets.

This is a script to identify sites of RT-mispriming from RNA-seq datasets. To use this script, you will need 4 files in the same folder as the python file.
1. Alignment file (bam) generated using a global aligner like BWA.
2. Genome fasta file from the same version used for alignment.
3. Sizes of the chromosomes in the genome as a bed file. For example if the length of chr1 is 249250617, the row corresponding to chr1 will be "chr1  0 249250617".
4. Bedfile of genomic features to filter out. This pipeline works best if reads from non-coding genome is filtered out.

There will be several intermediate files generated from this script. The file with suffix "Y_NAG" will contain the identified sites of mispriming. This file can be used as an input for Filter_misprimes.py script to extract all misprimed reads.  
