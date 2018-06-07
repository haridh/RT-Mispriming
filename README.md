# RT-Mispriming
Pipeline to remove mispriming artifacts from RNA-seq datasets.

## Mispriming_finder_with_dups.py:

This is a script to identify sites of RT-mispriming from RNA-seq datasets. To use this script, you will need 4 files in the same folder as the python file.
1. Alignment file (bam) generated using a global aligner like BWA.
2. Genome fasta file from the same version used for alignment.
3. Sizes of the chromosomes in the genome as a bed file. For example if the length of chr1 is 249250617, the row corresponding to chr1 will be "chr1  0 249250617".
4. Bedfile of genomic features to filter out. This pipeline works best if reads from non-coding genome is filtered out.

There will be several intermediate files generated from this script. The file with suffix "Y_NAG" will contain the identified sites of mispriming. This file can be used as an input for Filter_misprimes.py.
```
usage: Mispriming_finder_with_dups.py [-h] fasta dnuc genome nprot
```
## Mispriming_finder_7mer_2mm.py:

This script works exactly like the above script except that it uses a more stringent cut-off. Instead of relying on 2 nucleotides ("dnuc" parameter), this requires you to provide 7 bases and the script will look for mispriming events flanking a sequence that matches the 7 bases with 2 mismatches allowed.

## Filter_misprimes.py:

This script is to extract all misprimed reads from mispriming sites identified using "Mispriming_finder_with_dups.py". You will need 2 input files:
1. Alignment (Bam) file converted to bed. This is one of the intermediate files generated from "Mispriming_finder_with_dups.py".
2. Bedfile of misprimed sites: The file with suffix "Y_NAG".
```
usage: Filter_misprimes.py [-h] Input_file Mispriming_events Output_file
```
