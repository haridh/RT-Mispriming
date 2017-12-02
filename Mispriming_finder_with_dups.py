import argparse
parser=argparse.ArgumentParser(
    description='''Python script to identify RT-mispriming artifacts from RNA-seq datasets''')
parser.add_argument('fasta', help='Fasta file for the genome used to align the dataset', nargs = 1)
parser.add_argument('dnuc', help='First two nucleotides of the adapter', nargs = 1)
parser.add_argument('genome', help='Chromosome sizes as bed file', nargs = 1)
parser.add_argument('nprot', help='Bedfile of genomic regions to filter out', nargs = 1)
args=parser.parse_args()
import glob
import numpy
import pandas as pd
import pybedtools
import sys
fl = glob.glob("*.bam")
fasta = sys.argv[1]
inp1 = "Genome fasta file:" + " " + fasta
print (inp1)
dnuc = sys.argv[2]
inp2 = "Looking for stops at:" + " " + dnuc
print(inp2)
Nuc = int(len(dnuc))
inp3 = "Adapter match count:" + " " + str(Nuc)
print (inp3)
genome = sys.argv[3]
inp4 = "Edited genome sizes:" + " " + genome
print(inp4)
nprot = sys.argv[4]
inp5 = "Non-protein coding genes bed file:" + " " + nprot
print(inp5)
for f in fl:   
   fl2 = f.strip('.bam')
   fl2a = fl2 + '.bed'
   bam = pybedtools.BedTool(f)
   bed = bam.bam_to_bed(output=fl2a)
   nonprot = pybedtools.BedTool(nprot)
   fl2b = fl2 + '_protcod.bed'
   prot = bed.intersect(nonprot, v=True, s = True, output = fl2b)
   fl3 = fl2 + '_protcod_edit.bed'
   df = pd.read_csv(fl2b, sep = "\t", header = None)
   df.columns = ['chr', 'start', 'stop', 'name', 'score', 'strand']
   df.loc[:,['name']] = 1
   df.loc[:,['score']] = 2
   df2a = df
   cp1 = df2a.loc[(df2a.strand == '-') & (df2a.start >= 10)].copy()
   cp1.loc[:,'stop'] = cp1.loc[:,'start']
   cp1.loc[:,'start'] = cp1.loc[:,'stop'] - int(Nuc)
   cp2 = df2a.loc[(df2a.strand == '+') & (df2a.stop >= 10)].copy()
   cp2.loc[:,'start'] = cp2.loc[:,'stop']
   cp2.loc[:,'stop'] = cp2.loc[:,'start'] + int(Nuc)
   df2 = cp1.append(cp2)
   df2.loc[:,('name')] = df2.loc[:,('chr')].astype(str) + ':' + df2.loc[:,('start')].astype(str) + '-' + df2.loc[:,('stop')].astype(str)
   df2.to_csv(fl3, sep = "\t", header = False, quoting = False, doublequote = False, index = False)
   df2b = pybedtools.BedTool(fl3)
   df2c = pybedtools.BedTool(genome)
   k = df2b.intersect(df2c, u=True)   
   fl4 = fl2 + '_protcod_edit_fasta.bed'
   p = k.sequence(fi = fasta, s = True, tab = True, name = True)
   b = p.save_seqs(fl4)
   df4 = pd.read_csv(fl4, sep = "\t", header = None)
   df4.columns = ['name', 'dinuc']
   df4.loc[:,('dinuc')] = df4.loc[:,('dinuc')].str.upper()
   df4.to_csv(fl4, sep = "\t", header = False, quoting = False, doublequote = False, index = False)
   Match = "Y"
   df5 = df4.ix[df4.dinuc == dnuc, ['name', 'dinuc']]
   df6 = df5.groupby(['name', 'dinuc']).size().reset_index(name = "count")
   df7 = df6.ix[df6['count'] >= 10, ['name', 'dinuc']]
   li = df7['name'].tolist()
   df8 = df2.ix[df2.name.isin(li)]
   df9 = df8.drop_duplicates()
   fl5 = fl2 + '_protcod_mispriming_fasta_coor' + '_' + Match + dnuc + '.bed'
   cp3 = df9.loc[df9.strand == '-'].copy()
   cp3.loc[:,'stop'] = cp3.loc[:,'stop'] + 10
   cp3.loc[:,'start'] = cp3.loc[:,'start'] - 8
   cp4 = df9.loc[df9.strand == '+'].copy()
   cp4.loc[:,'start'] = cp4.loc[:,'start'] - 10
   cp4.loc[:,'stop'] = cp4.loc[:,'stop'] + 8
   cpb = cp3.append(cp4)
   cpb.to_csv(fl5, sep = "\t", header = False, quoting = False, doublequote = False, index = False)
   misprimeY = pybedtools.BedTool(fl5)
   Match = "N"
   df5 = df4.ix[df4.dinuc != dnuc, ['name', 'dinuc']]
   df6 = df5.groupby(['name', 'dinuc']).size().reset_index(name = "count")
   df7 = df6.ix[df6['count'] >= 10, ['name', 'dinuc']]
   li = df7['name'].tolist()
   df8 = df2.ix[df2.name.isin(li)]
   df9 = df8.drop_duplicates()
   fl5 = fl2 + '_protcod_mispriming_fasta_coor' + '_' + Match + dnuc + '.bed'
   cp3 = df9.loc[df9.strand == '-'].copy()
   cp3.loc[:,'stop'] = cp3.loc[:,'stop'] + 10
   cp3.loc[:,'start'] = cp3.loc[:,'start'] - 8
   cp4 = df9.loc[df9.strand == '+'].copy()
   cp4.loc[:,'start'] = cp4.loc[:,'start'] - 10
   cp4.loc[:,'stop'] = cp4.loc[:,'stop'] + 8
   cpb = cp3.append(cp4)
   cpb.to_csv(fl5, sep = "\t", header = False, quoting = False, doublequote = False, index = False)
   misprimeN = pybedtools.BedTool(fl5)
   fl6 = fl2 + '_protcod_mispriming_fasta_coor' + '_' + 'Y_N' + dnuc + '.bed'
   mis = misprimeY.intersect(misprimeN, v=True, output = fl6)
   misfasta = mis.sequence(fi = fasta, s = True, name = True)
   fl6 = fl2 + '_mispriming_'+ dnuc + 'fasta.txt'
   misfastasave = misfasta.save_seqs(fl6)

