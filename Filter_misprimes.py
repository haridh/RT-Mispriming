import argparse
parser=argparse.ArgumentParser(
    description='''Python script to identify all RT-misprimed reads''')
parser.add_argument('Input_file', help='Alignment bed file', nargs = 1)
parser.add_argument('Mispriming_events', help='Mispriming events file (File with suffix Y_NAG)', nargs = 1)
parser.add_argument('Output_file', help='Output file name without file type extension', nargs = 1)
args=parser.parse_args()
import pandas as pd
import pybedtools
import sys
orig = sys.argv[1]
inp1 = "Original bed file:" + " " + orig
print (inp1)
mis = sys.argv[2]
inp2 = "Mispriming events bed file:" + " " + mis
print (inp2)
out = sys.argv[3] + ".bed"
inp3 = "File name for misprimed reads to filter out:" + " " + out 
print (inp3)
df1 = pd.read_csv(orig, sep = "\t", header = None)
df1.columns = ['chr', 'start', 'stop', 'name', 'score', 'strand']
min1 = df1.loc[(df1.strand == '-')].copy()
min1.loc[:,('name')] = min1.loc[:,('chr')].astype(str) + ':' + min1.loc[:,('start')].astype(str)
plus1 = df1.loc[(df1.strand == '+')].copy()
plus1.loc[:,('name')] = plus1.loc[:,('chr')].astype(str) + ':' + plus1.loc[:,('stop')].astype(str)
df2 = min1.append(plus1)
df3 = pd.read_csv(mis, sep = "\t", header = None)
df3.columns = ['chr', 'start', 'stop', 'name', 'score', 'strand']
min2 = df3.loc[(df3.strand == '-')].copy()
min2.loc[:,'start'] = min2.loc[:,'start'] + 10
min2.loc[:,('name')] = min2.loc[:,('chr')].astype(str) + ':' + min2.loc[:,('start')].astype(str)
plus2 = df3.loc[(df3.strand == '+')].copy()
plus2.loc[:,'stop'] = plus2.loc[:,'stop'] - 10
plus2.loc[:,('name')] = plus2.loc[:,('chr')].astype(str) + ':' + plus2.loc[:,('stop')].astype(str)
df4 = min2.append(plus2)
li = df4['name'].tolist()
df5 = df2.ix[df2.name.isin(li)]
df5.to_csv(out, sep = "\t", header = False, quoting = False, doublequote = False, index = False)
