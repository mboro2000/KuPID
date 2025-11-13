import pandas as pd
import numpy as np
import sys, getopt
from Bio import SeqIO

input = None
method = None
output = None
scale_path = None
l = None

argv = sys.argv[1:]
options = 'amosl'
long_options = ['abundances', 'method', 'output', 'scale', 'l']
try:
    opts, args = getopt.getopt(args, options, long_options)
except getopt.error as err:
    print(str(err))

for opt, arg in opts:
    if opt in ('-a', '--abundances'):
        input = arg
    else if opt in ('-m', '--method'):
        method = arg
    else if opt in ('-o', '--output'):
        output = arg
    else if opt in ('-s', '--scale'):
        scale = arg
    else if opt in ('l', '--l'):
        l = arg

if method == 'stringtie2':
    dict = {}
    rows = 0

    f = open(input, "r")
    for line in f.readlines():
        data = line.split("\t")
        if len(data) > 2:
        if data[2] == 'transcript':
            if len(data[8].split("reference_id ")) > 1:
                rows += 1
                tid = data[8].split("reference_id ")[1].split('"')[1]
                tpm = float(data[8].split("TPM ")[1].split('"')[1])
                dict['row_' + str(rows)] = [tid, tpm]
    df = pd.DataFrame.from_dict(dict, orient='index', columns=['Transcript', 'TPM_all'])
    total_r = df['TPM_all'].sum()
    df['TPM'] = 1000000 * df['TPM_all'] / total_r
    df = df[df['TPM'] > 0]

if method == 'flair':
    df = pd.read_csv(input, sep='\t')
    df['TPM_all'] = df['sample1_condition1_batch1']
    total_r = df[df['ids'].str.contains('ENST')]['TPM_all'].sum()
    df['TPM'] = 1000000 * df['TPM_all'] / total_r
    df = df[df['ids'].str.contains('ENST')]
    df['Transcript'] = df['ids'].str.split("_").apply(lambda x: x[0])
    df = df[df['TPM'] > 0]

if method == 'IsoQuant':
    df = pd.read_csv(input, sep='\t')
    df = df[df['#feature_id'].str.contains('ENST')]
    df['Transcript'] = df['#feature_id']
    df = df[df['TPM'] > 0]

scale = pd.read_csv(scale_path)
scale['Transcript'] = scale['Transcript'].str.split("|").apply(lambda x: x[0])
scaled = pd.merge(df, scale, on='Transcript', how='outer')
scaled.fillna(0, inplace=True)
li = scaled[(scaled['TPM'] > 0) & (scaled['Group Count'] >= l)]['TPM'].min()
scaled['Scaled Ai'] = np.where(scaled['Group Count'] >= l, (scaled['TPM'] - li) + li*scaled['Scale'], scaled['TPM'])
scaled['Scaled Ai'] = np.where(scaled['Scaled Ai'] < 0, 0, scaled['Scaled Ai'])
scaled['Scaled Ai'] = np.where(scaled['TPM'] == 0, 0, scaled['Scaled Ai'])
total_Ai = scaled['Scaled Ai'].sum()
scaled['Scaled TPM'] = 1000000 * scaled['Scaled Ai'] / total_Ai
scaled.to_csv(output)
