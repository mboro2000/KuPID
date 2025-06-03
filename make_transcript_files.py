import pandas as pd
import numpy as np
import math
import sys
import fastaparser

depth_file = sys.argv[1]
fa = sys.argv[2]
output = sys.argv[3]

df = pd.read_csv(depth_file, sep="\t")

o_novel = open(output + "novel.transcript", "w")
o_annot = open(output + "annotated.transcript", "w")
o_novel_poly = open(output + "novel.polyA.transcript", "w")
o_annot_poly = open(output + "annotated.polyA.transcript", "w")

tail = ""
for i in range(30):
  tail += "A"

fa_dict = {}
novel = 0
annot = 0

with open(fa) as fasta_file:
  parser = fastaparser.Reader(fasta_file)
    for seq in parser:
      fa_dict[seq.id] = seq.sequence_as_string()

for row in df.itertuples(index=True, name='Row'):
  if "ENSG" in row[1]:
    o_novel.write(row[1] + "\t" + str(math.ceil(row[2])) + "\t" + str(0) + "\t" + fa_dict[row[1]] + "\n")
    o_novel_poly.write(row[1] + "\t" + str(math.ceil(row[2])) + "\t" + str(0) + "\t" + fa_dict[row[1]] + tail + "\n")
    novel += 1
  else:
    o_annot.write(row[1] + "\t" + str(math.ceil(row[2])) + "\t" + str(0) + "\t" + fa_dict[row[1]] + "\n")
    o_annot_poly.write(row[1] + "\t" + str(math.ceil(row[2])) + "\t" + str(0) + "\t" + fa_dict[row[1]] + tail + "\n")
    annot += 1

