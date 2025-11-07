from Bio import SeqIO

ref_fa = open("GRCh38.chr1-22.fa", "w")
g = open("gencode.v48.annotation.gtf", "r")
ref_gtf = open("gencode.v48.annotation.chr1-22.gtf", "w")

chr_ids = set()
for i in range(22):
  with open("Homo_sapiens.GRCh38.dna.chromosome." + str(i+1) + ".fa") as fasta_file:
    for record in SeqIO.parse(fasta_file, "fasta"):
      rid = record.id
      record.id = "chr" + str(i+1)
      record.description = "chr" + str(i+1)
      ref_fa.write(record.format("fasta"))
  
  chr_ids.add("chr" + str(i+1))

for line in g.readlines():
  chr = line.split("\t")[0]
  if chr in chr_ids:
    ref_gtf.write(line)
