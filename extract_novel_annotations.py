import sys

ase = sys.argv[1]
orig = sys.argv[2]
output = sys.argv[3]

gtf_as = open(ase, "r")
gtf_orig = open(orig, "r")

novel_only = open(output + "/gencode.v48.novel_only.gtf", "w")

gtf_orig_set = set()
gtf_orig_tids = set()
gtf_annot_tids = set()
for line in gtf_orig.readlines():
#       gtf_orig_set.add(line)
        if "transcript_id" in line:
                info = line.split("\t")[8]
                tid = info.split("transcript_id")[1].split('"')[1]
                gtf_orig_tids.add(tid)

novel_only_tids = set()
as_tids = set()
for line in gtf_as.readlines():
        if "transcript_id" in line:
                info = line.split("\t")[8]
          tid = info.split("transcript_id")[1].split('"')[1]
                as_tids.add(tid)
                if tid not in gtf_orig_tids:
                        novel_only_tids.add(tid)
                        print(tid)
                        novel_only.write(line)
                else:
                        gtf_annot_tids.add(tid)
