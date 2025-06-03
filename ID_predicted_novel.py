import sys

ID_tracking = sys.argv[1]
ID_annotations = sys.argv[2]
output_folder = sys.argv[3]

f = open(ID_tracking, "r")
g = open(ID_annotations, "r")
q = open(output_folder, "w")

matches = set()

for line in f.readlines():
        info = line.split("\t")
        if info[3] == "=":
                tid = info[4].split("|")[1]
                matches.add(tid)

for line in g.readlines():
        feature_info = line.split("transcript_id ")
        if len(feature_info) > 1:
                tid = feature_info[1].split('"')[1]
                if tid not in matches:
                        q.write(line)
