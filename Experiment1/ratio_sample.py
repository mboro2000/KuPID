import random
import sys

novel_input = sys.argv[1]
annot_input = sys.argv[2]
a = int(sys.argv[3])
output_folder = sys.argv[4]

f = open(novel_input, "r")
g = open(annot_input, "r")
o = open(output_folder + "/1:" + str(a) + ".css.fa", 'w')

n_lines = f.readlines()
a_lines = g.readlines()

for i in range(len(n_lines)):
    o.write(n_lines[i])
    if n_lines[i][0] == ">":
      num_novel += 1

to_select = a*num_novel  
for i in range(len(a_lines)):
    if a_lines[i][0] == ">":
        if random.randrange(num_annot) < to_select:
                o.write(a_lines[i])
                o.write(a_lines[i+1])
