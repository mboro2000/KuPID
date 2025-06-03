import random
import sys

novel_input = sys.argv[1]
s = int(sys.argv[2])
output_folder = sys.argv[3]

f = open(novel_input, "r")
o = open(output_folder + "/novel_" + str(s) + "%.css.fa", 'w')

lines = f.readlines()
for i in range(len(lines)):
        if lines[i][0] == ">":
                if random.randrange(100) < s:
                        o.write(lines[i])
                        o.write(lines[i+1])
