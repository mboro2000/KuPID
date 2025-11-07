import sys
input = sys.argv[1]
header = sys.argv[2]
f = open(input, "r")
lines = f.readlines()
f.close()
f = open(input, "w")
for line in lines:
        if line[0] == ">":
                f.write(">" + header + "\n")
        else:
                f.write(line)
