import sys, os

in_file = sys.argv[1]
out_file = sys.argv[2]

file = open(in_file)
N = 0
v = {}
line = file.readline()
while line:
    l = line.rstrip().split('\t')
    N += int(l[3])
    if l[3] not in v:
        v[l[3]] = int(l[3])
    line = file.readline()
file.close()

v_n = {}
for i in v:
    v_n[i] = str(round(v[i]*1000000/N, 3))

file = open(in_file)
out = open(out_file, "w")
line = file.readline()
while line:
    l = line.rstrip().split('\t')
    l[3] = v_n[l[3]]
    out.write('\t'.join(l) + '\n')
    line = file.readline()
file.close()
out.close()
