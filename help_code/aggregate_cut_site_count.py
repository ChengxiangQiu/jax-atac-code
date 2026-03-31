import sys, os
import gzip

in_file = sys.argv[1]
out_file = sys.argv[2]

file = gzip.open(in_file)
out = open(out_file, "w")

curr_chr = 'chr1'
dat = {}

line = file.readline()
while line:
    l = line.decode("utf-8").rstrip().split('\t')
    if l[0] == curr_chr:
        dat[l[1]] = dat.get(l[1], 0)
        dat[l[1]] += 1
    else:
        for i in dat:
            out.write('\t'.join([curr_chr, i, str(int(i)+1), str(dat[i])]) + '\n')
        curr_chr = l[0]
        dat = {}
        dat[l[1]] = dat.get(l[1], 0)
        dat[l[1]] += 1
    line = file.readline()
file.close()

if dat:
    for i in dat:
        out.write('\t'.join([curr_chr, i, str(int(i)+1), str(dat[i])]) + '\n')

out.close()
