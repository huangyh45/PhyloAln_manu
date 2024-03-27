#!/usr/bin/env python3

import sys
import os
import shutil
import random

def get_reads(sp, num, out1, out2, ifc=''):
	reads1 = {}
	reads2 = {}
	seqid = None
	for line in open(sp + ifc + 'trans_10X_1.fq'):
		line = line.rstrip()
		if line.startswith('@'):
			seqid = line.replace('@', '', 1).replace('/1', '')
		elif seqid:
			reads1[seqid] = line
			seqid = None
	for line in open(sp + ifc + 'trans_10X_2.fq'):
		line = line.rstrip()
		if line.startswith('@'):
			seqid = line.replace('@', '', 1).replace('/2', '')
		elif seqid:
			reads2[seqid] = line
			seqid = None
	reads = []
	for seqid, seqstr in reads1.items():
		reads.append([seqstr, reads2[seqid]])
	new_reads = []
	for i in range(num):
		idnum = random.randint(0, len(reads)-1)
		out1.write(">{}_{}/1\n{}\n".format(sp, i+1, reads[idnum][0]))
		out2.write(">{}_{}/2\n{}\n".format(sp, i+1, reads[idnum][1]))
	print(sp + ': ' + str(num))

sps = ['DMELA', 'DSIMU', 'DWILL']
outsps = ['AAEGY', 'TCAST', 'HSAPI', 'DRERI', 'ATHAL', 'SCERE', 'ECOLI']

for sp in sps:
	outfile1 = open(sp + 'contam_1.fa', 'w')
	outfile2 = open(sp + 'contam_2.fa', 'w')
	get_reads(sp, random.randint(15000000, 20000000), outfile1, outfile2, 'c')
	print("Contamination:")
	for sp1 in sps:
		if sp == sp1:
			continue
		get_reads(sp1, random.randint(10000, 500000), outfile1, outfile2)
	for sp1 in outsps:
		get_reads(sp1, random.randint(10000, 500000), outfile1, outfile2)
	outfile1.close()
	outfile2.close()
	print("\n")
