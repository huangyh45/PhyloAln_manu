#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys

id2name = {'DMELA': 'Drosophila melanogaster', 'DSIMU': 'Drosophila simulans', 'DWILL': 'Drosophila willistoni', 'DBUSC': 'Drosophila busckii', 'DMOJA': 'Drosophila mojavensis', 'DPSEU': 'Drosophila pseudoobscura', 'DYAKU': 'Drosophila yakuba'}
sps = ['DMELA', 'DSIMU', 'DWILL']
outs = ['DBUSC', 'DMOJA', 'DPSEU', 'DYAKU']
removes = {'DMELA': ['OG0003615.fa', 'OG0003803.fa', 'OG0003845.fa', 'OG0004076.fa', 'OG0004221.fa'], 'DSIMU': ['OG0004212.fa', 'OG0003977.fa', 'OG0003820.fa', 'OG0003709.fa', 'OG0003531.fa'], 'DWILL': ['OG0003636.fa', 'OG0003729.fa', 'OG0003830.fa', 'OG0003956.fa', 'OG0004227.fa']}
outfile = open('sum_contam_out.tsv', 'w')
outfile.write("species\toutgroup\tassembly\taverage completeness\taverage pident\n")
for sp in sps:
	for out in outs:
		com = 0
		pident = 0
		for line in open(sp + 'contam.' + out + '.PhyloAln.tsv'):
			arr = line.rstrip("\n").split("\t")
			if arr[0] in ['file', 'total'] or arr[0] in removes[sp]:
				continue
			com += float(arr[4])
			pident += float(arr[5])
		com = int(com / 41 * 100 + 0.5) / 100
		pident = int(pident / 41 * 100 + 0.5) / 100
		outfile.write("{}\t{}\tyes\t{}\t{}\n".format(id2name.get(sp, sp), id2name.get(out, out), com, pident))
		com = 0
		pident = 0
		for line in open(sp + 'contam_b.' + out + '.PhyloAln.tsv'):
			arr = line.rstrip("\n").split("\t")
			if arr[0] in ['file', 'total'] or arr[0] in removes[sp]:
				continue
			com += float(arr[4])
			pident += float(arr[5])
		com = int(com / 41 * 100 + 0.5) / 100
		pident = int(pident / 41 * 100 + 0.5) / 100
		outfile.write("{}\t{}\tno\t{}\t{}\n".format(id2name.get(sp, sp), id2name.get(out, out), com, pident))
outfile.close()
for out in outs:
	for ifb in ['', '_b']:
		outfile = open(f"sum_contam{ifb}.{out}.PR.tsv", 'w')
		outfile.write("species\ttype\torthogroup\tremove\tmapped reads\tTP\tFP\tFN\tTN\taccuracy\tprecision\trecall\tF-Measure\n")
		for sp in sps:
			for filename in os.listdir(f"/public/lihaosen/PhyloAln/PhyloAln_contam{ifb}_{out}/aa_out"):
				og = filename.replace('.fa', '')
				ifr = 'no'
				if filename in removes[sp]:
					ifr = 'yes'
				reads = {'clean': {}, 'outgroup': {}, 'cross': {}}
				for line in open(f"/public/lihaosen/PhyloAln/PhyloAln_contam{ifb}_{out}/stat_info/" + og + '.clean_info.tsv'):
					arr = line.rstrip("\n").split("\t")
					if arr[0] == sp:
						for read in arr[1].split(','):
							if read == 'mapped':
								continue
							species = read.split('_', 1)[0]
							if reads['clean'].get(species) is None:
								reads['clean'][species] = 0
							reads['clean'][species] += 1
				for line in open(f"/public/lihaosen/PhyloAln/PhyloAln_contam{ifb}_{out}/stat_info/" + og + '.cross_clean_info.tsv'):
					arr = line.rstrip("\n").split("\t")
					if arr[1] == sp:
						for read in arr[2].split(','):
							if read == 'mapped':
								continue
							species = read.split('_', 1)[0]
							if reads['cross'].get(species) is None:
								reads['cross'][species] = 0
							reads['cross'][species] += 1
				for line in open(f"/public/lihaosen/PhyloAln/PhyloAln_contam{ifb}_{out}/stat_info/" + og + '_' + sp + '.out_clean_info.tsv'):
					arr = line.rstrip("\n").split("\t")
					for read in arr[5].split(','):
						if read == 'mapped':
							continue
						species = read.split('_', 1)[0]
						if reads['outgroup'].get(species) is None:
							reads['outgroup'][species] = 0
						reads['outgroup'][species] += 1
				total = {'clean': 0, 'outgroup': 0, 'cross': 0, 'sp': 0, 'all': 0}
				for contam_type, sp2count in reads.items():
					total['sp'] += sp2count.get(sp, 0)
					for species, count in sp2count.items():
						total['all'] += count
						total[contam_type] += count
				TP = reads['clean'].get(sp, 0)
				FP = total['clean'] - reads['clean'].get(sp, 0)
				FN = reads['outgroup'].get(sp, 0)
				TN = total['outgroup'] - reads['outgroup'].get(sp, 0) + total['cross'] - reads['cross'].get(sp, 0)
				A = 'invalid'
				P = 'invalid'
				R = 'invalid'
				F = 'invalid'
				if total['all'] > 0:
					A = (TP + TN) / total['all']
					if TP + FP > 0:
						P = TP / (TP + FP)
					if TP + FN > 0:
						R = TP / (TP + FN)
					if str(P) != 'invalid' and str(R) != 'invalid':
						if P + R > 0:
							F = 2 * P * R / (P + R)
				outfile.write("{}\tclean\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sp, og, ifr, total['all'], TP, FP, FN, TN, A, P, R, F))
				TP = total['outgroup'] - reads['outgroup'].get('DMELA', 0) - reads['outgroup'].get('DSIMU', 0) - reads['outgroup'].get('DWILL', 0)
				FP = reads['outgroup'].get('DMELA', 0) + reads['outgroup'].get('DSIMU', 0) + reads['outgroup'].get('DWILL', 0)
				FN = total['cross'] - reads['cross'].get('DMELA', 0) - reads['cross'].get('DSIMU', 0) - reads['cross'].get('DWILL', 0) + total['clean'] - reads['clean'].get('DMELA', 0) - reads['clean'].get('DSIMU', 0) - reads['clean'].get('DWILL', 0)
				TN = reads['cross'].get('DMELA', 0) + reads['cross'].get('DSIMU', 0) + reads['cross'].get('DWILL', 0) + reads['clean'].get('DMELA', 0) + reads['clean'].get('DSIMU', 0) + reads['clean'].get('DWILL', 0)
				A = 'invalid'
				P = 'invalid'
				R = 'invalid'
				F = 'invalid'
				if total['all'] > 0:
					A = (TP + TN) / total['all']
					if TP + FP > 0:
						P = TP / (TP + FP)
					if TP + FN > 0:
						R = TP / (TP + FN)
					if str(P) != 'invalid' and str(R) != 'invalid':
						if P + R > 0:
							F = 2 * P * R / (P + R)
				outfile.write("{}\toutgroup\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sp, og, ifr, total['all'], TP, FP, FN, TN, A, P, R, F))
				TP = 0
				FP = 0
				FN = 0
				TN = 0
				for s in sps:
					if s == sp:
						FP += reads['cross'].get(s, 0)
						TN += reads['clean'].get(s, 0)
					else:
						TP += reads['cross'].get(s, 0)
						FN += reads['clean'].get(s, 0)
				A = 'invalid'
				P = 'invalid'
				R = 'invalid'
				F = 'invalid'
				if TP + TN + FP + FN > 0:
					A = (TP + TN) / (TP + TN + FP + FN)
					if TP + FP > 0:
						P = TP / (TP + FP)
					if TP + FN > 0:
						R = TP / (TP + FN)
					if str(P) != 'invalid' and str(R) != 'invalid':
						if P + R > 0:
							F = 2 * P * R / (P + R)
				outfile.write("{}\tcross contamination\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sp, og, ifr, total['all'], TP, FP, FN, TN, A, P, R, F))
		outfile.close()
