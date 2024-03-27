#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys

id2name = {'DMELA': 'Drosophila melanogaster', 'DSIMU': 'Drosophila simulans', 'DWILL': 'Drosophila willistoni', 'read2tree': 'Read2Tree'}
sps = ['DMELA', 'DSIMU', 'DWILL']
outfile = open('sum_contam_b.tsv', 'w')
outfile.write("species\tmethod\taverage completeness\taverage pident\n")
for sp in sps:
	for method in ['PhyloAln', 'read2tree']:
		com = 0
		pident = 0
		contam_name = 'contam'
		if method == 'PhyloAln':
			contam_name = 'contam_b'
		for line in open(sp + contam_name + '.' + method + '.tsv'):
			arr = line.rstrip("\n").split("\t")
			if arr[0] in ['file', 'total']:
				continue
			com += float(arr[4])
			pident += float(arr[5])
		com = int(com / 46 * 100 + 0.5) / 100
		pident = int(pident / 46 * 100 + 0.5) / 100
		outfile.write("{}\t{}\t{}\t{}\n".format(id2name.get(sp, sp), id2name.get(method, method), com, pident))
outfile.close()
outfile = open('sum_contam_b.PR.tsv', 'w')
outfile.write("species\ttype\torthogroup\tmapped reads\tTP\tFP\tFN\tTN\taccuracy\tprecision\trecall\tF-Measure\n")
for sp in sps:
	for filename in os.listdir('/public/lihaosen/PhyloAln/PhyloAln_contam_b/aa_out'):
		og = filename.replace('.fa', '')
		reads = {'clean': {}, 'outgroup': {}, 'cross': {}}
		for line in open('/public/lihaosen/PhyloAln/PhyloAln_contam_b/stat_info/' + og + '.clean_info.tsv'):
			arr = line.rstrip("\n").split("\t")
			if arr[0] == sp:
				for read in arr[1].split(','):
					if read == 'mapped':
						continue
					species = read.split('_', 1)[0]
					if reads['clean'].get(species) is None:
						reads['clean'][species] = 0
					reads['clean'][species] += 1
		for line in open('/public/lihaosen/PhyloAln/PhyloAln_contam_b/stat_info/' + og + '.cross_clean_info.tsv'):
			arr = line.rstrip("\n").split("\t")
			if arr[1] == sp:
				for read in arr[2].split(','):
					if read == 'mapped':
						continue
					species = read.split('_', 1)[0]
					if reads['cross'].get(species) is None:
						reads['cross'][species] = 0
					reads['cross'][species] += 1
		for line in open('/public/lihaosen/PhyloAln/PhyloAln_contam_b/stat_info/' + og + '_' + sp + '.out_clean_info.tsv'):
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
		A = 0
		P = 0
		R = 0
		F = 0
		if total['all'] > 0:
			A = (TP + TN) / total['all']
			if TP + FP > 0:
				P = TP / (TP + FP)
			if TP + FN > 0:
				R = TP / (TP + FN)
			if P + R > 0:
				F = 2 * P * R / (P + R)
		outfile.write("{}\tclean\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sp, og, total['all'], TP, FP, FN, TN, A, P, R, F))
		TP = total['outgroup'] - reads['outgroup'].get('DMELA', 0) - reads['outgroup'].get('DSIMU', 0) - reads['outgroup'].get('DWILL', 0)
		FP = reads['outgroup'].get('DMELA', 0) + reads['outgroup'].get('DSIMU', 0) + reads['outgroup'].get('DWILL', 0)
		FN = total['cross'] - reads['cross'].get('DMELA', 0) - reads['cross'].get('DSIMU', 0) - reads['cross'].get('DWILL', 0) + total['clean'] - reads['clean'].get('DMELA', 0) - reads['clean'].get('DSIMU', 0) - reads['clean'].get('DWILL', 0)
		TN = reads['cross'].get('DMELA', 0) + reads['cross'].get('DSIMU', 0) + reads['cross'].get('DWILL', 0) + reads['clean'].get('DMELA', 0) + reads['clean'].get('DSIMU', 0) + reads['clean'].get('DWILL', 0)
		A = 0
		P = 0
		R = 0
		F = 0
		if total['all'] > 0:
			A = (TP + TN) / total['all']
			if TP + FP > 0:
				P = TP / (TP + FP)
			if TP + FN > 0:
				R = TP / (TP + FN)
			if P + R > 0:
				F = 2 * P * R / (P + R)
		outfile.write("{}\toutgroup\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sp, og, total['all'], TP, FP, FN, TN, A, P, R, F))
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
		A = 0
		P = 0
		R = 0
		F = 0
		if TP + TN + FP + FN > 0:
			A = (TP + TN) / (TP + TN + FP + FN)
			if TP + FP > 0:
				P = TP / (TP + FP)
			if TP + FN > 0:
				R = TP / (TP + FN)
			if P + R > 0:
				F = 2 * P * R / (P + R)
		outfile.write("{}\tcross contamination\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sp, og, total['all'], TP, FP, FN, TN, A, P, R, F))
outfile.close()
