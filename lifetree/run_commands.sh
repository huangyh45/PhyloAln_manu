# run OrthoFinder and generate the CDS alignments via custom scripts
time ./orthogroup.pl -g genome -f gff -p 20 -b -e -t -k
cd ortho
time ./orthomsa.pl -d OrthoFinder -p 20 -c cds -k -u ECOLI -y OrthoFinder/Orthogroups/Orthogroups_SingleCopyOrthologues.txt -m -a -e
cd ..

# remove the target sequences
PhyloAln/scripts/select_seqs.py ortho/OrthoFinder/CDS_trim ECOLI,SCERE,DRERI,TCAST,AAEGY ortho/ortho/OGs_treelife
PhyloAln/scripts/select_seqs.py ortho/OrthoFinder/Orthogroup_Sequences ECOLI,SCERE,DRERI,TCAST,AAEGY orthogroup/ortho/OGs_treelife_unaln # only retain 46 single-copy gene alignments

# simulated the reads of genomic and transcript sequences
for file in genome/$sp.fna
do
	sp=`basename $file .fna`
	cp genome/$sp.fna dataset/$sp.genome.fasta
	gffread -g genome/$sp.fna -w dataset/$sp.transcript.fasta gff/$sp.gff
done
cd dataset
for file in $sp.transcript.fasta
do
	sp=`basename $file .transcript.fasta`
	covs=(2 5 10 20 30)
	for i in "${covs[@]}"
	do
		art_illumina -ss HS25 -i $sp.genome.fasta -p -l 100 -f $i -m 200 -s 10 -o ${sp}_${i}X
		art_illumina -ss HS25 -i $sp.transcript.fasta -p -l 100 -f $i -m 200 -s 10 -o ${sp}trans_${i}X
		./readsim_multi.py $sp.genome.fasta ${sp}_pacbio_${i}X pacbio $i 20
		./readsim_multi.py $sp.transcript.fasta ${sp}trans_pacbio_${i}X pacbio $i 20
		./readsim_multi.py $sp.genome.fasta ${sp}_nanopore_${i}X nanopore $i 20
		./readsim_multi.py $sp.transcript.fasta ${sp}trans_nanopore_${i}X nanopore $i 20
	done
done
cd ..

# run PhyloAln for Illumina reads
covs=(2 5 10 20 30)
for i in "${covs[@]}"
do
	time PhyloAln/PhyloAln -d orthogroup/ortho/OGs_treelife -s ATHAL_${i}X -i dataset/ATHAL_${i}X_1.fq dataset/ATHAL_${i}X_2.fq -f fastq -o PhyloAln_ATHAL_${i}X -p 20 -m codon -u ECOLI
	time PhyloAln/PhyloAln -d orthogroup/ortho/OGs_treelife -s ATHALtrans_${i}X -i dataset/ATHALtrans_${i}X_1.fq dataset/ATHALtrans_${i}X_2.fq -f fastq -o PhyloAln_ATHALtrans_${i}X -p 20 -m codon -u ECOLI
	time PhyloAln/PhyloAln -d orthogroup/ortho/OGs_treelife -s DMELA_${i}X -i dataset/DMELA_${i}X_1.fq dataset/DMELA_${i}X_2.fq -f fastq -o PhyloAln_DMELA_${i}X -p 20 -m codon -u ECOLI
	time PhyloAln/PhyloAln -d orthogroup/ortho/OGs_treelife -s DMELAtrans_${i}X -i dataset/DMELAtrans_${i}X_1.fq dataset/DMELAtrans_${i}X_2.fq -f fastq -o PhyloAln_DMELAtrans_${i}X -p 20 -m codon -u ECOLI
	time PhyloAln/PhyloAln -d orthogroup/ortho/OGs_treelife -s HSAPItrans_${i}X -i dataset/HSAPItrans_${i}X_1.fq dataset/HSAPItrans_${i}X_2.fq -f fastq -o PhyloAln_HSAPItrans_${i}X -p 20 -m codon -u ECOLI
done
covs=(2 5 10)
for i in "${covs[@]}"
do
	time PhyloAln/PhyloAln -d orthogroup/ortho/OGs_treelife -s HSAPI_${i}X -i dataset/HSAPI_${i}X_1.fq dataset/HSAPI_${i}X_2.fq -f fastq -o PhyloAln_HSAPI_${i}X -p 20 -m codon -u ECOLI
done
covs=(20 30)
for i in "${covs[@]}"
do
	time PhyloAln/PhyloAln -d orthogroup/ortho/OGs_treelife -s HSAPI_${i}X -i dataset/HSAPI_${i}X_1.fq dataset/HSAPI_${i}X_2.fq -f fastq -o PhyloAln_HSAPI_${i}X -p 20 -m codon -u ECOLI --low_mem
done

# run PhyloAln for genomic long reads
covs=(2 5 10 20 30)
for i in "${covs[@]}"
do
	time PhyloAln/PhyloAln -d orthogroup/ortho/OGs_treelife -s ATHAL_pacbio_${i}X -i dataset/ATHAL_pacbio_${i}X.fasta -f fasta -o PhyloAln_ATHAL_pacbio_${i}X -p 20 -m dna_codon -u ECOLI -l 200
	time PhyloAln/PhyloAln -d orthogroup/ortho/OGs_treelife -s ATHAL_nanopore_${i}X -i dataset/ATHAL_nanopore_${i}X.fasta -f fasta -o PhyloAln_ATHAL_nanopore_${i}X -p 20 -m dna_codon -u ECOLI -l 200
done
covs=(2 5 10 20 30)
for i in "${covs[@]}"
do
	time PhyloAln/PhyloAln -d orthogroup/ortho/OGs_treelife -s DMELA_pacbio_${i}X -i dataset/DMELA_pacbio_${i}X.fasta -f fasta -o PhyloAln_DMELA_pacbio_${i}X -p 20 -m dna_codon -u ECOLI -l 200
	time PhyloAln/PhyloAln -d orthogroup/ortho/OGs_treelife -s DMELA_nanopore_${i}X -i dataset/DMELA_nanopore_${i}X.fasta -f fasta -o PhyloAln_DMELA_nanopore_${i}X -p 20 -m dna_codon -u ECOLI -l 200
done
covs=(2 5 10)
for i in "${covs[@]}"
do
	time PhyloAln/PhyloAln -d orthogroup/ortho/OGs_treelife -s HSAPI_pacbio_${i}X -i dataset/HSAPI_pacbio_${i}X.fasta -f fasta -o PhyloAln_HSAPI_pacbio_${i}X -p 20 -m dna_codon -u ECOLI -l 200
	time PhyloAln/PhyloAln -d orthogroup/ortho/OGs_treelife -s HSAPI_nanopore_${i}X -i dataset/HSAPI_nanopore_${i}X.fasta -f fasta -o PhyloAln_HSAPI_nanopore_${i}X -p 20 -m dna_codon -u ECOLI -l 200
done
covs=(20 30)
for i in "${covs[@]}"
do
	time PhyloAln/PhyloAln -d orthogroup/ortho/OGs_treelife -s HSAPI_pacbio_${i}X -i dataset/HSAPI_pacbio_${i}X.fasta -f fasta -o PhyloAln_HSAPI_pacbio_${i}X -p 20 -m dna_codon -u ECOLI -l 200 --low_mem
	time PhyloAln/PhyloAln -d orthogroup/ortho/OGs_treelife -s HSAPI_nanopore_${i}X -i dataset/HSAPI_nanopore_${i}X.fasta -f fasta -o PhyloAln_HSAPI_nanopore_${i}X -p 20 -m dna_codon -u ECOLI -l 200 --low_mem
done

# run PhyloAln for transcriptomic long reads
covs=(2 5 10 20 30)
for i in "${covs[@]}"
do
	time PhyloAln/PhyloAln -d orthogroup/ortho/OGs_treelife -s ATHALtrans_pacbio_${i}X -i dataset/ATHALtrans_pacbio_${i}X.fasta -f fasta -o PhyloAln_ATHALtrans_pacbio_${i}X -p 20 -m dna_codon -u ECOLI
	time PhyloAln/PhyloAln -d orthogroup/ortho/OGs_treelife -s ATHALtrans_nanopore_${i}X -i dataset/ATHALtrans_nanopore_${i}X.fasta -f fasta -o PhyloAln_ATHALtrans_nanopore_${i}X -p 20 -m dna_codon -u ECOLI
done
covs=(2 5 10 20 30)
for i in "${covs[@]}"
do
	time PhyloAln/PhyloAln -d orthogroup/ortho/OGs_treelife -s DMELAtrans_pacbio_${i}X -i dataset/DMELAtrans_pacbio_${i}X.fasta -f fasta -o PhyloAln_DMELAtrans_pacbio_${i}X -p 20 -m dna_codon -u ECOLI
	time PhyloAln/PhyloAln -d orthogroup/ortho/OGs_treelife -s DMELAtrans_nanopore_${i}X -i dataset/DMELAtrans_nanopore_${i}X.fasta -f fasta -o PhyloAln_DMELAtrans_nanopore_${i}X -p 20 -m dna_codon -u ECOLI
done
covs=(2 5 10 20 30)
for i in "${covs[@]}"
do
	time PhyloAln/PhyloAln -d orthogroup/ortho/OGs_treelife -s HSAPItrans_pacbio_${i}X -i dataset/HSAPItrans_pacbio_${i}X.fasta -f fasta -o PhyloAln_HSAPItrans_pacbio_${i}X -p 20 -m dna_codon -u ECOLI
	time PhyloAln/PhyloAln -d orthogroup/ortho/OGs_treelife -s HSAPItrans_nanopore_${i}X -i dataset/HSAPItrans_nanopore_${i}X.fasta -f fasta -o PhyloAln_HSAPItrans_nanopore_${i}X -p 20 -m dna_codon -u ECOLI
done

# run PhyloAln for assemblies
time PhyloAln/PhyloAln -d orthogroup/ortho/OGs_treelife -s ATHAL_seq -i dataset/ATHAL.genome.fasta -f large_fasta -o PhyloAln_ATHAL_seq -p 20 -m codon -u ECOLI -b -r -l 200
time PhyloAln/PhyloAln -d orthogroup/ortho/OGs_treelife -s ATHALtrans_seq -i dataset/ATHAL.transcript.fasta -f fasta -o PhyloAln_ATHALtrans_seq -p 20 -m codon -u ECOLI -b -r
time PhyloAln/PhyloAln -d orthogroup/ortho/OGs_treelife -s DMELA_seq -i dataset/DMELA.genome.fasta -f large_fasta -o PhyloAln_DMELA_seq -p 20 -m codon -u ECOLI -b -r -l 200
time PhyloAln/PhyloAln -d orthogroup/ortho/OGs_treelife -s DMELAtrans_seq -i dataset/DMELA.transcript.fasta -f fasta -o PhyloAln_DMELAtrans_seq -p 20 -m codon -u ECOLI -b -r
time PhyloAln/PhyloAln -d orthogroup/ortho/OGs_treelife -s HSAPI_seq -i dataset/HSAPI.genome.fasta -f large_fasta -o PhyloAln_HSAPI_seq -p 20 -m codon -u ECOLI -b -r -l 200 --low_mem
time PhyloAln/PhyloAln -d orthogroup/ortho/OGs_treelife -s HSAPItrans_seq -i dataset/HSAPI.transcript.fasta -f fasta -o PhyloAln_HSAPItrans_seq -p 20 -m codon -u ECOLI -b -r

# calculate the percent completeness and identity of PhyloAln results
PhyloAln/scripts/test_effect.py orthogroup/ortho/OrthoFinder/CDS_trim:ATHAL PhyloAln_ATHAL_seq/nt_out:ATHAL_seq ATHAL_seq.PhyloAln.tsv N . .fa ECOLI,SCERE,DRERI,TCAST,AAEGY
PhyloAln/scripts/test_effect.py orthogroup/ortho/OrthoFinder/CDS_trim:HSAPI PhyloAln_HSAPI_seq/nt_out:HSAPI_seq HSAPI_seq.PhyloAln.tsv N . .fa ECOLI,SCERE,DRERI,TCAST,AAEGY
PhyloAln/scripts/test_effect.py orthogroup/ortho/OrthoFinder/CDS_trim:DMELA PhyloAln_DMELA_seq/nt_out:DMELA_seq DMELA_seq.PhyloAln.tsv N . .fa ECOLI,SCERE,DRERI,TCAST,AAEGY
PhyloAln/scripts/test_effect.py orthogroup/ortho/OrthoFinder/CDS_trim:ATHAL PhyloAln_ATHALtrans_seq/nt_out:ATHALtrans_seq ATHALtrans_seq.PhyloAln.tsv N . .fa ECOLI,SCERE,DRERI,TCAST,AAEGY
PhyloAln/scripts/test_effect.py orthogroup/ortho/OrthoFinder/CDS_trim:HSAPI PhyloAln_HSAPItrans_seq/nt_out:HSAPItrans_seq HSAPItrans_seq.PhyloAln.tsv N . .fa ECOLI,SCERE,DRERI,TCAST,AAEGY
PhyloAln/scripts/test_effect.py orthogroup/ortho/OrthoFinder/CDS_trim:DMELA PhyloAln_DMELAtrans_seq/nt_out:DMELAtrans_seq DMELAtrans_seq.PhyloAln.tsv N . .fa ECOLI,SCERE,DRERI,TCAST,AAEGY
covs=(2 5 10 20 30)
for i in "${covs[@]}"
do
	PhyloAln/scripts/test_effect.py orthogroup/ortho/OrthoFinder/CDS_trim:ATHAL PhyloAln_ATHAL_${i}X/nt_out:ATHAL_${i}X ATHAL_${i}X.PhyloAln.tsv N . .fa ECOLI,SCERE,DRERI,TCAST,AAEGY
	PhyloAln/scripts/test_effect.py orthogroup/ortho/OrthoFinder/CDS_trim:HSAPI PhyloAln_HSAPI_${i}X/nt_out:HSAPI_${i}X HSAPI_${i}X.PhyloAln.tsv N . .fa ECOLI,SCERE,DRERI,TCAST,AAEGY
	PhyloAln/scripts/test_effect.py orthogroup/ortho/OrthoFinder/CDS_trim:DMELA PhyloAln_DMELA_${i}X/nt_out:DMELA_${i}X DMELA_${i}X.PhyloAln.tsv N . .fa ECOLI,SCERE,DRERI,TCAST,AAEGY
	PhyloAln/scripts/test_effect.py orthogroup/ortho/OrthoFinder/CDS_trim:ATHAL PhyloAln_ATHALtrans_${i}X/nt_out:ATHALtrans_${i}X ATHALtrans_${i}X.PhyloAln.tsv N . .fa ECOLI,SCERE,DRERI,TCAST,AAEGY
	PhyloAln/scripts/test_effect.py orthogroup/ortho/OrthoFinder/CDS_trim:HSAPI PhyloAln_HSAPItrans_${i}X/nt_out:HSAPItrans_${i}X HSAPItrans_${i}X.PhyloAln.tsv N . .fa ECOLI,SCERE,DRERI,TCAST,AAEGY
	PhyloAln/scripts/test_effect.py orthogroup/ortho/OrthoFinder/CDS_trim:DMELA PhyloAln_DMELAtrans_${i}X/nt_out:DMELAtrans_${i}X DMELAtrans_${i}X.PhyloAln.tsv N . .fa ECOLI,SCERE,DRERI,TCAST,AAEGY
	PhyloAln/scripts/test_effect.py orthogroup/ortho/OrthoFinder/CDS_trim:ATHAL PhyloAln_ATHAL_pacbio_${i}X/nt_out:ATHAL_pacbio_${i}X ATHAL_pacbio_${i}X.PhyloAln.tsv N . .fa ECOLI,SCERE,DRERI,TCAST,AAEGY
	PhyloAln/scripts/test_effect.py orthogroup/ortho/OrthoFinder/CDS_trim:HSAPI PhyloAln_HSAPI_pacbio_${i}X/nt_out:HSAPI_pacbio_${i}X HSAPI_pacbio_${i}X.PhyloAln.tsv N . .fa ECOLI,SCERE,DRERI,TCAST,AAEGY
	PhyloAln/scripts/test_effect.py orthogroup/ortho/OrthoFinder/CDS_trim:DMELA PhyloAln_DMELA_pacbio_${i}X/nt_out:DMELA_pacbio_${i}X DMELA_pacbio_${i}X.PhyloAln.tsv N . .fa ECOLI,SCERE,DRERI,TCAST,AAEGY
	PhyloAln/scripts/test_effect.py orthogroup/ortho/OrthoFinder/CDS_trim:ATHAL PhyloAln_ATHALtrans_pacbio_${i}X/nt_out:ATHALtrans_pacbio_${i}X ATHALtrans_pacbio_${i}X.PhyloAln.tsv N . .fa ECOLI,SCERE,DRERI,TCAST,AAEGY
	PhyloAln/scripts/test_effect.py orthogroup/ortho/OrthoFinder/CDS_trim:HSAPI PhyloAln_HSAPItrans_pacbio_${i}X/nt_out:HSAPItrans_pacbio_${i}X HSAPItrans_pacbio_${i}X.PhyloAln.tsv N . .fa ECOLI,SCERE,DRERI,TCAST,AAEGY
	PhyloAln/scripts/test_effect.py orthogroup/ortho/OrthoFinder/CDS_trim:DMELA PhyloAln_DMELAtrans_pacbio_${i}X/nt_out:DMELAtrans_pacbio_${i}X DMELAtrans_pacbio_${i}X.PhyloAln.tsv N . .fa ECOLI,SCERE,DRERI,TCAST,AAEGY
	PhyloAln/scripts/test_effect.py orthogroup/ortho/OrthoFinder/CDS_trim:ATHAL PhyloAln_ATHAL_nanopore_${i}X/nt_out:ATHAL_nanopore_${i}X ATHAL_nanopore_${i}X.PhyloAln.tsv N . .fa ECOLI,SCERE,DRERI,TCAST,AAEGY
	PhyloAln/scripts/test_effect.py orthogroup/ortho/OrthoFinder/CDS_trim:HSAPI PhyloAln_HSAPI_nanopore_${i}X/nt_out:HSAPI_nanopore_${i}X HSAPI_nanopore_${i}X.PhyloAln.tsv N . .fa ECOLI,SCERE,DRERI,TCAST,AAEGY
	PhyloAln/scripts/test_effect.py orthogroup/ortho/OrthoFinder/CDS_trim:DMELA PhyloAln_DMELA_nanopore_${i}X/nt_out:DMELA_nanopore_${i}X DMELA_nanopore_${i}X.PhyloAln.tsv N . .fa ECOLI,SCERE,DRERI,TCAST,AAEGY
	PhyloAln/scripts/test_effect.py orthogroup/ortho/OrthoFinder/CDS_trim:ATHAL PhyloAln_ATHALtrans_nanopore_${i}X/nt_out:ATHALtrans_nanopore_${i}X ATHALtrans_nanopore_${i}X.PhyloAln.tsv N . .fa ECOLI,SCERE,DRERI,TCAST,AAEGY
	PhyloAln/scripts/test_effect.py orthogroup/ortho/OrthoFinder/CDS_trim:HSAPI PhyloAln_HSAPItrans_nanopore_${i}X/nt_out:HSAPItrans_nanopore_${i}X HSAPItrans_nanopore_${i}X.PhyloAln.tsv N . .fa ECOLI,SCERE,DRERI,TCAST,AAEGY
	PhyloAln/scripts/test_effect.py orthogroup/ortho/OrthoFinder/CDS_trim:DMELA PhyloAln_DMELAtrans_nanopore_${i}X/nt_out:DMELAtrans_nanopore_${i}X DMELAtrans_nanopore_${i}X.PhyloAln.tsv N . .fa ECOLI,SCERE,DRERI,TCAST,AAEGY
done

# run Read2Tree for Illumina reads
time read2tree --standalone_path orthogroup/ortho/OGs_treelife_unaln --output_path read2tree_lifetree --reference --dna_reference orthogroup/all.cds.fasta
covs=(2 5 10 20 30)
for i in "${covs[@]}"
do
	time read2tree --output_path read2tree_lifetree --standalone_path orthogroup/ortho/OGs_treelife_unaln --reads dataset/ATHALtrans_${i}X_1.fq dataset/ATHALtrans_${i}X_2.fq --threads 20 --dna_reference orthogroup/all.cds.fasta -s ATHALtrans_${i}X
	time read2tree --output_path read2tree_lifetree --standalone_path orthogroup/ortho/OGs_treelife_unaln --reads dataset/ATHAL_${i}X_1.fq dataset/ATHAL_${i}X_2.fq --threads 20 --dna_reference orthogroup/all.cds.fasta -s ATHAL_${i}X
	time read2tree --output_path read2tree_lifetree --standalone_path orthogroup/ortho/OGs_treelife_unaln --reads dataset/DMELAtrans_${i}X_1.fq dataset/DMELAtrans_${i}X_2.fq --threads 20 --dna_reference orthogroup/all.cds.fasta -s DMELAtrans_${i}X
	time read2tree --output_path read2tree_lifetree --standalone_path orthogroup/ortho/OGs_treelife_unaln --reads dataset/DMELA_${i}X_1.fq dataset/DMELA_${i}X_2.fq --threads 20 --dna_reference orthogroup/all.cds.fasta -s DMELA_${i}X
	time read2tree --output_path read2tree_lifetree --standalone_path orthogroup/ortho/OGs_treelife_unaln --reads dataset/HSAPItrans_${i}X_1.fq dataset/HSAPItrans_${i}X_2.fq --threads 20 --dna_reference orthogroup/all.cds.fasta -s HSAPItrans_${i}X
	time read2tree --output_path read2tree_lifetree --standalone_path orthogroup/ortho/OGs_treelife_unaln --reads dataset/HSAPI_${i}X_1.fq dataset/HSAPI_${i}X_2.fq --threads 20 --dna_reference orthogroup/all.cds.fasta -s HSAPI_${i}X
done

# run Read2Tree for long reads, fa2fq is a tool in IDBA
covs=(2 5 10 20 30)
for i in "${covs[@]}"
do
	fa2fq dataset/ATHALtrans_pacbio_${i}X.fasta dataset/ATHALtrans_pacbio_${i}X.fastq
	time read2tree --output_path read2tree_lifetree --standalone_path orthogroup/ortho/OGs_treelife_unaln --reads dataset/ATHALtrans_pacbio_${i}X.fastq --read_type long --threads 20 --dna_reference orthogroup/all.cds.fasta -s ATHALtrans_pacbio_${i}X --split_reads --ngmlr_parameters pacbio,256,0.25
	fa2fq dataset/ATHALtrans_nanopore_${i}X.fasta dataset/ATHALtrans_nanopore_${i}X.fastq
	time read2tree --output_path read2tree_lifetree --standalone_path orthogroup/ortho/OGs_treelife_unaln --reads dataset/ATHALtrans_nanopore_${i}X.fastq --read_type long --threads 20 --dna_reference orthogroup/all.cds.fasta -s ATHALtrans_nanopore_${i}X --split_reads
	fa2fq dataset/ATHAL_pacbio_${i}X.fasta dataset/ATHAL_pacbio_${i}X.fastq
	time read2tree --output_path read2tree_lifetree --standalone_path orthogroup/ortho/OGs_treelife_unaln --reads dataset/ATHAL_pacbio_${i}X.fastq --read_type long --threads 20 --dna_reference orthogroup/all.cds.fasta -s ATHAL_pacbio_${i}X --split_reads --ngmlr_parameters pacbio,256,0.25
	fa2fq dataset/ATHAL_nanopore_${i}X.fasta dataset/ATHAL_nanopore_${i}X.fastq
	time read2tree --output_path read2tree_lifetree --standalone_path orthogroup/ortho/OGs_treelife_unaln --reads dataset/ATHAL_nanopore_${i}X.fastq --read_type long --threads 20 --dna_reference orthogroup/all.cds.fasta -s ATHAL_nanopore_${i}X --split_reads
	fa2fq dataset/DMELAtrans_pacbio_${i}X.fasta dataset/DMELAtrans_pacbio_${i}X.fastq
	time read2tree --output_path read2tree_lifetree --standalone_path orthogroup/ortho/OGs_treelife_unaln --reads dataset/DMELAtrans_pacbio_${i}X.fastq --read_type long --threads 20 --dna_reference orthogroup/all.cds.fasta -s DMELAtrans_pacbio_${i}X --split_reads --ngmlr_parameters pacbio,256,0.25
	fa2fq dataset/DMELAtrans_nanopore_${i}X.fasta dataset/DMELAtrans_nanopore_${i}X.fastq
	time read2tree --output_path read2tree_lifetree --standalone_path orthogroup/ortho/OGs_treelife_unaln --reads dataset/DMELAtrans_nanopore_${i}X.fastq --read_type long --threads 20 --dna_reference orthogroup/all.cds.fasta -s DMELAtrans_nanopore_${i}X --split_reads
	fa2fq dataset/DMELA_pacbio_${i}X.fasta dataset/DMELA_pacbio_${i}X.fastq
	time read2tree --output_path read2tree_lifetree --standalone_path orthogroup/ortho/OGs_treelife_unaln --reads dataset/DMELA_pacbio_${i}X.fastq --read_type long --threads 20 --dna_reference orthogroup/all.cds.fasta -s DMELA_pacbio_${i}X --split_reads --ngmlr_parameters pacbio,256,0.25
	fa2fq dataset/DMELA_nanopore_${i}X.fasta dataset/DMELA_nanopore_${i}X.fastq
	time read2tree --output_path read2tree_lifetree --standalone_path orthogroup/ortho/OGs_treelife_unaln --reads dataset/DMELA_nanopore_${i}X.fastq --read_type long --threads 20 --dna_reference orthogroup/all.cds.fasta -s DMELA_nanopore_${i}X --split_reads
	fa2fq dataset/HSAPItrans_pacbio_${i}X.fasta dataset/HSAPItrans_pacbio_${i}X.fastq
	time read2tree --output_path read2tree_lifetree --standalone_path orthogroup/ortho/OGs_treelife_unaln --reads dataset/HSAPItrans_pacbio_${i}X.fastq --read_type long --threads 20 --dna_reference orthogroup/all.cds.fasta -s HSAPItrans_pacbio_${i}X --split_reads --ngmlr_parameters pacbio,256,0.25
	fa2fq dataset/HSAPItrans_nanopore_${i}X.fasta dataset/HSAPItrans_nanopore_${i}X.fastq
	time read2tree --output_path read2tree_lifetree --standalone_path orthogroup/ortho/OGs_treelife_unaln --reads dataset/HSAPItrans_nanopore_${i}X.fastq --read_type long --threads 20 --dna_reference orthogroup/all.cds.fasta -s HSAPItrans_nanopore_${i}X --split_reads
	fa2fq dataset/HSAPI_pacbio_${i}X.fasta dataset/HSAPI_pacbio_${i}X.fastq
	time read2tree --output_path read2tree_lifetree --standalone_path orthogroup/ortho/OGs_treelife_unaln --reads dataset/HSAPI_pacbio_${i}X.fastq --read_type long --threads 20 --dna_reference orthogroup/all.cds.fasta -s HSAPI_pacbio_${i}X --split_reads --ngmlr_parameters pacbio,256,0.25
	fa2fq dataset/HSAPI_nanopore_${i}X.fasta dataset/HSAPI_nanopore_${i}X.fastq
	time read2tree --output_path read2tree_lifetree --standalone_path orthogroup/ortho/OGs_treelife_unaln --reads dataset/HSAPI_nanopore_${i}X.fastq --read_type long --threads 20 --dna_reference orthogroup/all.cds.fasta -s HSAPI_nanopore_${i}X --split_reads
done

# calculate the percent completeness and identity of Read2Tree results
for file in read2tree_lifetree/06_align_merge_dna/*.fa
do
	name=`basename $file`
	readal -in $file -out read2tree_result/$name -fasta
done
PhyloAln/scripts/select_seqs.py ortho/OrthoFinder/Orthogroup_Sequences ATHAL,DMELA,HSAPI ogs
for file in ortho/single_copy_ref/*.fa
do
	name=`basename $file .fa`
	readal -in read2tree_lifetree/03_align_aa/$name.phy -out read2tree_test/$name.ref_aa.fas -fasta
	mafft --add ogs/$name.fa --localpair --maxiterate 1000 --thread 10 read2tree_test/$name.ref_aa.fas > read2tree_test/$name.aa_aln.fas
	readal -in read2tree_lifetree/03_align_dna/$name.phy -out read2tree_test/$name.dna.fas -fasta
	sed -i 's/-//g' $name.dna.fas
	cat ortho/OrthoFinder/CDS_seq/$name.fa >> read2tree_test/$name.dna.fas
	trimal -in read2tree_test/$name.aa_aln.fas -out read2tree_test/$name.fa -keepheader -backtrans read2tree_test/$name.dna.fas
done
for file in read2tree_lifetree/05_ogs_map_*_dna
do
	sp=`basename $file _dna | sed 's/^05_ogs_map_//'`
	sp0=`echo $sp | sed 's/_[^_]*//g'`
	PhyloAln/scripts/test_effect.py read2tree_test:$sp0 read2tree_result:$sp $sp.read2tree.tsv N . .fa ECOLI,SCERE,DRERI,TCAST,AAEGY
done
