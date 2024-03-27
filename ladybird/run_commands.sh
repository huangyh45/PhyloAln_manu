# assemble the transcriptomes via a custom script for 00-04 batches respectively
time ./transassemble.pl -n 20 -p 1 -T -C 1 -R 2 -E

# run OrthoFinder and generate the CDS alignments via custom scripts
time ./orthogroup.pl -g genome -f gff -y .gff3 -a trans_aa -c trans_cds -p 20 -b -e -t -k
cd ortho
time ./orthomsa.pl -d OrthoFinder -p 20 -c cds -k -u DHELO -y singlecopy.list -m -a -e
cd ..

# remove the target sequences
PhyloAln/scripts/select_seqs.py ortho/OrthoFinder/CDS_trim ABIPU,CIMPU,CMONT,CSEPT,DHELO,HAXYR,HSEDE,HVIGI,HVIMA,MDISC,NPUMI,PJAPO ortho/genome_ortho/CDS_trim
PhyloAln/scripts/select_seqs.py ortho/OrthoFinder/Orthogroup_Sequences ABIPU,CIMPU,CMONT,CSEPT,DHELO,HAXYR,HSEDE,HVIGI,HVIMA,MDISC,NPUMI,PJAPO ortho/genome_ortho/unaln # only retain 1290 single-copy genes

# run PhyloAln to obtain result alignment from read files of each batch
time PhyloAln/PhyloAln -d ortho/genome_ortho/CDS_trim -c 00.config -f fastq -o PhyloAln_00 -p 20 -m codon -b -u DHELO
time PhyloAln/PhyloAln -d ortho/genome_ortho/CDS_trim -c 01.config -f fastq -o PhyloAln_01 -p 20 -m codon -b -u DHELO
time PhyloAln/PhyloAln -d ortho/genome_ortho/CDS_trim -c 02.config -f fastq -o PhyloAln_02 -p 20 -m codon -b -u DHELO
time PhyloAln/PhyloAln -d ortho/genome_ortho/CDS_trim -c 03.config -f fastq -o PhyloAln_03 -p 20 -m codon -b -u DHELO
time PhyloAln/PhyloAln -d ortho/genome_ortho/CDS_trim -c 04.config -f fastq -o PhyloAln_04 -p 20 -m codon -b -u DHELO

# merge the result alignments
time PhyloAln/scripts/merge_seqs.py PhyloAln_all PhyloAln_0*

# calculate the percent completeness and identity of PhyloAln results from reads
for file in ortho/0*.faa
do 
	sp=`basename $file .faa`
	PhyloAln/scripts/test_effect.py ortho/single_copy_ref:$sp PhyloAln_all/nt_out:$sp $sp.PhyloAln.tsv N . .fa ABIPU,CIMPU,CMONT,CSEPT,DHELO,HAXYR,HSEDE,HVIGI,HVIMA,MDISC,NPUMI,PJAPO
done

# phylogenetic reconstruction for PhyloAln based on reads
cd PhyloAln_all
time PhyloAln/scripts/connect.pl -i aa_out -f X -b all.block -n
time iqtree -s all.fas -p all.block -m MFP+MERGE -B 1000 -T AUTO --threads-max 20 --prefix speciestree
cd ..

# run Read2Tree
time read2tree --standalone_path ortho/genome_ortho/unaln --output_path read2tree_out --reference --dna_reference orthogroup/all.cds.fasta
for file in ortho/0*.faa
do
	sp=`basename $file .faa`
	time read2tree --output_path read2tree_out --standalone_path ortho/genome_ortho/unaln --reads dataset/${sp}_1.fq.gz dataset/${sp}_2.fq.gz --threads 20 --dna_reference orthogroup/all.cds.fasta -s $sp
done
time read2tree --standalone_path ortho/genome_ortho/unaln --output_path read2tree_out --merge_all_mappings

# calculate the percent completeness and identity of Read2Tree results
for file in read2tree_out/06_align_merge_dna/*.fa
do
	name=`basename $file`
	readal -in $file -out read2tree_result/$name -fasta
done
PhyloAln/scripts/select_seqs.py ortho/OrthoFinder/Orthogroup_Sequences ABIPU,CIMPU,CMONT,CSEPT,DHELO,HAXYR,HSEDE,HVIGI,HVIMA,MDISC,NPUMI,PJAPO ogs .fa . 1
for file in ortho/single_copy_ref/*.fa
do
	name=`basename $file .fa`
	readal -in read2tree_out/03_align_aa/$name.phy -out read2tree_test/$name.ref_aa.fas -fasta
	mafft --add ogs/$name.fa --localpair --maxiterate 1000 --thread 10 read2tree_test/$name.ref_aa.fas > read2tree_test/$name.aa_aln.fas
	readal -in read2tree_out/03_align_dna/$name.phy -out read2tree_test/$name.dna.fas -fasta
	sed -i 's/-//g' $name.dna.fas
	cat ortho/OrthoFinder/CDS_seq/$name.fa >> read2tree_test/$name.dna.fas
	trimal -in read2tree_test/$name.aa_aln.fas -out read2tree_test/$name.fa -keepheader -backtrans read2tree_test/$name.dna.fas
done
for file in ortho/0*.faa
do
	sp=`basename $file .faa`
	PhyloAln/scripts/test_effect.py read2tree_test:$sp read2tree_result:$sp $sp.read2tree.tsv N . .fa ABIPU,CIMPU,CMONT,CSEPT,DHELO,HAXYR,HSEDE,HVIGI,HVIMA,MDISC,NPUMI,PJAPO
done

# phylogenetic reconstruction for Read2Tree
cd read2tree_out
for file in 06_align_merge_aa/*.fa
do
	name=`basename $file`
	readal -in $file -out align_all/$name -fasta
done
time PhyloAln/scripts/connect.pl -i align_all -f X -b all.block -n
time iqtree -s all.fas -p all.block -m MFP+MERGE -B 1000 -T AUTO --threads-max 20 --prefix speciestree
cd ..

# run PhyloAln to obtain result alignments from the assembly files
time PhyloAln/PhyloAln -d ortho/genome_ortho/CDS_trim -c seq.config -f fasta -o PhyloAln_seq -p 20 -m codon -b -r -u DHELO

# calculate the percent completeness and identity of PhyloAln results from assemblies
for file in ortho/0*.faa
do
	sp=`basename $file .faa`
	PhyloAln/scripts/test_effect.py ortho/single_copy_ref:$sp PhyloAln_seq/nt_out:$sp $sp.PhyloAln_seq.tsv N . .fa ABIPU,CIMPU,CMONT,CSEPT,DHELO,HAXYR,HSEDE,HVIGI,HVIMA,MDISC,NPUMI,PJAPO
done

# phylogenetic reconstruction for PhyloAln based on assemblies
cd PhyloAln_seq
time PhyloAln/scripts/connect.pl -i aa_out -f X -b all.block -n
time iqtree -s all.fas -p all.block -m MFP+MERGE -B 1000 -T AUTO --threads-max 20 --prefix speciestree
cd ..

# run Orthograph
orthograph-manager --create -c orthograph.conf
for sp in `cat genomes.list`
do
	orthograph-manager --load-ogs-peptide ortho/$sp.faa -c orthograph.conf
	orthograph-manager --load-ogs-nucleotide ortho/cds/$sp.cds -c orthograph.conf
done
orthograph-manager orthograph.og.tsv -c orthograph.conf
orthograph-analyzer --prepare -c orthograph.conf
dir=orthograph
for file in 0*.fasta; do
	basename=`basename $file .fasta`
	mkdir $dir/$basename
	cp $file $dir/$basename/$basename.fas
	./rewriteconf.pl $dir $dir/$basename
	time orthograph-analyzer --prepare -c $dir/$basename/$basename.conf
	time orthograph-analyzer -c $dir/$basename/$basename.conf
	time orthograph-reporter -c $dir/$basename/$basename.conf
done
time summarize_orthograph_results.pl -i orthograph -o orthograph_test -d genome.list -c

# calculate the percent completeness and identity of Orthogragh results
for file in ortho/single_copy_ref/*.fa
do
	name=`basename $file .fa`
	mafft --add orthograph_test/aa/$name.aa.summarized.fa --localpair --maxiterate 1000 --thread 20 PhyloAln_00/ref_hmm/$name.aa.fas > orthograph_test/$name.aa_aln.fas
	sed 's/-//g' $file >> orthograph_test/nt/$name.nt.summarized.fa
	trimal -in orthograph_test/$name.aa_aln.fas -out orthograph_test/$name.fa -keepheader -backtrans orthograph_test/nt/$name.nt.summarized.fa
done
for file in ortho/0*.faa
do
	sp=`basename $file .faa`
	PhyloAln/scripts/test_effect.py ortho/single_copy_ref:$sp orthograph_test:$sp $sp.orthograph.tsv N . .fa ABIPU,CIMPU,CMONT,CSEPT,DHELO,HAXYR,HSEDE,HVIGI,HVIMA,MDISC,NPUMI,PJAPO
done

# phylogenetic reconstruction for Orthograph
time summarize_orthograph_results.pl -i orthograph -o orthograph_all -c
cd orthograph_all
for file in aa_summarized/*.aa.summarized.fa
do
	name=`basename $file .aa.summarized.fa`
	time mafft-linsi --thread 20 $file > aa_aln/$name.aln.fas
	time trimal -in aa_aln/$name.aln.fas -out aa_trim/$name.fa -automated1 -keepheader
done
time PhyloAln/scripts/connect.pl -i aa_trim -t orthograph -f X -b all.block -n
time iqtree -s all.fas -p all.block -m MFP+MERGE -B 1000 -T AUTO --threads-max 20 --prefix speciestree
cd ..

# root the trees
PhyloAln/scripts/root_tree.py PhyloAln_all/speciestree.tre PhyloAln_all/speciestree.rooted.tre DHELO
PhyloAln/scripts/root_tree.py read2tree_out/speciestree.tre read2tree_out/speciestree.rooted.tre DHELO
PhyloAln/scripts/root_tree.py PhyloAln_seq/speciestree.tre PhyloAln_seq/speciestree.rooted.tre DHELO
PhyloAln/scripts/root_tree.py orthograph_all/speciestree.tre orthograph_all/speciestree.rooted.tre DHELO
