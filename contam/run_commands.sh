# simulate the contaminated reads, other preparation steps has been conducted in the commands of the lifetree dataset
# the transcripts of five respectively selected single-copy genes in three target species were removed mannually from ${target_species}.transcript.fasta first
cd dataset
for file in *.transcript.fasta
do
	sp=`basename $file .transcript.fasta`
	art_illumina -ss HS25 -i $sp.transcript.fasta -p -l 100 -f 10 -m 200 -s 10 -o ${sp}trans_10X
done
./sample_reads_contam.py
cd ..

# remove the target sequences
PhyloAln/scripts/select_seqs.py ortho/OrthoFinder/CDS_trim DBUSC,DMOJA,DPSEU,DYAKU,SLEBA ortho/ortho/OGs_fly
PhyloAln/scripts/select_seqs.py ortho/OrthoFinder/Orthogroup_Sequences DBUSC,DMOJA,DPSEU,DYAKU,SLEBA orthogroup/ortho/OGs_fly_unaln # only retain 46 single-copy gene alignments

# run PhyloAln with assembly step
time PhyloAln/PhyloAln -d orthogroup/ortho/OGs_fly -c contam.config -f fasta -o PhyloAln_contam -p 20 -m codon -u SLEBA

# calculate the percent completeness and identity of PhyloAln results with assembly step
PhyloAln/scripts/test_effect.py orthogroup/ortho/OrthoFinder/CDS_trim:DWILL PhyloAln_contam/nt_out:DWILL DWILLcontam.PhyloAln.tsv N . .fa DBUSC,DMOJA,DPSEU,DYAKU,SLEBA
PhyloAln/scripts/test_effect.py orthogroup/ortho/OrthoFinder/CDS_trim:DSIMU PhyloAln_contam/nt_out:DSIMU DSIMUcontam.PhyloAln.tsv N . .fa DBUSC,DMOJA,DPSEU,DYAKU,SLEBA
PhyloAln/scripts/test_effect.py orthogroup/ortho/OrthoFinder/CDS_trim:DMELA PhyloAln_contam/nt_out:DMELA DMELAcontam.PhyloAln.tsv N . .fa DBUSC,DMOJA,DPSEU,DYAKU,SLEBA

# run PhyloAln without assembly step
time PhyloAln/PhyloAln -d orthogroup/ortho/OGs_fly -c contam.config -f fasta -o PhyloAln_contam_b -p 20 -m codon -u SLEBA -b

# calculate the percent completeness and identity of PhyloAln results without assembly step
PhyloAln/scripts/test_effect.py orthogroup/ortho/OrthoFinder/CDS_trim:DWILL PhyloAln_contam_b/nt_out:DWILL DWILLcontam_b.PhyloAln.tsv N . .fa DBUSC,DMOJA,DPSEU,DYAKU,SLEBA
PhyloAln/scripts/test_effect.py orthogroup/ortho/OrthoFinder/CDS_trim:DSIMU PhyloAln_contam_b/nt_out:DSIMU DSIMUcontam_b.PhyloAln.tsv N . .fa DBUSC,DMOJA,DPSEU,DYAKU,SLEBA
PhyloAln/scripts/test_effect.py orthogroup/ortho/OrthoFinder/CDS_trim:DMELA PhyloAln_contam_b/nt_out:DMELA DMELAcontam_b.PhyloAln.tsv N . .fa DBUSC,DMOJA,DPSEU,DYAKU,SLEBA

# run Read2Tree
time read2tree --standalone_path orthogroup/ortho/OGs_fly_unaln --output_path read2tree_contam --reference --dna_reference orthogroup/all.cds.fasta
time read2tree --output_path read2tree_contam --standalone_path orthogroup/ortho/OGs_fly_unaln --reads dataset/DMELAcontam_1.fq dataset/DMELAcontam_2.fq --threads 20 --dna_reference orthogroup/all.cds.fasta -s DMELA
time read2tree --output_path read2tree_contam --standalone_path orthogroup/ortho/OGs_fly_unaln --reads dataset/DSIMUcontam_1.fq dataset/DSIMUcontam_2.fq --threads 20 --dna_reference orthogroup/all.cds.fasta -s DSIMU
time read2tree --output_path read2tree_contam --standalone_path orthogroup/ortho/OGs_fly_unaln --reads dataset/DWILLcontam_1.fq dataset/DWILLcontam_2.fq --threads 20 --dna_reference orthogroup/all.cds.fasta -s DWILL
time read2tree --standalone_path orthogroup/ortho/OGs_fly_unaln --output_path read2tree_contam --merge_all_mappings

# calculate the percent completeness and identity of Read2Tree results
for file in read2tree_contam/06_align_merge_dna/*.fa
do
	name=`basename $file`
	readal -in $file -out read2tree_result/$name -fasta
done
PhyloAln/scripts/select_seqs.py ortho/OrthoFinder/Orthogroup_Sequences DWILL,DMELA,DSIMU ogs
for file in ortho/single_copy_ref/*.fa
do
	name=`basename $file .fa`
	readal -in read2tree_contam/03_align_aa/$name.phy -out read2tree_test/$name.ref_aa.fas -fasta
	mafft --add ogs/$name.fa --keeplength --localpair --maxiterate 1000 --thread 10 read2tree_test/$name.ref_aa.fas > read2tree_test/$name.aa_aln.fas
	readal -in read2tree_contam/03_align_dna/$name.phy -out read2tree_test/$name.dna.fas -fasta
	sed -i 's/-//g' $name.dna.fas
	cat ortho/OrthoFinder/CDS_seq/$name.fa >> read2tree_test/$name.dna.fas
	trimal -in read2tree_test/$name.aa_aln.fas -out read2tree_test/$name.fa -keepheader -backtrans read2tree_test/$name.dna.fas
done
for file in read2tree_contam/05_ogs_map_*_dna
do
	sp=`basename $file _dna | sed 's/^05_ogs_map_//'`
	PhyloAln/scripts/test_effect.py read2tree_test:$sp read2tree_result:$sp $sp.read2tree.tsv N . .fa DBUSC,DMOJA,DPSEU,DYAKU,SLEBA
done

# run PhyloAln to obtain the information of species source of the reads removed ot retained during decontamination, need to manually replace PhyloAln/lib/assemble.py with assemble.debug.py first
PhyloAln/PhyloAln -d orthogroup/ortho/OGs_fly -c contam.config -f fasta -o PhyloAln_contam -p 20 -m codon -u SLEBA
PhyloAln/PhyloAln -d orthogroup/ortho/OGs_fly -c contam.config -f fasta -o PhyloAln_contam_b -p 20 -m codon -u SLEBA -b

# calculate the average percent completeness and identity, and precision and recall of clean, foreign and cross contamination reads of PhyloAln results with or without assembly step
./sum_contam.py
./sum_contam_b.py
# replace the outgroup SLEBA with DBUSC, DMOJA, DPSEU, DYAKU respecitively and run all the commands of PhyloAln and test_effect.py above, in order to test the impact of different outgroups
# calculate the average percent completeness and identity, and precision and recall of clean and foreign contamination reads of PhyloAln results using different outgroups
./sum_contam_outgroup.py
# the calculated precision and recall needs to be further managed in manual, such as removal of invalid zeros, calculation of average values, and summary in a new table

# run PhyloAln to generate the alignments of 'outgroup contamination' from the target species, need to manually replace PhyloAln/lib/assemble.py with assemble.add_out.py first
PhyloAln/PhyloAln -d orthogroup/ortho/OGs_fly -c contam.config -f fasta -o PhyloAln_contam_b_add -p 20 -m codon -u SLEBA -b

# calculate the percent identity of 'outgroup contamination' from the target species
PhyloAln/scripts/test_effect.py orthogroup/ortho/OrthoFinder/CDS_trim:DWILL PhyloAln_contam_b_add/nt_out:DWILL DWILLcontam_b_add.PhyloAln.tsv N . .fa DBUSC,DMOJA,DPSEU,DYAKU,SLEBA
PhyloAln/scripts/test_effect.py orthogroup/ortho/OrthoFinder/CDS_trim:DSIMU PhyloAln_contam_b_add/nt_out:DSIMU DSIMUcontam_b_add.PhyloAln.tsv N . .fa DBUSC,DMOJA,DPSEU,DYAKU,SLEBA
PhyloAln/scripts/test_effect.py orthogroup/ortho/OrthoFinder/CDS_trim:DMELA PhyloAln_contam_b_add/nt_out:DMELA DMELAcontam_b_add.PhyloAln.tsv N . .fa DBUSC,DMOJA,DPSEU,DYAKU,SLEBA

# phylogenetic reconstruction for PhyloAln results
cd PhyloAln_contam_b
PhyloAln/scripts/connect.pl -i aa_out -f X -b all.block -n
time iqtree -s all.fas -p all.block -m MFP+MERGE -B 1000 -T AUTO --threads-max 20 --prefix speciestree
PhyloAln/scripts/connect.pl -i nt_out -f N -b all.block -n -c 123 # need to rename the above results with prefix as 'aatree'
time iqtree -s all.fas -p all.block -m MFP+MERGE -B 1000 -T AUTO --threads-max 20 --prefix speciestree

# phylogenetic reconstruction for 'outgroup contamination' from the target species
cd PhyloAln_contam_b_add
PhyloAln/scripts/connect.pl -i aa_out -f X -b all.block -n
time iqtree -s all.fas -p all.block -m MFP+MERGE -B 1000 -T AUTO --threads-max 20 --prefix speciestree
PhyloAln/scripts/connect.pl -i nt_out -f N -b all.block -n -c 123 # need to rename the above results with prefix as 'aatree'
time iqtree -s all.fas -p all.block -m MFP+MERGE -B 1000 -T AUTO --threads-max 20 --prefix speciestree

# root the trees
PhyloAln/scripts/root_tree.py PhyloAln_contam_b/speciestree.tre PhyloAln_contam_b/speciestree.rooted.tre SLEBA
PhyloAln/scripts/root_tree.py PhyloAln_contam_b/aatree.tre PhyloAln_contam_b/aatree.rooted.tre SLEBA
PhyloAln/scripts/root_tree.py PhyloAln_contam_b_add/aatree.tre PhyloAln_contam_b_add/aatree.rooted.tre SLEBA
PhyloAln/scripts/root_tree.py PhyloAln_contam_b_add/speciestree.tre PhyloAln_contam_b_add/speciestree.rooted.tre SLEBA
