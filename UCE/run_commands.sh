# remove the target sequences
PhyloAln/scripts/select_seqs.py . Anolis_carolinensis,Chelonia_mydas,Chrysemys_picta,Dermatemys_sp,Homo_sapiens19,Pelodiscus_sinensis,Python_molurus,Sphenodon_punctatus,Trachemys_scripta,Crocodylus_porosus,Gallus_gallus ref .fasta
mv ref/all.fasta ref.fasta

# run PhyloAln
time /public/lihaosen/PhyloAln/PhyloAln -a ref.fasta -c UCE.config -f fastq -p 20 -u Homo_sapiens19 --ref_split_len 1000 -b

# calculate the percent completeness and identity of PhyloAln results
cp PhyloAln_out/nt_out/aln.concatenated.fa PhyloAln_out/nt_out/all.fasta
for sp in `cut -f 1 UCE.config`
do
	../scripts/test_effect.py .:$sp PhyloAln_out/nt_out:$sp $sp.UCE.tsv N . .fasta Anolis_carolinensis,Chelonia_mydas,Chrysemys_picta,Dermatemys_sp,Homo_sapiens19,Pelodiscus_sinensis,Python_molurus,Sphenodon_punctatus,Trachemys_scripta,Crocodylus_porosus,Gallus_gallus
done

# phylogenetic reconstruction for PhyloAln
time iqtree -s PhyloAln_out/nt_out/aln.concatenated.fa -B 1000 -T AUTO --threads-max 20 --prefix speciestree

# phylogenetic reconstruction for downloaded matrix
time iqtree -s all.fasta -B 1000 -T AUTO --threads-max 20 --prefix reftree

# root the trees
PhyloAln/scripts/root_tree.py reftree.tre reftree.rooted.tre Homo_sapiens19
PhyloAln/scripts/root_tree.py speciestree.tre speciestree.rooted.tre Homo_sapiens19

# assemble the reads of Graptemys pseudogeographica
fq2fa --merge --filter Graptemys_pseudogeographica_1.fq Graptemys_pseudogeographica_2.fq Graptemys_pseudogeographica.fa
idba_ud -r Graptemys_pseudogeographica.fa -o Graptemys_pseudogeographica_test --num_threads 35

# check the source of Graptemys pseudogeographica
sed 's/-//g' all.fasta > all.nogap.fa
sed -i 's/N//g' all.nogap.fa
makeblastdb -in all.nogap.fa -dbtype nucl
blastn -query Graptemys_pseudogeographica_test/scaffold.fa -db all.nogap.fa -evalue 1e-5 -max_target_seqs 1 -max_hsps 1 -outfmt 6 -num_threads 35 -out Graptemys_pseudogeographica.test.tsv
