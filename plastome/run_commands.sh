# align the codon sequences of the plastid genes  
for file in test_unaln/*.fa; do  
  name=`basename $file`  
  scripts/alignseq.pl -i $file -o aln/$name -a codon -n 20 -g 11
  mv aln/$name test/$name
done

# remove the target sequences
PhyloAln/scripts/select_seqs.py test Piper_sp._2_Munzinger_6536,Piper_borbonense_1_NCY_873616,Piper_borbonense_2_Ramahefaharivelo_191,Piper_caninum_1_Smith_6545,Piper_caninum_2_Gray_9430,Piper_hederaceum_Gray_9456,Piper_insectifugum_Munzinger_6571,Piper_mestonii_Gray_9458,Piper_sarmentosum_Smith_s.n.1,Piper_sp._3_Davidson_12685,Piper_commutatum_Angel_2356,Piper_darienense_Aguilar_11197,Piper_umbellatum_Davidson_12352,Piper_capense_1_Davidson_11009,Piper_capense_2_Smith_4925,Piper_guahamense_Flynn_6748 ref .fa . 1

# run PhyloAln
time /public/lihaosen/PhyloAln/PhyloAln -d ref -c plastid.config -f fastq -p 20 -m codon -g 11 -b -u Magnolia_grandiflora -r

# calculate the percent completeness and identity of PhyloAln results
PhyloAln/scripts/test_effect.py test:Piper_sp._2_Munzinger_6536 PhyloAln_out/nt_out:Piper_sp2 Piper_sp2.plastome.tsv N . .fa Drimys_granadensis,Aristolochia_contorta,Aristolochia_debilis,Magnolia_grandiflora,Piper_auritum,Piper_cenocladum,Piper_kadsura,Piper_laetispicum,Piper_nigrum,Saruma_henryi,Peperomia_berteroana_JBN_3340,Peperomia_congesta_HW_6593,Peperomia_fernandopoiana_Razoazanany_97,Peperomia_glabella_Stevens_34572,Peperomia_margaritifera_JBN_3342,Peperomia_tetraphylla_Lim_74,Piper_ponapense_Wood_13509,Piper_puberulum_Lorence_10610,Piper_sp._1_McPherson_19040
PhyloAln/scripts/test_effect.py test:Piper_borbonense_1_NCY_873616 PhyloAln_out/nt_out:Piper_borbonense1 Piper_borbonense1.plastome.tsv N . .fa Drimys_granadensis,Aristolochia_contorta,Aristolochia_debilis,Magnolia_grandiflora,Piper_auritum,Piper_cenocladum,Piper_kadsura,Piper_laetispicum,Piper_nigrum,Saruma_henryi,Peperomia_berteroana_JBN_3340,Peperomia_congesta_HW_6593,Peperomia_fernandopoiana_Razoazanany_97,Peperomia_glabella_Stevens_34572,Peperomia_margaritifera_JBN_3342,Peperomia_tetraphylla_Lim_74,Piper_ponapense_Wood_13509,Piper_puberulum_Lorence_10610,Piper_sp._1_McPherson_19040
PhyloAln/scripts/test_effect.py test:Piper_borbonense_2_Ramahefaharivelo_191 PhyloAln_out/nt_out:Piper_borbonense2 Piper_borbonense2.plastome.tsv N . .fa Drimys_granadensis,Aristolochia_contorta,Aristolochia_debilis,Magnolia_grandiflora,Piper_auritum,Piper_cenocladum,Piper_kadsura,Piper_laetispicum,Piper_nigrum,Saruma_henryi,Peperomia_berteroana_JBN_3340,Peperomia_congesta_HW_6593,Peperomia_fernandopoiana_Razoazanany_97,Peperomia_glabella_Stevens_34572,Peperomia_margaritifera_JBN_3342,Peperomia_tetraphylla_Lim_74,Piper_ponapense_Wood_13509,Piper_puberulum_Lorence_10610,Piper_sp._1_McPherson_19040
PhyloAln/scripts/test_effect.py test:Piper_caninum_1_Smith_6545 PhyloAln_out/nt_out:Piper_caninum1 Piper_caninum1.plastome.tsv N . .fa Drimys_granadensis,Aristolochia_contorta,Aristolochia_debilis,Magnolia_grandiflora,Piper_auritum,Piper_cenocladum,Piper_kadsura,Piper_laetispicum,Piper_nigrum,Saruma_henryi,Peperomia_berteroana_JBN_3340,Peperomia_congesta_HW_6593,Peperomia_fernandopoiana_Razoazanany_97,Peperomia_glabella_Stevens_34572,Peperomia_margaritifera_JBN_3342,Peperomia_tetraphylla_Lim_74,Piper_ponapense_Wood_13509,Piper_puberulum_Lorence_10610,Piper_sp._1_McPherson_19040
PhyloAln/scripts/test_effect.py test:Piper_caninum_2_Gray_9430 PhyloAln_out/nt_out:Piper_caninum2 Piper_caninum2.plastome.tsv N . .fa Drimys_granadensis,Aristolochia_contorta,Aristolochia_debilis,Magnolia_grandiflora,Piper_auritum,Piper_cenocladum,Piper_kadsura,Piper_laetispicum,Piper_nigrum,Saruma_henryi,Peperomia_berteroana_JBN_3340,Peperomia_congesta_HW_6593,Peperomia_fernandopoiana_Razoazanany_97,Peperomia_glabella_Stevens_34572,Peperomia_margaritifera_JBN_3342,Peperomia_tetraphylla_Lim_74,Piper_ponapense_Wood_13509,Piper_puberulum_Lorence_10610,Piper_sp._1_McPherson_19040
PhyloAln/scripts/test_effect.py test:Piper_capense_1_Davidson_11009 PhyloAln_out/nt_out:Piper_capense1 Piper_capense1.plastome.tsv N . .fa Drimys_granadensis,Aristolochia_contorta,Aristolochia_debilis,Magnolia_grandiflora,Piper_auritum,Piper_cenocladum,Piper_kadsura,Piper_laetispicum,Piper_nigrum,Saruma_henryi,Peperomia_berteroana_JBN_3340,Peperomia_congesta_HW_6593,Peperomia_fernandopoiana_Razoazanany_97,Peperomia_glabella_Stevens_34572,Peperomia_margaritifera_JBN_3342,Peperomia_tetraphylla_Lim_74,Piper_ponapense_Wood_13509,Piper_puberulum_Lorence_10610,Piper_sp._1_McPherson_19040
PhyloAln/scripts/test_effect.py test:Piper_capense_2_Smith_4925 PhyloAln_out/nt_out:Piper_capense2 Piper_capense2.plastome.tsv N . .fa Drimys_granadensis,Aristolochia_contorta,Aristolochia_debilis,Magnolia_grandiflora,Piper_auritum,Piper_cenocladum,Piper_kadsura,Piper_laetispicum,Piper_nigrum,Saruma_henryi,Peperomia_berteroana_JBN_3340,Peperomia_congesta_HW_6593,Peperomia_fernandopoiana_Razoazanany_97,Peperomia_glabella_Stevens_34572,Peperomia_margaritifera_JBN_3342,Peperomia_tetraphylla_Lim_74,Piper_ponapense_Wood_13509,Piper_puberulum_Lorence_10610,Piper_sp._1_McPherson_19040
PhyloAln/scripts/test_effect.py test:Piper_commutatum_Angel_2356 PhyloAln_out/nt_out:Piper_commutatum Piper_commutatum.plastome.tsv N . .fa Drimys_granadensis,Aristolochia_contorta,Aristolochia_debilis,Magnolia_grandiflora,Piper_auritum,Piper_cenocladum,Piper_kadsura,Piper_laetispicum,Piper_nigrum,Saruma_henryi,Peperomia_berteroana_JBN_3340,Peperomia_congesta_HW_6593,Peperomia_fernandopoiana_Razoazanany_97,Peperomia_glabella_Stevens_34572,Peperomia_margaritifera_JBN_3342,Peperomia_tetraphylla_Lim_74,Piper_ponapense_Wood_13509,Piper_puberulum_Lorence_10610,Piper_sp._1_McPherson_19040
PhyloAln/scripts/test_effect.py test:Piper_darienense_Aguilar_11197 PhyloAln_out/nt_out:Piper_darienense Piper_darienense.plastome.tsv N . .fa Drimys_granadensis,Aristolochia_contorta,Aristolochia_debilis,Magnolia_grandiflora,Piper_auritum,Piper_cenocladum,Piper_kadsura,Piper_laetispicum,Piper_nigrum,Saruma_henryi,Peperomia_berteroana_JBN_3340,Peperomia_congesta_HW_6593,Peperomia_fernandopoiana_Razoazanany_97,Peperomia_glabella_Stevens_34572,Peperomia_margaritifera_JBN_3342,Peperomia_tetraphylla_Lim_74,Piper_ponapense_Wood_13509,Piper_puberulum_Lorence_10610,Piper_sp._1_McPherson_19040
PhyloAln/scripts/test_effect.py test:Piper_guahamense_Flynn_6748 PhyloAln_out/nt_out:Piper_guahamense Piper_guahamense.plastome.tsv N . .fa Drimys_granadensis,Aristolochia_contorta,Aristolochia_debilis,Magnolia_grandiflora,Piper_auritum,Piper_cenocladum,Piper_kadsura,Piper_laetispicum,Piper_nigrum,Saruma_henryi,Peperomia_berteroana_JBN_3340,Peperomia_congesta_HW_6593,Peperomia_fernandopoiana_Razoazanany_97,Peperomia_glabella_Stevens_34572,Peperomia_margaritifera_JBN_3342,Peperomia_tetraphylla_Lim_74,Piper_ponapense_Wood_13509,Piper_puberulum_Lorence_10610,Piper_sp._1_McPherson_19040
PhyloAln/scripts/test_effect.py test:Piper_hederaceum_Gray_9456 PhyloAln_out/nt_out:Piper_hederaceum Piper_hederaceum.plastome.tsv N . .fa Drimys_granadensis,Aristolochia_contorta,Aristolochia_debilis,Magnolia_grandiflora,Piper_auritum,Piper_cenocladum,Piper_kadsura,Piper_laetispicum,Piper_nigrum,Saruma_henryi,Peperomia_berteroana_JBN_3340,Peperomia_congesta_HW_6593,Peperomia_fernandopoiana_Razoazanany_97,Peperomia_glabella_Stevens_34572,Peperomia_margaritifera_JBN_3342,Peperomia_tetraphylla_Lim_74,Piper_ponapense_Wood_13509,Piper_puberulum_Lorence_10610,Piper_sp._1_McPherson_19040
PhyloAln/scripts/test_effect.py test:Piper_insectifugum_Munzinger_6571 PhyloAln_out/nt_out:Piper_insectifugum Piper_insectifugum.plastome.tsv N . .fa Drimys_granadensis,Aristolochia_contorta,Aristolochia_debilis,Magnolia_grandiflora,Piper_auritum,Piper_cenocladum,Piper_kadsura,Piper_laetispicum,Piper_nigrum,Saruma_henryi,Peperomia_berteroana_JBN_3340,Peperomia_congesta_HW_6593,Peperomia_fernandopoiana_Razoazanany_97,Peperomia_glabella_Stevens_34572,Peperomia_margaritifera_JBN_3342,Peperomia_tetraphylla_Lim_74,Piper_ponapense_Wood_13509,Piper_puberulum_Lorence_10610,Piper_sp._1_McPherson_19040
PhyloAln/scripts/test_effect.py test:Piper_mestonii_Gray_9458 PhyloAln_out/nt_out:Piper_mestonii Piper_mestonii.plastome.tsv N . .fa Drimys_granadensis,Aristolochia_contorta,Aristolochia_debilis,Magnolia_grandiflora,Piper_auritum,Piper_cenocladum,Piper_kadsura,Piper_laetispicum,Piper_nigrum,Saruma_henryi,Peperomia_berteroana_JBN_3340,Peperomia_congesta_HW_6593,Peperomia_fernandopoiana_Razoazanany_97,Peperomia_glabella_Stevens_34572,Peperomia_margaritifera_JBN_3342,Peperomia_tetraphylla_Lim_74,Piper_ponapense_Wood_13509,Piper_puberulum_Lorence_10610,Piper_sp._1_McPherson_19040
PhyloAln/scripts/test_effect.py test:Piper_sarmentosum_Smith_s.n.1 PhyloAln_out/nt_out:Piper_sarmentosum Piper_sarmentosum.plastome.tsv N . .fa Drimys_granadensis,Aristolochia_contorta,Aristolochia_debilis,Magnolia_grandiflora,Piper_auritum,Piper_cenocladum,Piper_kadsura,Piper_laetispicum,Piper_nigrum,Saruma_henryi,Peperomia_berteroana_JBN_3340,Peperomia_congesta_HW_6593,Peperomia_fernandopoiana_Razoazanany_97,Peperomia_glabella_Stevens_34572,Peperomia_margaritifera_JBN_3342,Peperomia_tetraphylla_Lim_74,Piper_ponapense_Wood_13509,Piper_puberulum_Lorence_10610,Piper_sp._1_McPherson_19040
PhyloAln/scripts/test_effect.py test:Piper_sp._3_Davidson_12685 PhyloAln_out/nt_out:Piper_sp3 Piper_sp3.plastome.tsv N . .fa Drimys_granadensis,Aristolochia_contorta,Aristolochia_debilis,Magnolia_grandiflora,Piper_auritum,Piper_cenocladum,Piper_kadsura,Piper_laetispicum,Piper_nigrum,Saruma_henryi,Peperomia_berteroana_JBN_3340,Peperomia_congesta_HW_6593,Peperomia_fernandopoiana_Razoazanany_97,Peperomia_glabella_Stevens_34572,Peperomia_margaritifera_JBN_3342,Peperomia_tetraphylla_Lim_74,Piper_ponapense_Wood_13509,Piper_puberulum_Lorence_10610,Piper_sp._1_McPherson_19040
PhyloAln/scripts/test_effect.py test:Piper_umbellatum_Davidson_12352 PhyloAln_out/nt_out:Piper_umbellatum Piper_umbellatum.plastome.tsv N . .fa Drimys_granadensis,Aristolochia_contorta,Aristolochia_debilis,Magnolia_grandiflora,Piper_auritum,Piper_cenocladum,Piper_kadsura,Piper_laetispicum,Piper_nigrum,Saruma_henryi,Peperomia_berteroana_JBN_3340,Peperomia_congesta_HW_6593,Peperomia_fernandopoiana_Razoazanany_97,Peperomia_glabella_Stevens_34572,Peperomia_margaritifera_JBN_3342,Peperomia_tetraphylla_Lim_74,Piper_ponapense_Wood_13509,Piper_puberulum_Lorence_10610,Piper_sp._1_McPherson_19040
  
# phylogenetic reconstruction for PhyloAln
cd PhyloAln_out
time PhyloAln/scripts/connect.pl -i nt_out -f N -b all.block -n -c 123
time iqtree -s all.fas -p all.block -m MFP+MERGE -B 1000 -T AUTO --threads-max 20 --prefix speciestree
cd ..

# phylogenetic reconstruction for de novo alignments using predicted genes
time PhyloAln/scripts/connect.pl -i ref -f N -b all.block -n -c 123
time iqtree -s all.fas -p all.block -m MFP+MERGE -B 1000 -T AUTO --threads-max 20 --prefix speciestree

# root the trees
PhyloAln/scripts/root_tree.py PhyloAln_out/speciestree.tre PhyloAln_out/speciestree.rooted.tre Magnolia_grandiflora
PhyloAln/scripts/root_tree.py speciestree.tre speciestree.rooted.tre Magnolia_grandiflora