#!/bin/bash
set -o xtrace

sample_info=input/sample_info.txt
g1name=grp1
g2name=grp2
outdir=output
pkldir=output/intermediatefiles
gene_annotation=input/annotations/ENSG_to_genename.txt
transcript_annotation=input/annotations/ENST.txt

date

disco \
--outdir $outdir \
--pkldir $pkldir \
--group1color "#2BF0F0" \
--group2color "#FF8686" \
--geneannotationfile $gene_annotation \
--transcriptannotationfile $transcript_annotation \
--maxciwidth 1 \
--mininfreads 3 \
--mindefreads 0 \
--minavgpsi 0.1 \
--minnumcells 50 \
--minavgshift 0.1 \
--stattest KS \
--alpha 0.05 \
--multitestmethod "fdr_bh" \
$sample_info $g1name $g2name

date
