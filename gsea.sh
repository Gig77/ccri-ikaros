#!/bin/bash

mkdir -p /mnt/projects/ikaros/results/gsea
#cp /mnt/projects/ikaros/results/anduril/execute/rnk_IKNp_vs_IKCp/rnk.csv /mnt/projects/ikaros/results/gsea/IKNp.vs.IKCp.rnk
#cp /mnt/projects/ikaros/results/anduril/execute/rnk_IKMp_vs_IKCp/rnk.csv /mnt/projects/ikaros/results/gsea/IKMp.vs.IKCp.rnk
#cp /mnt/projects/ikaros/results/anduril/execute/rnk_IKD_vs_IKCp/rnk.csv /mnt/projects/ikaros/results/gsea/IKD.vs.IKCp.rnk
#cp /mnt/projects/ikaros/results/anduril_validation/execute/rnk_IKNv_vs_IKCv/rnk.csv /mnt/projects/ikaros/results/gsea/IKNv.vs.IKCv.rnk
#cp /mnt/projects/ikaros/results/anduril_validation/execute/rnk_IKDv_vs_IKCv/rnk.csv /mnt/projects/ikaros/results/gsea/IKDv.vs.IKCv.rnk
cp /mnt/projects/ikaros/results/anduril_validation/execute/rnk_IKMv_vs_IKCv/rnk.csv /mnt/projects/ikaros/results/gsea/IKMv.vs.IKCv.rnk

cd /mnt/projects/ikaros/results/gsea
java -cp /data_synology/software/gsea-2.0.13/gsea2-2.0.13.jar -Xmx5g xtools.gsea.GseaPreranked \
	-rpt_label IKMv.vs.IKCv \
	-rnk /mnt/projects/ikaros/results/gsea/IKMv.vs.IKCv.rnk \
	-gmx /mnt/projects/generic/data/ccri/ccri_custom_gene_sets.gmt,/mnt/projects/generic/data/ccri/ccri_literature_curated_genesets_gsea.gmt \
	-out /mnt/projects/ikaros/results/gsea \
	-plot_top_x 3000 \
	-collapse false \
	-mode Max_probe  \
	-norm meandiv \
	-scoring_scheme weighted \
	-include_only_symbols true \
	-make_sets true \
	-rnd_seed 149 \
	-zip_report false \
	-gui false \
	-nperm 1000 \
	-set_max 5000 \
	-set_min 5