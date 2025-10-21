wget "https://hgdownload.soe.ucsc.edu/gbdb/hg38/bbi/clinvar/clinvarMain.bb"

curl -O https://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/bigBedToBed

chmod +x bigBedToBed
sudo mv bigBedToBed /usr/local/bin/
which bigBedToBed

bigBedToBed clinvarMain.bb clinvarMain.bed

bedtools intersect -a "/Volumes/Lucky me/8_STAMP/29apr_mincov10/metrics_5kb/region.bed" \
-b "/Volumes/Lucky me/8_STAMP/bin/04.sample_metrics_final_v3/expression_level/disease_expression/clinvarMain.bed" \
-wao \
> "/Volumes/Lucky me/8_STAMP/29apr_mincov10/metrics_5kb/output_clinvarMain.tsv"

cut -f 2- "/Volumes/Lucky me/8_STAMP/bin/04.sample_metrics_final_v3/expression_level/disease_expression/GWAS" \
> "/Volumes/Lucky me/8_STAMP/bin/04.sample_metrics_final_v3/expression_level/disease_expression/GWAS_update.tsv"

bedtools intersect -a "/Volumes/Lucky me/8_STAMP/29apr_mincov10/metrics_5kb/region.bed" \
-b "/Volumes/Lucky me/8_STAMP/bin/04.sample_metrics_final_v3/expression_level/disease_expression/GWAS_update.tsv" \
-wao \
> "/Volumes/Lucky me/8_STAMP/29apr_mincov10/metrics_5kb/output_GWAS.tsv"
