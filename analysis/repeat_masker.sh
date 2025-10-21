awk 'BEGIN{OFS="\t"} {print $6, $7, $8, $12, $2, $10}' repeatmasker.tsv > repeatmasker_class.bed
tail -n +2 "/Volumes/Lucky me/8_STAMP/bin/04.sample_metrics_final_v3/expression_level/repeat/repeatmasker_class.bed" \
| cut -f1-6 \
> "/Volumes/Lucky me/8_STAMP/bin/04.sample_metrics_final_v3/expression_level/repeat/repeatmasker_class_clean.bed"
mv repeatmasker_class_clean.bed repeatmasker_class.bed

bedtools intersect \
-a "/Volumes/Lucky me/8_STAMP/29apr_mincov10/metrics_5kb/region.bed" \
-b "/Volumes/Lucky me/8_STAMP/bin/04.sample_metrics_final_v3/expression_level/repeat/repeatmasker_class.bed" \
-wao \
> "/Volumes/Lucky me/8_STAMP/29apr_mincov10/metrics_5kb/output_repeatmasker_class.tsv"
