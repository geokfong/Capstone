INPUT_DIR="/Volumes/Lucky me/8_STAMP/bin/03.comparative_analysis"
REGION_BED="/Volumes/Lucky me/8_STAMP/genemodels/temporary/To_map.bed"
OUTPUT_DIR="/Volumes/Lucky me/8_STAMP/29apr_mincov10/metrics"
mkdir -p "${OUTPUT_DIR}"

# Define group to type/gene mapping
declare -A GROUP_TYPE_MAP=(
    ["APO_FTO"]="apo_fto"
    ["APO_RIPK1"]="apo_ripk1"
    ["TAD_FTO"]="tad_fto"
)

declare -A GROUP_GENE_MAP=(
    ["APO_FTO"]="FTO"
    ["APO_RIPK1"]="RIPK1"
    ["TAD_FTO"]="FTO"
)

declare -A SAMPLE_GROUPS=(
    ["APO_FTO"]="SKOV3_Apo-FTO_BT1_ID10_R45 SKOV3_Apo-FTO_BT1_ID11_R46 SKOV3_Apo-FTO_UT_ID7_R43 SKOV3_Apo-FTO_UT_ID9_R44"
    ["APO_RIPK1"]="SKOV3_Apo-RIPK1_BT1_ID11_R11 SKOV3_Apo-RIPK1_BT1_ID12_R11_2 SKOV3_Apo-RIPK1_BT1_ID12_R12 SKOV3_Apo-RIPK1_BT1_ID13_R12_2 SKOV3_Apo-RIPK1_UT_ID5_R5 SKOV3_Apo-RIPK1_UT_ID6_R5_2 SKOV3_Apo-RIPK1_UT_ID6_R6 SKOV3_Apo-RIPK1_UT_ID7_R6_2"
    ["TAD_FTO"]="SKOV3_Tad-FTO_BT1_ID10_R64 SKOV3_Tad-FTO_BT1_ID9_R63 SKOV3_Tad-FTO_UT_ID7_R61 SKOV3_Tad-FTO_UT_ID8_R62"
)

for group in "${!SAMPLE_GROUPS[@]}"; do
    echo "Processing group: ${group}"
    mkdir -p "${OUTPUT_DIR}/${group}"
    
    # Get type and gene for this group
    type="${GROUP_TYPE_MAP[$group]}"
    gene="${GROUP_GENE_MAP[$group]}"
    
    for sample in ${SAMPLE_GROUPS[${group}]}; do
        echo "  Processing sample: ${sample}"
        
        input_bed="$INPUT_DIR/${type}/${gene}/${sample}_clean.bed"
        output_metrics_file="${OUTPUT_DIR}/${group}/${sample}_mincov10_metrics.bed"

        if [ ! -f "$input_bed" ]; then
            echo "Warning: Input file not found - $input_bed"
            continue
        fi

        if [ ! -f "${REGION_BED}" ]; then
            echo "Error: Region BED file not found - ${REGION_BED}"
            exit 1
        fi

        echo -e "chrom\tstart\tend\tfeature_info\tscore\tstrand\tedit_count\tedit_mean\tedit_sum" > "$output_metrics_file"

        if ! bedtools map \
            -a "${REGION_BED}" \
            -b "$input_bed" \
            -c 5 \
            -o count,mean,sum \
            >> "$output_metrics_file"; then
            echo "Error processing $sample"
            continue
        fi

        echo "  Completed processing $sample"
    done
done

echo "Processing complete. Results saved to: $OUTPUT_DIR"

INPUT_DIR="/Volumes/Lucky me/8_STAMP/29apr_mincov10/metrics"
OUTPUT_FILE="/Volumes/Lucky me/8_STAMP/29apr_mincov10/metrics/all_metrics_mincov10_v2.tsv"
TMP_DIR="tmp_processing"

# Clean up temp directory if it exists and create new
rm -rf "${TMP_DIR}"
mkdir -p "${TMP_DIR}"

# Define sample groups
declare -A SAMPLE_GROUPS=(
    ["APO_FTO"]="SKOV3_Apo-FTO_BT1_ID10_R45 SKOV3_Apo-FTO_BT1_ID11_R46 SKOV3_Apo-FTO_UT_ID7_R43 SKOV3_Apo-FTO_UT_ID9_R44"
    ["APO_RIPK1"]="SKOV3_Apo-RIPK1_BT1_ID11_R11 SKOV3_Apo-RIPK1_BT1_ID12_R11_2 SKOV3_Apo-RIPK1_BT1_ID12_R12 SKOV3_Apo-RIPK1_BT1_ID13_R12_2 SKOV3_Apo-RIPK1_UT_ID5_R5 SKOV3_Apo-RIPK1_UT_ID6_R5_2 SKOV3_Apo-RIPK1_UT_ID6_R6 SKOV3_Apo-RIPK1_UT_ID7_R6_2"
    ["TAD_FTO"]="SKOV3_Tad-FTO_BT1_ID10_R64 SKOV3_Tad-FTO_BT1_ID9_R63 SKOV3_Tad-FTO_UT_ID7_R61 SKOV3_Tad-FTO_UT_ID8_R62"
)

first_group=$(echo "${!SAMPLE_GROUPS[@]}" | awk '{print $1}')
first_sample="${SAMPLE_GROUPS[$first_group]%% *}"
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6}' \
    "${INPUT_DIR}/${first_group}/${first_sample}_mincov10_metrics.bed" \
    | sort -k1,1 -k2,2n > "${TMP_DIR}/base_regions.bed"

cp "${TMP_DIR}/base_regions.bed" "${TMP_DIR}/combined.bed"

HEADER="chr\tstart\tend\tfeature\tscore\tstrand"
for group in "${!SAMPLE_GROUPS[@]}"; do
    for sample in ${SAMPLE_GROUPS[$group]}; do
        HEADER+="\t${sample}_count\t${sample}_mean\t${sample}_sum"
    done
done
echo -e "$HEADER" > "$OUTPUT_FILE"

for group in "${!SAMPLE_GROUPS[@]}"; do
    for sample in ${SAMPLE_GROUPS[$group]}; do
        input_file="${INPUT_DIR}/${group}/${sample}_mincov10_metrics.bed"
        
        echo "Processing $input_file"
        
        # Sort and extract metrics columns
        sort -k1,1 -k2,2n "$input_file" | \
        awk 'BEGIN{OFS="\t"} {print $7,$8,$9}' > "${TMP_DIR}/${sample}_metrics.tmp"
        
        # Paste with combined file
        paste "${TMP_DIR}/combined.bed" "${TMP_DIR}/${sample}_metrics.tmp" \
            > "${TMP_DIR}/new_combined.bed"
        
        mv "${TMP_DIR}/new_combined.bed" "${TMP_DIR}/combined.bed"
    done
done

cat "${TMP_DIR}/combined.bed" >> "$OUTPUT_FILE"

echo "Combined file created: $OUTPUT_FILE"
