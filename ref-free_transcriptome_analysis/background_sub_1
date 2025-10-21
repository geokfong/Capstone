# Background subtraction A:	Enzyme-Specific Filtering - Opposite Enzyme Background Removal
APO_FTO_CTRL=(
    "SKOV3_Apo-Cnt-FTO_BT1_ID3_R39"
    "SKOV3_Apo-Cnt-FTO_BT1_ID4_R40"
    "SKOV3_Apo-Cnt-FTO_UT_ID1_R37"
    "SKOV3_Apo-Cnt-FTO_UT_ID2_R38"
)

APO_RIPK1_CTRL=(
    "SKOV3_Apo-Cnt-RIPK1_BT1_ID10_R10"
    "SKOV3_Apo-Cnt-RIPK1_BT1_ID10_R9_2"
    "SKOV3_Apo-Cnt-RIPK1_BT1_ID11_R10_2"
    "SKOV3_Apo-Cnt-RIPK1_BT1_ID9_R9"
    "SKOV3_Apo-Cnt-RIPK1_UT_ID1_R1"
    "SKOV3_Apo-Cnt-RIPK1_UT_ID2_R1_2"
    "SKOV3_Apo-Cnt-RIPK1_UT_ID2_R2"
    "SKOV3_Apo-Cnt-RIPK1_UT_ID3_R2_2"
)

APO_FTO_SAMPLES=(
    "SKOV3_Apo-FTO_BT1_ID10_R45"
    "SKOV3_Apo-FTO_BT1_ID11_R46"
    "SKOV3_Apo-FTO_UT_ID7_R43"
    "SKOV3_Apo-FTO_UT_ID9_R44"
)

APO_RIPK1_SAMPLES=(
    "SKOV3_Apo-RIPK1_BT1_ID11_R11"
    "SKOV3_Apo-RIPK1_BT1_ID12_R11_2"
    "SKOV3_Apo-RIPK1_BT1_ID12_R12"
    "SKOV3_Apo-RIPK1_BT1_ID13_R12_2"
    "SKOV3_Apo-RIPK1_UT_ID5_R5"
    "SKOV3_Apo-RIPK1_UT_ID6_R5_2"
    "SKOV3_Apo-RIPK1_UT_ID6_R6"
    "SKOV3_Apo-RIPK1_UT_ID7_R6_2"
)

TAD_FTO_CTRL=(
    "SKOV3_Tad-Cnt-FTO_BT1_ID3_R57"
    "SKOV3_Tad-Cnt-FTO_BT1_ID4_R58"
    "SKOV3_Tad-Cnt-FTO_UT_ID1_R55"
    "SKOV3_Tad-Cnt-FTO_UT_ID2_R56"
)

TAD_FTO_SAMPLES=(
    "SKOV3_Tad-FTO_BT1_ID10_R64"
    "SKOV3_Tad-FTO_BT1_ID9_R63"
    "SKOV3_Tad-FTO_UT_ID7_R61"
    "SKOV3_Tad-FTO_UT_ID8_R62"
)

#parameters
MIN_COVERAGE=1
#MIN_EDIT_RATE=0.01

mkdir -p /Volumes/Lucky_me/9_STAMP_transcript/05.transcript_analysis/{apo_fto,apo_ripk1,tad_fto}/{control,FTO,RIPK1,final}
output_dir="/Volumes/Lucky_me/9_STAMP_transcript/05.transcript_analysis"
input_dir_tad="/Volumes/Lucky_me/9_STAMP_transcript/pileup2var"
input_dir_apo="/Volumes/Lucky_me/9_STAMP_transcript/pileup2var"

for sample in "${APO_FTO_CTRL[@]}" "${APO_RIPK1_CTRL[@]}" "${APO_FTO_SAMPLES[@]}" "${APO_RIPK1_SAMPLES[@]}"; do
    if [ ! -f "${input_dir_apo}/${sample}.Transcript_CT.pileup2var.flt.txt" ]; then
        echo "Error: Input file ${input_dir_apo}/${sample}.Transcript_CT.pileup2var.flt.txt not found"
    fi
done

for sample in "${TAD_FTO_CTRL[@]}" "${TAD_FTO_SAMPLES[@]}"; do
    if [ ! -f "${input_dir_tad}/${sample}.Transcript_AG.pileup2var.flt.txt" ]; then
        echo "Error: Input file ${input_dir_tad}/${sample}.Transcript_AG.pileup2var.flt.txt not found"
    fi
done

# Function to process background variants
process_background() {
    local sample=$1
    local output_loc=$2
    local conv_type=$3
    
    if [ "$conv_type" = "apo" ]; then
        awk -v min_cov=$MIN_COVERAGE '
        NR>1 && $4=="A" && $5>=min_cov && $9>0 {
            print $1"\t"$2-1"\t"$2"\t"$4"\t"$9/$5"\t"$3
        }' $input_dir_apo/${sample}.Transcript_CT.pileup2var.flt.txt > $output_loc/${sample}_background.bed
    else
        awk -v min_cov=$MIN_COVERAGE '
        NR>1 && $4=="C" && $5>=min_cov && $7>0 {
            print $1"\t"$2-1"\t"$2"\t"$4"\t"$7/$5"\t"$3
        }' $input_dir_tad/${sample}.Transcript_AG.pileup2var.flt.txt > $output_loc/${sample}_background.bed
    fi
}

echo "Processing Apo-fto ctrl..."
for sample in "${APO_FTO_CTRL[@]}"; do
    process_background $sample "$output_dir/apo_fto/control" "apo"
done

echo "Processing Apo-ripk1 ctrl..."
for sample in "${APO_RIPK1_CTRL[@]}"; do
    process_background $sample "$output_dir/apo_ripk1/control" "apo"
done

echo "Processing Tad-fto ctrl..."
for sample in "${TAD_FTO_CTRL[@]}"; do
    process_background $sample "$output_dir/tad_fto/control" "tad"
done


MIN_COVERAGE=10
process_sample() {
    local sample=$1
    local output_loc=$2
    local conv_type=$3
    
    if [ "$conv_type" = "apo" ]; then
        # C>T conversion for Apo
        awk -v min_cov=$MIN_COVERAGE '
        NR>1 && $4=="C" && $5>=min_cov && $7>0 {
            edit_rate = $7/$5
            print $1"\t"$2-1"\t"$2"\t"$4"\t"edit_rate"\t"$5"\t"$7"\t"$3
        }' $input_dir_apo/${sample}.Transcript_CT.pileup2var.flt.txt > $output_loc/${sample}_edits.bed
    else
        # A>G conversion for Tad
        awk -v min_cov=$MIN_COVERAGE '
        NR>1 && $4=="A" && $5>=min_cov && $9>0 {
            edit_rate = $9/$5
            print $1"\t"$2-1"\t"$2"\t"$4"\t"edit_rate"\t"$5"\t"$9"\t"$3
        }' $input_dir_tad/${sample}.Transcript_AG.pileup2var.flt.txt > $output_loc/${sample}_edits.bed
    fi
}


for sample in "${APO_FTO_SAMPLES[@]}"; do
    process_sample $sample "$output_dir/apo_fto/FTO" "apo"
done

for sample in "${APO_RIPK1_SAMPLES[@]}"; do
    process_sample $sample "$output_dir/apo_ripk1/RIPK1" "apo"
done

for sample in "${TAD_FTO_SAMPLES[@]}"; do
    process_sample $sample "$output_dir/tad_fto/FTO" "tad"
done

#cat apo and tad separately
cat "$output_dir"/apo_fto/control/*_background.bed \
    "$output_dir"/apo_ripk1/control/*_background.bed | \
    sort -k1,1 -k2,2n | \
    bedtools merge -i stdin > "$output_dir"/SKOV3_Apo-Combined_background_variants.bed

cat "$output_dir"/tad_fto/control/*_background.bed | \
    sort -k1,1 -k2,2n | \
    bedtools merge -i stdin > "$output_dir"/SKOV3_Tad-Cnt-FTO_background_variants.bed

process_group() {
    local sample=$1
    local type=$2
    local gene=$3
    
    # Remove background variants for the sample
    echo "Processing ${type} sample: $sample"
    bedtools subtract \
        -a "$output_dir/${type}/${gene}/${sample}_edits.bed" \
        -b "$output_dir/SKOV3_Apo-Combined_background_variants.bed" | \
    sort -k1,1 -k2,2n | \
    awk '$5 != 1' > "$output_dir/${type}/${gene}/${sample}_clean_v2.bed"
}

for sample in "${TAD_FTO_SAMPLES[@]}"; do
    process_group "$sample" "tad_fto" "FTO"
done

process_group() {
    local sample=$1
    local type=$2
    local gene=$3
    
    # Remove background variants for the sample
    echo "Processing ${type} sample: $sample"
    bedtools subtract \
        -a "$output_dir/${type}/${gene}/${sample}_edits.bed" \
        -b "$output_dir/SKOV3_Tad-Cnt-FTO_background_variants.bed" | \
    sort -k1,1 -k2,2n | \
    awk '$5 != 1' > "$output_dir/${type}/${gene}/${sample}_clean_v2.bed"
}

for sample in "${APO_FTO_SAMPLES[@]}"; do
    process_group "$sample" "apo_fto" "FTO"
done

for sample in "${APO_RIPK1_SAMPLES[@]}"; do
    process_group "$sample" "apo_ripk1" "RIPK1"
done

echo -e "chromosome\tstart\tend\tapo_fto_BT_vs_CTRL_edit_rate\tapo_fto_UT_vs_CTRL_edit_rate\tapo_ripk1_BT_vs_CTRL_edit_rate\tapo_ripk1_UT_vs_CTRL_edit_rate\ttad_fto_BT_vs_CTRL_edit_rate\ttad_fto_UT_vs_CTRL_edit_rate" > combined_edit_rates_v2.tsv

awk '
BEGIN {
    OFS="\t"
}

{
    sample = FILENAME
    if (sample ~ /BT/) bt = 1
    else if (sample ~ /UT/) bt = 0
    
    # Create cluster key
    key = $1 "\t" $2 "\t" $3
    
    # Store values by condition
    if (FILENAME ~ /apo_fto/) {
        if (bt) apo_fto_bt[key] += $5
        else apo_fto_ut[key] += $5
    }
    else if (FILENAME ~ /apo_ripk1/) {
        if (bt) apo_ripk1_bt[key] += $5
        else apo_ripk1_ut[key] += $5
    }
    else if (FILENAME ~ /tad_fto/) {
        if (bt) tad_fto_bt[key] += $5
        else tad_fto_ut[key] += $5
    }
    
    positions[key] = 1
}

END {
    # Print combined results
    for (pos in positions) {
        printf "%s\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",
            pos,
            (pos in apo_fto_bt) ? apo_fto_bt[pos] : 0,
            (pos in apo_fto_ut) ? apo_fto_ut[pos] : 0,
            (pos in apo_ripk1_bt) ? apo_ripk1_bt[pos] : 0,
            (pos in apo_ripk1_ut) ? apo_ripk1_ut[pos] : 0,
            (pos in tad_fto_bt) ? tad_fto_bt[pos] : 0,
            (pos in tad_fto_ut) ? tad_fto_ut[pos] : 0
    }
}' \
"$output_dir"/*/FTO/*_clean_v2.bed \
"$output_dir"/*/RIPK1/*_clean_v2.bed | \
sort -k1,1 -k2,2n >> combined_edit_rates_v2.tsv

echo "First few lines of output:"
head -n 5 combined_edit_rates_v2.tsv

# Summarize edit rates and calculate log2fc by cluster
awk -v OFS='\t' '
NR == 1 {
    for (i = 4; i <= NF; i++) {
        header[i] = $i
        col_count = NF
    }
    # Print header with all columns
    printf "chromosome\tstart\tend\tsize"
    for (i = 4; i <= col_count; i++) {
        printf "\t%s_count\t%s_sum", header[i], header[i]
    }
    printf "\tapo_fto_log2fc\tapo_ripk1_log2fc\ttad_fto_log2fc\tcombined_score\n"
    next
}

{
    cluster = $1
    clusters[cluster] = 1
    
    # Track start and end positions
    if (!(cluster in cluster_start)) {
        cluster_start[cluster] = $2  # First start position
    }
    cluster_end[cluster] = $3        # Last end position
    
    # Count non-zero edits and sum values for each condition
    for (i = 4; i <= col_count; i++) {
        key_sum = cluster "_" i "_sum"
        key_count = cluster "_" i "_count"
        if ($i > 0) {
            count[key_count]++
        }
        sum[key_sum] += $i
    }
}

END {
    n = 0
    for (c in clusters) {
        size = cluster_end[c] - cluster_start[c]
        
        # Calculate log2fc using sums
        apo_fto_fc = log2((sum[c "_4_sum"] + 1)/(sum[c "_5_sum"] + 1))
        apo_ripk1_fc = log2((sum[c "_6_sum"] + 1)/(sum[c "_7_sum"] + 1))
        tad_fto_fc = log2((sum[c "_8_sum"] + 1)/(sum[c "_9_sum"] + 1))
        
        # Calculate combined score
        combined_score = 0
        if (apo_fto_fc > 0 && apo_ripk1_fc > 0 && tad_fto_fc > 0) {
            combined_score = (apo_fto_fc + apo_ripk1_fc + tad_fto_fc) / 3
        }
        
        lines[++n] = sprintf("%s\t%d\t%d\t%d", c, cluster_start[c], cluster_end[c], size)
        for (j = 4; j <= col_count; j++) {
            key_sum = c "_" j "_sum"
            key_count = c "_" j "_count"
            lines[n] = lines[n] sprintf("\t%d\t%.6f", 
                  (key_count in count) ? count[key_count] : 0,
                  (key_sum in sum) ? sum[key_sum] : 0)
        }
        lines[n] = lines[n] sprintf("\t%.4f\t%.4f\t%.4f\t%.4f", 
               apo_fto_fc, apo_ripk1_fc, tad_fto_fc, combined_score)
        scores[n] = combined_score
    }
    
    # Sort by combined score (descending)
    for (i = 1; i <= n; i++) {
        for (j = i + 1; j <= n; j++) {
            if (scores[j] > scores[i]) {
                # Swap scores
                temp = scores[i]
                scores[i] = scores[j]
                scores[j] = temp
                # Swap lines
                temp = lines[i]
                lines[i] = lines[j]
                lines[j] = temp
            }
        }
    }
    
    # Print sorted results
    for (i = 1; i <= n; i++) {
        print lines[i]
    }
}

function log2(x) {
    return log(x)/log(2)
}' combined_edit_rates_v2.tsv > summed_edit_rates_v2.tsv

echo "Top results by combined score:"
head -n 5 summed_edit_rates_v2.tsv

genes_file="/Volumes/Lucky_me/9_STAMP_transcript/cluster_unique_genes.tsv"

awk -F'\t' -v OFS='\t' -v genes_file="$genes_file" '
FILENAME == genes_file {
    if (FNR > 1) {  # Skip header of genes file
        genes[$1] = $2
    }
    next
}

NR == 1 {
    print $0 "\tgenes"
    next
}

{   
    cluster = $1
    print $0, (genes[cluster] ? genes[cluster] : "NA")
}' "$genes_file" summed_edit_rates_v2.tsv > summed_edit_rates_with_genes_v2.tsv

echo -e "\nResults with genes:"
if [ -s summed_edit_rates_with_genes.tsv ]; then
    echo "First 5 lines of output:"
    head -n 5 summed_edit_rates_with_genes.tsv
    echo -e "\nTotal lines:"
    wc -l summed_edit_rates_with_genes.tsv
else
    echo "Error: Output file is empty"
fi
