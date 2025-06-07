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

# Parameters
MIN_COVERAGE=10

mkdir -p /home/onggf/data/7_STAMP/data/3Mar_data/SKOV3/03.comparative_analysis/{apo_fto,apo_ripk1,tad_fto}/{control,FTO,RIPK1,final}

output_dir="/home/onggf/data/7_STAMP/data/3Mar_data/SKOV3/03.comparative_analysis"
input_dir_tad="/home/onggf/data/7_STAMP/data/3Mar_data/SKOV3/01.Tad_analysis"
input_dir_apo="/home/onggf/data/7_STAMP/data/3Mar_data/SKOV3/02.Apo_analysis"

# Check for input files presence
for sample in "${APO_FTO_CTRL[@]}" "${APO_RIPK1_CTRL[@]}" "${APO_FTO_SAMPLES[@]}" "${APO_RIPK1_SAMPLES[@]}"; do
    if [ ! -f "${input_dir_apo}/${sample}/${sample}.GRCh38_Apo_CT.pileup2var.flt.txt" ]; then
        echo "Error: Input file ${input_dir_apo}/${sample}/${sample}.GRCh38_Apo_CT.pileup2var.flt.txt not found"
    fi
done

for sample in "${TAD_FTO_CTRL[@]}" "${TAD_FTO_SAMPLES[@]}"; do
    if [ ! -f "${input_dir_tad}/${sample}/${sample}.GRCh38_Tad_AG.pileup2var.flt.txt" ]; then
        echo "Error: Input file ${input_dir_tad}/${sample}/${sample}.GRCh38_Tad_AG.pileup2var.flt.txt not found"
    fi
done

# Step 1: Process background variants
process_background() {
    local sample=$1
    local output_loc=$2
    local conv_type=$3
    
    if [ "$conv_type" = "apo" ]; then
        # C>T background for Apo
        awk -v min_cov=$MIN_COVERAGE '
        NR>1 && $4=="C" && $5>=min_cov && $7>0 {
            print $1"\t"$2-1"\t"$2"\t"$4"\t"$7/$5"\t"$3
        }' "$input_dir_apo/${sample}/${sample}.GRCh38_Apo_CT.pileup2var.flt.txt" > "$output_loc/${sample}_background.bed"
    else
        # A>G background for Tad
        awk -v min_cov=$MIN_COVERAGE '
        NR>1 && $4=="A" && $5>=min_cov && $9>0 {
            print $1"\t"$2-1"\t"$2"\t"$4"\t"$9/$5"\t"$3
        }' "$input_dir_tad/${sample}/${sample}.GRCh38_Tad_AG.pileup2var.flt.txt" > "$output_loc/${sample}_background.bed"
    fi
}

echo "Processing Apo-fto ctrl..."
for sample in "${APO_FTO_CTRL[@]}"; do
    process_background "$sample" "$output_dir/apo_fto/control" "apo"
done

echo "Processing Apo-ripk1 ctrl..."
for sample in "${APO_RIPK1_CTRL[@]}"; do
    process_background "$sample" "$output_dir/apo_ripk1/control" "apo"
done

echo "Processing Tad-fto ctrl..."
for sample in "${TAD_FTO_CTRL[@]}"; do
    process_background "$sample" "$output_dir/tad_fto/control" "tad"
done

# Merge background variants into single bed files per group
cat "$output_dir/apo_fto/control/"*_background.bed | sort -k1,1 -k2,2n | bedtools merge -i stdin > "$output_dir/apo_fto/control/SKOV3_Apo-Cnt-FTO_background_variants.bed"
cat "$output_dir/apo_ripk1/control/"*_background.bed | sort -k1,1 -k2,2n | bedtools merge -i stdin > "$output_dir/apo_ripk1/control/SKOV3_Apo-Cnt-RIPK1_background_variants.bed"
cat "$output_dir/tad_fto/control/"*_background.bed | sort -k1,1 -k2,2n | bedtools merge -i stdin > "$output_dir/tad_fto/control/SKOV3_Tad-Cnt-FTO_background_variants.bed"

# Step 2: Process sample variants with respective base changes
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
        }' "$input_dir_apo/${sample}/${sample}.GRCh38_Apo_CT.pileup2var.flt.txt" > "$output_loc/${sample}_edits.bed"
    else
        # A>G conversion for Tad
        awk -v min_cov=$MIN_COVERAGE '
        NR>1 && $4=="A" && $5>=min_cov && $9>0 {
            edit_rate = $9/$5
            print $1"\t"$2-1"\t"$2"\t"$4"\t"edit_rate"\t"$5"\t"$9"\t"$3
        }' "$input_dir_tad/${sample}/${sample}.GRCh38_Tad_AG.pileup2var.flt.txt" > "$output_loc/${sample}_edits.bed"
    fi
}

for sample in "${APO_FTO_SAMPLES[@]}"; do
    process_sample "$sample" "$output_dir/apo_fto/FTO" "apo"
done

for sample in "${APO_RIPK1_SAMPLES[@]}"; do
    process_sample "$sample" "$output_dir/apo_ripk1/RIPK1" "apo"
done

for sample in "${TAD_FTO_SAMPLES[@]}"; do
    process_sample "$sample" "$output_dir/tad_fto/FTO" "tad"
done

# Step 3: Remove background variants and create cleaned bedgraph
process_group() {
    local sample=$1
    local type=$2
    local gene=$3
    
    echo "Processing ${type} sample: $sample"
    
    local bg_bed=""
    if [[ "$type" == "apo_fto" ]]; then
        bg_bed="$output_dir/${type}/control/SKOV3_Apo-Cnt-FTO_background_variants.bed"
    elif [[ "$type" == "apo_ripk1" ]]; then
        bg_bed="$output_dir/${type}/control/SKOV3_Apo-Cnt-RIPK1_background_variants.bed"
    elif [[ "$type" == "tad_fto" ]]; then
        bg_bed="$output_dir/${type}/control/SKOV3_Tad-Cnt-FTO_background_variants.bed"
    fi

    local sample_bed="$output_dir/${type}/${gene}/${sample}_edits.bed"
    local final_bed="$output_dir/${type}/final/${sample}_clean.bed"

    bedtools subtract -a "$sample_bed" -b "$bg_bed" > "$final_bed"

    awk '{ print $1"\t"$2"\t"$3"\t"$5 }' "$final_bed" > "${final_bed%.bed}.bedgraph"
}

for sample in "${APO_FTO_SAMPLES[@]}"; do
    process_group "$sample" "apo_fto" "FTO"
done

for sample in "${APO_RIPK1_SAMPLES[@]}"; do
    process_group "$sample" "apo_ripk1" "RIPK1"
done

for sample in "${TAD_FTO_SAMPLES[@]}"; do
    process_group "$sample" "tad_fto" "FTO"
done

echo "All processing complete."
