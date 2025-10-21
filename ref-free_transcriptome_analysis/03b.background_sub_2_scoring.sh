#background subtraction B: Site-Level Correction
output_dir="/Volumes/Lucky_me/9_STAMP_transcript/06.transcript_analysis_v2"

combine_and_sum() {
    local input_dir=$1
    local output_file=$2
    local sample_type=$3  # "BT" or "UT"
    
    if ! ls ${input_dir}/*${sample_type}*_background.bed >/dev/null 2>&1; then
        echo "No ${sample_type} files found in ${input_dir}"
        return 1
    fi
    
    mkdir -p "$(dirname "$output_file")"
    
    cat ${input_dir}/*${sample_type}*_background.bed | sort -k1,1 -k2,2n > ${input_dir}/temp_sorted_${sample_type}.bed
    
    # Calculate sum of edit rates for each position
    awk 'BEGIN {OFS="\t"}
    {
        # Ensure fields are properly separated
        gsub(/[[:space:]]+/, "\t")
        
        # Skip malformed lines
        if (NF != 6) {
            print "Skipping malformed line:", $0 > "/dev/stderr"
            next
        }
        
        # Create key with explicit tab separation
        key = $1 "\t" $2 "\t" $3 "\t" $4 "\t" $6
        sums[key] += $5
    }
    END {
        for (key in sums) {
            split(key, fields, "\t")
            printf "%s\t%s\t%s\t%s\t%.6f\t%s\n", 
                fields[1], fields[2], fields[3], fields[4], sums[key], fields[5]
        }
    }' ${input_dir}/temp_sorted_${sample_type}.bed | sort -k1,1 -k2,2n > "$output_file"
    
    rm ${input_dir}/temp_sorted_${sample_type}.bed
    
    echo "Created summed file: ${output_file}"
}

echo "Processing Apo-FTO controls..."
combine_and_sum "$output_dir/apo_fto/control" "$output_dir/apo_fto/final/apo_fto_ctrl_BT_summed.bed" "BT"
combine_and_sum "$output_dir/apo_fto/control" "$output_dir/apo_fto/final/apo_fto_ctrl_UT_summed.bed" "UT"

echo "Processing Apo-RIPK1 controls..."
combine_and_sum "$output_dir/apo_ripk1/control" "$output_dir/apo_ripk1/final/apo_ripk1_ctrl_BT_summed.bed" "BT"
combine_and_sum "$output_dir/apo_ripk1/control" "$output_dir/apo_ripk1/final/apo_ripk1_ctrl_UT_summed.bed" "UT"

echo "Processing TAD-FTO controls..."
combine_and_sum "$output_dir/tad_fto/control" "$output_dir/tad_fto/final/tad_fto_ctrl_BT_summed.bed" "BT"
combine_and_sum "$output_dir/tad_fto/control" "$output_dir/tad_fto/final/tad_fto_ctrl_UT_summed.bed" "UT"


input_dir="/Volumes/Lucky_me/9_STAMP_transcript/05.transcript_analysis"
# Function to process sample edits with control subtraction
process_edits_with_subtraction() {
    local sample=$1
    local type=$2
    local gene=$3
    local sample_treatment
    local control_file
    
    if [[ $sample == *"BT"* ]]; then
        sample_treatment="BT"
    else
        sample_treatment="UT"
    fi
    
    # Set appropriate control file based on type and treatment
    control_file="$output_dir/${type}/final/${type}_ctrl_${sample_treatment}_summed.bed"
    
    echo "Processing ${type} ${sample_treatment} sample: $sample with control subtraction"
    echo "Using control file: $control_file"
    mkdir -p "$output_dir/${type}/${gene}/stats"
    local stats_file="$output_dir/${type}/${gene}/stats/${sample}_filtering_stats.txt"
    
    echo "Processing ${type} sample: $sample with control subtraction"
    mkdir -p "$output_dir/${type}/${gene}/stats"
    local stats_file="$output_dir/${type}/${gene}/stats/${sample}_filtering_stats.txt"
    
    local initial_count=$(wc -l < "$input_dir/${type}/${gene}/${sample}_edits_trim.bed")
    echo "Initial sites: $initial_count" > "$stats_file"
    
    bedtools intersect -a "$input_dir/${type}/${gene}/${sample}_edits_trim.bed" \
                      -b "$control_file" \
                      -wb | \
    awk 'BEGIN {OFS="\t"
           print "=== Filtering Statistics ===" > "'$stats_file'"
    }
    {
        # $5 is sample edit rate, $11 is control edit rate
        adjusted_rate = $5 - $13
        
        # Count different filtering cases
        if ($5 == 1) {
            perfect_edits++
            print $1, $2, $3, $4, $5, $13 > "'$output_dir/${type}/${gene}/stats/${sample}_perfect_edits.txt'"
        }
        else if (adjusted_rate <= 0) {
            removed_by_control++
            print $1, $2, $3, $4, $5, $13 > "'$output_dir/${type}/${gene}/stats/${sample}_removed_by_control.txt'"
        }
        else {
            # Output adjusted sites
            print $1, $2, $3, $4, adjusted_rate, $6, $7, $8 > "'$output_dir/${type}/${gene}/${sample}_clean_subtracted.bed'"
            print $1, $2, $3, $4, $5, $13, adjusted_rate > "'$output_dir/${type}/${gene}/stats/${sample}_adjusted_sites.txt'"
            kept_after_adjustment++
        }
    }
    END {
        print "\nFiltering Results:" >> "'$stats_file'"
        print "Perfect edit sites (rate=1):", perfect_edits >> "'$stats_file'"
        print "Sites removed by control:", removed_by_control >> "'$stats_file'"
        print "Sites kept after adjustment:", kept_after_adjustment >> "'$stats_file'"
    }' 
    
    # Handle non-overlapping regions
    bedtools intersect -a "$input_dir/${type}/${gene}/${sample}_edits_trim.bed" \
                      -b "$control_file" \
                      -v | \
    awk '$5 != 1' >> "$output_dir/${type}/${gene}/${sample}_clean_subtracted.bed"
    
    # Count non-overlapping sites
    local non_overlapping=$(bedtools intersect -a "$input_dir/${type}/${gene}/${sample}_edits_trim.bed" \
                                             -b "$control_file" -v | wc -l)
    echo "Non-overlapping sites: $non_overlapping" >> "$stats_file"
    
    # Sort and create final output
    sort -k1,1 -k2,2n "$output_dir/${type}/${gene}/${sample}_clean_subtracted.bed" \
        > "$output_dir/${type}/${gene}/${sample}_clean_v2.bed"
    
    # Final count
    local final_count=$(wc -l < "$output_dir/${type}/${gene}/${sample}_clean_v2.bed")
    echo -e "\nFinal sites after all filtering: $final_count" >> "$stats_file"
    
    echo "  Found $final_count sites after processing"
    echo "  Detailed statistics saved to: $stats_file"
    
    # Cleanup
    rm "$output_dir/${type}/${gene}/${sample}_clean_subtracted.bed"
}

# Process each enzyme group with appropriate control files
echo "Processing Apo-FTO samples..."
for sample in "${APO_FTO_SAMPLES[@]}"; do
    process_edits_with_subtraction "$sample" "apo_fto" "FTO"
done

echo "Processing Apo-RIPK1 samples..."
for sample in "${APO_RIPK1_SAMPLES[@]}"; do
    process_edits_with_subtraction "$sample" "apo_ripk1" "RIPK1"
done

echo "Processing TAD-FTO samples..."
for sample in "${TAD_FTO_SAMPLES[@]}"; do
    process_edits_with_subtraction "$sample" "tad_fto" "FTO"
done

output_dir="/Volumes/Lucky_me/9_STAMP_transcript/06.transcript_analysis_v2"

awk -v OFS='\t' '
BEGIN {
    OFS="\t"
}

{
    sample = FILENAME
    if (sample ~ /BT/) bt = 1
    else if (sample ~ /UT/) bt = 0
    
    if (match(sample, /R[0-9]+(_2)?/)) {
        rep_num = substr(sample, RSTART+1, RLENGTH-1)
        is_rep2 = match(sample, /_2$/)
    }
    
    # Create cluster key
    key = $1 "\t" $2 "\t" $3
    
    # Store values by condition
    if (FILENAME ~ /apo_ripk1/) {
        if (FILENAME ~ /control/) {
            if (bt) {
                if (rep_num ~ /^9/) {          # Changed from (9|10) to just 9
                    apo_ripk1_ctrl_bt_1[key] += $5 # R9 and R9_2 go to bt_1
                } else if (rep_num ~ /^10/) {  # Handle R10 separately
                    apo_ripk1_ctrl_bt_2[key] += $5 # R10 and R10_2 go to bt_2
                }
            } else {
                if (rep_num ~ /^1/) {
                    apo_ripk1_ctrl_ut_1[key] += $5  # R1 and R1_2
                } else if (rep_num ~ /^2/) {
                    apo_ripk1_ctrl_ut_2[key] += $5  # R2 and R2_2
                }
            }
        } else {
            if (bt) apo_ripk1_bt[key] += $5
            else apo_ripk1_ut[key] += $5
        }
    } else if (FILENAME ~ /apo_fto/) {
        if (FILENAME ~ /control/) {
            if (bt) apo_fto_ctrl_bt[key] += $5
            else apo_fto_ctrl_ut[key] += $5
        } else {
            if (bt) apo_fto_bt[key] += $5
            else apo_fto_ut[key] += $5
        }
    } else if (FILENAME ~ /tad_fto/) {
        if (FILENAME ~ /control/) {
            if (bt) tad_fto_ctrl_bt[key] += $5
            else tad_fto_ctrl_ut[key] += $5
        } else {
            if (bt) tad_fto_bt[key] += $5
            else tad_fto_ut[key] += $5
        }
    }
    positions[key] = 1
}

END {
    print "chromosome\tstart\tend\t" \
          "apo_fto_bt\tapo_fto_ut\t" \
          "apo_ripk1_bt\tapo_ripk1_ut\t" \
          "tad_fto_bt\ttad_fto_ut\t" \
          "apo_fto_ctrl_bt\tapo_fto_ctrl_ut\t" \
          "apo_ripk1_ctrl_bt_1\tapo_ripk1_ctrl_ut_1\t" \
          "apo_ripk1_ctrl_bt_2\tapo_ripk1_ctrl_ut_2\t" \
          "tad_fto_ctrl_bt\ttad_fto_ctrl_ut"
    
    for (pos in positions) {
        printf "%s\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",
            pos,
            (pos in apo_fto_bt) ? apo_fto_bt[pos] : 0,
            (pos in apo_fto_ut) ? apo_fto_ut[pos] : 0,
            (pos in apo_ripk1_bt) ? apo_ripk1_bt[pos] : 0,
            (pos in apo_ripk1_ut) ? apo_ripk1_ut[pos] : 0,
            (pos in tad_fto_bt) ? tad_fto_bt[pos] : 0,
            (pos in tad_fto_ut) ? tad_fto_ut[pos] : 0,
            (pos in apo_fto_ctrl_bt) ? apo_fto_ctrl_bt[pos] : 0,
            (pos in apo_fto_ctrl_ut) ? apo_fto_ctrl_ut[pos] : 0,
            (pos in apo_ripk1_ctrl_bt_1) ? apo_ripk1_ctrl_bt_1[pos] : 0,
            (pos in apo_ripk1_ctrl_ut_1) ? apo_ripk1_ctrl_ut_1[pos] : 0,
            (pos in apo_ripk1_ctrl_bt_2) ? apo_ripk1_ctrl_bt_2[pos] : 0,
            (pos in apo_ripk1_ctrl_ut_2) ? apo_ripk1_ctrl_ut_2[pos] : 0,
            (pos in tad_fto_ctrl_bt) ? tad_fto_ctrl_bt[pos] : 0,
            (pos in tad_fto_ctrl_ut) ? tad_fto_ctrl_ut[pos] : 0
    }
}' \
"$output_dir"/*/control/*_background.bed \
"$output_dir"/*/FTO/*_clean_v2.bed \
"$output_dir"/*/RIPK1/*_clean_v2.bed | \
sort -k1,1 -k2,2n > combined_edit_rates_with_controls.tsv

awk -v OFS='\t' '
NR == 1 {
    for (i = 4; i <= NF; i++) {
        header[i] = $i
        printf "Column %d: %s\n", i, $i > "/dev/stderr"
        col_count = NF
    }
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
    # Output results (stored in array for sorting)
    n = 0
    for (c in clusters) {
        size = cluster_end[c] - cluster_start[c]

        apo_fto_fc = log2((sum[c "_4_sum"] + 1)/(sum[c "_5_sum"] + 1))
        apo_ripk1_fc = log2((sum[c "_6_sum"] + 1)/(sum[c "_7_sum"] + 1))
        tad_fto_fc = log2((sum[c "_8_sum"] + 1)/(sum[c "_9_sum"] + 1))

        combined_score = 0
        if ((apo_fto_fc > 0 && apo_ripk1_fc > 0 && tad_fto_fc > 0) ||
            (apo_ripk1_fc > 0 && (apo_fto_fc > 0 || tad_fto_fc > 0)))
        {
            # Primary score based on fold change magnitudes
            fc_score = apo_fto_fc * apo_ripk1_fc * tad_fto_fc

            # Calculate differences from RIPK1
            apo_fto_ripk1_diff = abs(apo_fto_fc - apo_ripk1_fc)
            ripk1_tad_diff = abs(apo_ripk1_fc - tad_fto_fc)

            # Weight by similarity to RIPK1
            sim_score = 1 / (1 + apo_fto_ripk1_diff + ripk1_tad_diff)

            # Set weight
            weight = 1
            if (apo_fto_fc > 0 && apo_ripk1_fc > 0 && tad_fto_fc > 0) {
                weight = 1
            } else if (apo_ripk1_fc > 0 && (apo_fto_fc > 0 || tad_fto_fc > 0)) {
                weight = 1
            }

            # Final combined score
            combined_score = fc_score * sim_score * weight
        }

        # Store line and score for sorting
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
}

function abs(x) {
    return x < 0 ? -x : x
}' combined_edit_rates_with_controls.tsv > summed_edit_rates_trim_with_controls.tsv

echo "Top results by combined score:"
head -n 5 summed_edit_rates_trim_with_controls.tsv

awk -v OFS='\t' '
NR == 1 {
    for (i = 4; i <= NF; i++) {
        header[i] = $i
        printf "Column %d: %s\n", i, $i > "/dev/stderr"
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
    # Output results (stored in array for sorting)
    n = 0
    for (c in clusters) {
        size = cluster_end[c] - cluster_start[c]

        # Calculate log2fc using sums
        # Get fold change columns by name
        apo_fto_fc = log2((sum[c "_4_sum"] + 1)/(sum[c "_5_sum"] + 1))
        apo_ripk1_fc = log2((sum[c "_6_sum"] + 1)/(sum[c "_7_sum"] + 1))
        tad_fto_fc = log2((sum[c "_8_sum"] + 1)/(sum[c "_9_sum"] + 1))

        combined_score = 0
        if (apo_ripk1_fc > 0 && apo_fto_fc > 0) {
            # Simplified score using only apo_fto and apo_ripk1
            fc_score = apo_fto_fc * apo_ripk1_fc
            
            # Calculate similarity between apo_fto and apo_ripk1
            apo_fto_ripk1_diff = abs(apo_fto_fc - apo_ripk1_fc)
            sim_score = 1 / (1 + apo_fto_ripk1_diff)

            weight = 1
            if (apo_fto_fc >= 1 && apo_ripk1_fc >= 1) {
                weight = 1
            }

            # Final combined score
            combined_score = fc_score * sim_score * weight
        }

        # Store line and score for sorting
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
}

function abs(x) {
    return x < 0 ? -x : x
}' combined_edit_rates_with_controls.tsv > summed_edit_rates_trim_with_controls_apo.tsv

echo "Top results by combined score:"
head -n 5 summed_edit_rates_trim_with_controls_apo.tsv

awk -v OFS='\t' '
NR == 1 {
    for (i = 4; i <= NF; i++) {
        header[i] = $i
        printf "Column %d: %s\n", i, $i > "/dev/stderr"
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

        combined_score = 0
        if (apo_ripk1_fc > 0 && tad_fto_fc > 0) {
            # Simplified score using only tad_fto and apo_ripk1
            fc_score = tad_fto_fc * apo_ripk1_fc
            
            # Calculate similarity between tad_fto and apo_ripk1
            ripk1_tad_diff = abs(apo_ripk1_fc - tad_fto_fc)
            sim_score = 1 / (1 + ripk1_tad_diff)

            # Set weight based on magnitude
            weight = 1
            if (tad_fto_fc >= 1 && apo_ripk1_fc >= 1) {
                weight = 1
            }

            # Final combined score
            combined_score = fc_score * sim_score * weight
        }

        # Store line and score for sorting
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
}

function abs(x) {
    return x < 0 ? -x : x
}' combined_edit_rates_with_controls.tsv > summed_edit_rates_trim_with_controls_tad.tsv

echo "Top results by combined score:"
head -n 5 summed_edit_rates_trim_with_controls_tad.tsv

#annotate with genes name
genes_file="/Volumes/Lucky_me/9_STAMP_transcript/cluster_unique_genes.tsv"

awk -F'\t' -v OFS='\t' -v genes_file="$genes_file" '
FILENAME == genes_file {
    if (FNR > 1) {  # Skip header of genes file
        genes[$1] = $2
    }
    next
}

NR == 1 {
    # Print header with genes column
    print $0 "\tgenes"
    next
}

{   
    cluster = $1
    print $0, (genes[cluster] ? genes[cluster] : "NA")
}' "$genes_file" summed_edit_rates_trim_with_controls_tad.tsv > summed_edit_rates_with_genes_trim_with_controls_tad.tsv

echo -e "\nResults with genes:"
if [ -s summed_edit_rates_with_genes_trim_with_controls_tad.tsv ]; then
    echo "First 5 lines of output:"
    head -n 5 summed_edit_rates_with_genes_trim_with_controls_tad.tsv
    echo -e "\nTotal lines:"
    wc -l summed_edit_rates_with_genes_trim_with_controls_tad.tsv
else
    echo "Error: Output file is empty"
fi

awk -v OFS='\t' '
FILENAME == ARGV[1] {
    seq_lengths[$1] = $2
    next
}

NR == 1 {
    print $0
    next
}

{
    cluster = $1
    if (cluster in seq_lengths) {
        # Replace start and end positions with 1 and sequence length
        $2 = 1
        $3 = seq_lengths[cluster]
        $4 = seq_lengths[cluster]
    }
    print $0
}' /Volumes/Lucky_me/9_STAMP_transcript/transcriptome_clean.fa.fai \
   summed_edit_rates_with_genes_trim_with_controls_tad.tsv > summed_edit_rates_with_genes_trim_with_controls_tad_cl.tsv 

