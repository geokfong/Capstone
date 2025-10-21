#background subtraction C: Cluster-Level Correction
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

output_dir="/Volumes/Lucky_me/9_STAMP_transcript/07.transcript_analysis_v3"
input_dir_tad="/Volumes/Lucky_me/9_STAMP_transcript/pileup2var"
input_dir_apo="/Volumes/Lucky_me/9_STAMP_transcript/pileup2var"
mkdir -p "$output_dir/apo_fto/FTO"
mkdir -p "$output_dir/apo_ripk1/RIPK1"
mkdir -p "$output_dir/tad_fto/FTO"

MIN_COVERAGE=10

process_sample() {
    local sample=$1
    local output_loc=$2
    local conv_type=$3
    
    if [ "$conv_type" = "apo" ]; then
        awk -v min_cov=$MIN_COVERAGE '
        FILENAME == ARGV[1] {
            lengths[$1] = $2
            next
        }
        FILENAME == ARGV[2] && NR>1 && $4=="C" && $5>=min_cov && $7>0 {
            if ($2 > 100 && $2 < (lengths[$1] - 100)) {
                edit_rate = $7/$5
                print $1"\t"$2-1"\t"$2"\t"$4"\t"edit_rate"\t"$5"\t"$7"\t"$3
            }
        }' /Volumes/Lucky_me/9_STAMP_transcript/transcriptome_clean.fa.fai \
           $input_dir_apo/${sample}.Transcript_CT.pileup2var.flt.txt > $output_loc/${sample}_edits_trim.bed
    else
        awk -v min_cov=$MIN_COVERAGE '
        FILENAME == ARGV[1] {
            lengths[$1] = $2
            next
        }
        FILENAME == ARGV[2] && NR>1 && $4=="A" && $5>=min_cov && $9>0 {
            if ($2 > 100 && $2 < (lengths[$1] - 100)) {
                edit_rate = $9/$5
                print $1"\t"$2-1"\t"$2"\t"$4"\t"edit_rate"\t"$5"\t"$9"\t"$3
            }
        }' /Volumes/Lucky_me/9_STAMP_transcript/transcriptome_clean.fa.fai \
           $input_dir_tad/${sample}.Transcript_AG.pileup2var.flt.txt > $output_loc/${sample}_edits_trim.bed
    fi
}

for sample in "${APO_FTO_SAMPLES[@]}"; do
    echo "Processing APO FTO sample: $sample"
    process_sample "$sample" "$output_dir/apo_fto/FTO" "apo"
done

for sample in "${APO_RIPK1_SAMPLES[@]}"; do
    echo "Processing APO RIPK1 sample: $sample"
    process_sample "$sample" "$output_dir/apo_ripk1/RIPK1" "apo"
done

for sample in "${TAD_FTO_SAMPLES[@]}"; do
    echo "Processing TAD FTO sample: $sample"
    process_sample "$sample" "$output_dir/tad_fto/FTO" "tad"
done

mkdir -p "$output_dir/apo_fto/control"
mkdir -p "$output_dir/apo_ripk1/control"
mkdir -p "$output_dir/tad_fto/control"  

process_background() {
    local sample=$1
    local output_loc=$2
    local conv_type=$3
    
    if [ "$conv_type" = "apo" ]; then
        # C>T background for Apo with length filtering
        awk -v min_cov=$MIN_COVERAGE '
        # First pass to store sequence lengths from fai
        FILENAME == ARGV[1] {
            lengths[$1] = $2
            next
        }
        # Second pass for pileup2var processing
        FILENAME == ARGV[2] && NR>1 && $4=="C" && $5>=min_cov && $7>0 {
            # Only process if position is beyond first/last 100bp
            if ($2 > 100 && $2 < (lengths[$1] - 100)) {
                edit_rate = $7/$5
                print $1"\t"$2-1"\t"$2"\t"$4"\t"edit_rate"\t"$3
            }
        }' /Volumes/Lucky_me/9_STAMP_transcript/transcriptome_clean.fa.fai \
           $input_dir_apo/${sample}.Transcript_CT.pileup2var.flt.txt > $output_loc/${sample}_background_trim.bed
    else
        # A>G background for Tad with length filtering
        awk -v min_cov=$MIN_COVERAGE '
        # First pass to store sequence lengths from fai
        FILENAME == ARGV[1] {
            lengths[$1] = $2
            next
        }
        # Second pass for pileup2var processing
        FILENAME == ARGV[2] && NR>1 && $4=="A" && $5>=min_cov && $9>0 {
            # Only process if position is beyond first/last 100bp
            if ($2 > 100 && $2 < (lengths[$1] - 100)) {
                edit_rate = $9/$5
                print $1"\t"$2-1"\t"$2"\t"$4"\t"edit_rate"\t"$3
            }
        }' /Volumes/Lucky_me/9_STAMP_transcript/transcriptome_clean.fa.fai \
           $input_dir_tad/${sample}.Transcript_AG.pileup2var.flt.txt > $output_loc/${sample}_background_trim.bed
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

awk '
BEGIN {
    OFS = "\t"
}

FNR == 1 {
    printed = 0
}

{
    sample = FILENAME

    if (sample ~ /BT/) bt = 1
    else if (sample ~ /UT/) bt = 0

    if (!printed) {
        print "Processing file:", FILENAME, "| BT =", bt > "/dev/stderr"
        printed = 1
    }

    key = $1 "\t" $2 "\t" $3

    if (sample ~ /apo_ripk1_2/) {
        if (sample ~ /control/) {
            if (bt) apo_ripk1_ctrl_bt_2[key] += $5
            else apo_ripk1_ctrl_ut_2[key] += $5
        } else if (sample ~ /RIPK1/) {
            if (bt) apo_ripk1_bt_2[key] += $5
            else apo_ripk1_ut_2[key] += $5
        }
    } else if (sample ~ /apo_ripk1/) {
        if (sample ~ /control/) {
            if (bt) apo_ripk1_ctrl_bt_1[key] += $5
            else apo_ripk1_ctrl_ut_1[key] += $5
        } else if (sample ~ /RIPK1/) {
            if (bt) apo_ripk1_bt_1[key] += $5
            else apo_ripk1_ut_1[key] += $5
        }
    } else if (sample ~ /apo_fto/) {
        if (sample ~ /control/) {
            if (bt) apo_fto_ctrl_bt[key] += $5
            else apo_fto_ctrl_ut[key] += $5
        } else {
            if (bt) apo_fto_bt[key] += $5
            else apo_fto_ut[key] += $5
        }
    } else if (sample ~ /tad_fto/) {
        if (sample ~ /control/) {
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
          "apo_ripk1_bt_1\tapo_ripk1_bt_2\t" \
          "apo_ripk1_ut_1\tapo_ripk1_ut_2\t" \
          "tad_fto_bt\ttad_fto_ut\t" \
          "apo_fto_ctrl_bt\tapo_fto_ctrl_ut\t" \
          "apo_ripk1_ctrl_bt_1\tapo_ripk1_ctrl_ut_1\t" \
          "apo_ripk1_ctrl_bt_2\tapo_ripk1_ctrl_ut_2\t" \
          "tad_fto_ctrl_bt\ttad_fto_ctrl_ut"

    for (pos in positions) {
        printf "%s\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",
            pos,
            (pos in apo_fto_bt) ? apo_fto_bt[pos] : 0,
            (pos in apo_fto_ut) ? apo_fto_ut[pos] : 0,
            (pos in apo_ripk1_bt_1) ? apo_ripk1_bt_1[pos] : 0,
            (pos in apo_ripk1_bt_2) ? apo_ripk1_bt_2[pos] : 0,
            (pos in apo_ripk1_ut_1) ? apo_ripk1_ut_1[pos] : 0,
            (pos in apo_ripk1_ut_2) ? apo_ripk1_ut_2[pos] : 0,
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
}
' "$output_dir"/*/control/*_background_trim.bed "$output_dir"/*/FTO/*_edits_trim.bed "$output_dir"/*/RIPK1/*_edits_trim.bed \
| sort -k1,1 -k2,2n > combined_edit_rates_with_controls.tsv

awk -v OFS='\t' '
NR == 1 {
    print "chromosome\tstart\tend\t" \
          "apo_fto_bt_sum\tapo_fto_ut_sum\t" \
          "apo_ripk1_bt_1_sum\tapo_ripk1_ut_1_sum\t" \
          "apo_ripk1_bt_2_sum\tapo_ripk1_ut_2_sum\t" \
          "tad_fto_bt_sum\ttad_fto_ut_sum\t" \
          "apo_fto_ctrl_bt_sum\tapo_fto_ctrl_ut_sum\t" \
          "apo_ripk1_ctrl_bt_1_sum\tapo_ripk1_ctrl_ut_1_sum\t" \
          "apo_ripk1_ctrl_bt_2_sum\tapo_ripk1_ctrl_ut_2_sum\t" \
          "tad_fto_ctrl_bt_sum\ttad_fto_ctrl_ut_sum\t" \
          "apo_fto_bt_diff\tapo_fto_ut_diff\t" \
          "apo_ripk1_bt_1_diff\tapo_ripk1_ut_1_diff\t" \
          "apo_ripk1_bt_2_diff\tapo_ripk1_ut_2_diff\t" \
          "tad_fto_bt_diff\ttad_fto_ut_diff\t" \
          "log2fc_apo_fto\tlog2fc_ripk1_1\tlog2fc_ripk1_2\tlog2fc_tad_fto\t" \
          "combined_score"
    next
}

{
    cluster = $1
    clusters[cluster] = 1

    # Track start and end positions
    if (!(cluster in cluster_start)) {
        cluster_start[cluster] = $2
    }
    cluster_end[cluster] = $3

    # Sum values for each condition
    apo_fto_bt_sum[cluster] += $4
    apo_fto_ut_sum[cluster] += $5
    apo_ripk1_bt_1_sum[cluster] += $6
    apo_ripk1_ut_1_sum[cluster] += $7
    apo_ripk1_bt_2_sum[cluster] += $8
    apo_ripk1_ut_2_sum[cluster] += $9
    tad_fto_bt_sum[cluster] += $10
    tad_fto_ut_sum[cluster] += $11
    apo_fto_ctrl_bt_sum[cluster] += $12
    apo_fto_ctrl_ut_sum[cluster] += $13
    apo_ripk1_ctrl_bt_1_sum[cluster] += $14
    apo_ripk1_ctrl_ut_1_sum[cluster] += $15
    apo_ripk1_ctrl_bt_2_sum[cluster] += $16
    apo_ripk1_ctrl_ut_2_sum[cluster] += $17
    tad_fto_ctrl_bt_sum[cluster] += $18
    tad_fto_ctrl_ut_sum[cluster] += $19
}

END {
    for (cluster in clusters) {
        # Calculate differences using summed values
        apo_fto_bt_diff = (apo_fto_bt_sum[cluster] - apo_fto_ctrl_bt_sum[cluster] > 0) ? \
            apo_fto_bt_sum[cluster] - apo_fto_ctrl_bt_sum[cluster] : 0
        apo_fto_ut_diff = (apo_fto_ut_sum[cluster] - apo_fto_ctrl_ut_sum[cluster] > 0) ? \
            apo_fto_ut_sum[cluster] - apo_fto_ctrl_ut_sum[cluster] : 0
        
        apo_ripk1_bt_1_diff = (apo_ripk1_bt_1_sum[cluster] - apo_ripk1_ctrl_bt_1_sum[cluster] > 0) ? \
            apo_ripk1_bt_1_sum[cluster] - apo_ripk1_ctrl_bt_1_sum[cluster] : 0
        apo_ripk1_ut_1_diff = (apo_ripk1_ut_1_sum[cluster] - apo_ripk1_ctrl_ut_1_sum[cluster] > 0) ? \
            apo_ripk1_ut_1_sum[cluster] - apo_ripk1_ctrl_ut_1_sum[cluster] : 0
        
        apo_ripk1_bt_2_diff = (apo_ripk1_bt_2_sum[cluster] - apo_ripk1_ctrl_bt_2_sum[cluster] > 0) ? \
            apo_ripk1_bt_2_sum[cluster] - apo_ripk1_ctrl_bt_2_sum[cluster] : 0
        apo_ripk1_ut_2_diff = (apo_ripk1_ut_2_sum[cluster] - apo_ripk1_ctrl_ut_2_sum[cluster] > 0) ? \
            apo_ripk1_ut_2_sum[cluster] - apo_ripk1_ctrl_ut_2_sum[cluster] : 0
        
        tad_fto_bt_diff = (tad_fto_bt_sum[cluster] - tad_fto_ctrl_bt_sum[cluster] > 0) ? \
            tad_fto_bt_sum[cluster] - tad_fto_ctrl_bt_sum[cluster] : 0
        tad_fto_ut_diff = (tad_fto_ut_sum[cluster] - tad_fto_ctrl_ut_sum[cluster] > 0) ? \
            tad_fto_ut_sum[cluster] - tad_fto_ctrl_ut_sum[cluster] : 0

        # Calculate log2 fold changes using summed values
        log2fc_apo_fto = log2((apo_fto_bt_diff + 1)/(apo_fto_ut_diff + 1))
        log2fc_ripk1_1 = log2((apo_ripk1_bt_1_diff + 1)/(apo_ripk1_ut_1_diff + 1))
        log2fc_ripk1_2 = log2((apo_ripk1_bt_2_diff + 1)/(apo_ripk1_ut_2_diff + 1))
        log2fc_tad_fto = log2((tad_fto_bt_diff + 1)/(tad_fto_ut_diff + 1))

        log2fc_ripk1_avg = (log2fc_ripk1_1 + log2fc_ripk1_2)/2

        # Calculate combined score with similarity weighting
        combined_score = 0
        if (log2fc_apo_fto > 0 && log2fc_ripk1_1 > 0 && log2fc_ripk1_2 > 0 && log2fc_tad_fto > 0) {
            # Calculate fold change score
            fc_score = log2fc_apo_fto * log2fc_ripk1_avg * log2fc_tad_fto
            
            # Calculate differences from RIPK1 average
            apo_fto_ripk1_diff = abs(log2fc_apo_fto - log2fc_ripk1_avg)
            ripk1_tad_diff = abs(log2fc_ripk1_avg - log2fc_tad_fto)
            
            # Calculate similarity score (inversely proportional to differences)
            sim_score = 1 / (1 + apo_fto_ripk1_diff + ripk1_tad_diff)
            
            # Calculate final score incorporating both fold changes and similarity
            combined_score = fc_score * sim_score
        }

        printf "%s\t%d\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
            cluster, cluster_start[cluster], cluster_end[cluster],
            apo_fto_bt_sum[cluster], apo_fto_ut_sum[cluster],
            apo_ripk1_bt_1_sum[cluster], apo_ripk1_ut_1_sum[cluster],
            apo_ripk1_bt_2_sum[cluster], apo_ripk1_ut_2_sum[cluster],
            tad_fto_bt_sum[cluster], tad_fto_ut_sum[cluster],
            apo_fto_ctrl_bt_sum[cluster], apo_fto_ctrl_ut_sum[cluster],
            apo_ripk1_ctrl_bt_1_sum[cluster], apo_ripk1_ctrl_ut_1_sum[cluster],
            apo_ripk1_ctrl_bt_2_sum[cluster], apo_ripk1_ctrl_ut_2_sum[cluster],
            tad_fto_ctrl_bt_sum[cluster], tad_fto_ctrl_ut_sum[cluster],
            apo_fto_bt_diff, apo_fto_ut_diff,
            apo_ripk1_bt_1_diff, apo_ripk1_ut_1_diff,
            apo_ripk1_bt_2_diff, apo_ripk1_ut_2_diff,
            tad_fto_bt_diff, tad_fto_ut_diff,
            log2fc_apo_fto, log2fc_ripk1_1, log2fc_ripk1_2, log2fc_tad_fto,
            combined_score
    }
}

function log2(x) {
    return log(x)/log(2)
}

function abs(x) {
    return x < 0 ? -x : x
}' combined_edit_rates_with_controls.tsv | \
(head -n 1; tail -n +2 | sort -k32,32nr) > edit_rates_differences.tsv

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
}' "$genes_file" edit_rates_differences.tsv > summed_edit_rates_with_genes_trim_with_controls.tsv

echo -e "\nResults with genes:"
if [ -s summed_edit_rates_with_genes_trim_with_controls.tsv ]; then
    echo "First 5 lines of output:"
    head -n 5 summed_edit_rates_with_genes_trim_with_controls.tsv
    echo -e "\nTotal lines:"
    wc -l summed_edit_rates_with_genes_trim_with_controls.tsv
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
    }
    print $0
}' /Volumes/Lucky_me/9_STAMP_transcript/transcriptome_clean.fa.fai \
   summed_edit_rates_with_genes_trim_with_controls.tsv > summed_edit_rates_with_genes_trim_with_controls_cl.tsv 

R --vanilla --slave << 'EOF'
library(pheatmap)
library(RColorBrewer)

ht <- read.delim("summed_edit_rates_with_genes_trim_with_controls_cl.tsv", 
                 sep="\t", header=TRUE, stringsAsFactors=FALSE)
ht = head(ht, 100)  # Take top 100 rows

# Calculate control sums and order
ht$ctrbtsum <- rowSums(ht[,grep("ctrl_bt",colnames(ht))[c(FALSE,TRUE)]])
ht = ht[order(ht$ctrbtsum),]
write.table(ht, 
            file = "processed_top100_ordered.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)
# Create row labels with chromosome and gene information
row_labels <- ifelse(ht[, ncol(ht)-1] == "NA", 
                    ht$chromosome,
                    paste0(ht$chromosome, "_", ht[, ncol(ht)-1]))

# Prepare matrix data
newht <- ht[,grep("sum",colnames(ht))]
# Prepare matrix data with explicit column order
matrix_data <- as.matrix(ht[,c(
    # BT samples
    "apo_fto_bt_sum", "apo_ripk1_bt_1_sum", "apo_ripk1_bt_2_sum", "tad_fto_bt_sum",
    # UT samples
    "apo_fto_ut_sum", "apo_ripk1_ut_1_sum", "apo_ripk1_ut_2_sum", "tad_fto_ut_sum",
    # BT controls
    "apo_fto_ctrl_bt_sum", "apo_ripk1_ctrl_bt_1_sum", "apo_ripk1_ctrl_bt_2_sum", "tad_fto_ctrl_bt_sum",
    # UT controls
    "apo_fto_ctrl_ut_sum", "apo_ripk1_ctrl_ut_1_sum", "apo_ripk1_ctrl_ut_2_sum", "tad_fto_ctrl_ut_sum"
)])

rownames(matrix_data) <- row_labels

# Clean up column names
colnames(matrix_data) <- gsub("_sum", "", colnames(matrix_data))

# Create heatmap with columns in specified order
png("heatmap_output.png", width=2400, height=3000, res=300)
pheatmap(matrix_data,
         scale="row",
         cluster_cols=FALSE,
         cluster_rows=FALSE,
         show_rownames=TRUE,
         color=colorRampPalette(rev(brewer.pal(n=7, name="RdYlBu")))(100),
         main="Edit Rate Heatmap",
         fontsize_row=7,
         fontsize_col=10)
dev.off()
EOF

#apo-fto and apo-ripk1
awk -v OFS='\t' '
NR == 1 {
    print "chromosome\tstart\tend\t" \
          "apo_fto_bt_sum\tapo_fto_ut_sum\t" \
          "apo_ripk1_bt_1_sum\tapo_ripk1_ut_1_sum\t" \
          "apo_ripk1_bt_2_sum\tapo_ripk1_ut_2_sum\t" \
          "tad_fto_bt_sum\ttad_fto_ut_sum\t" \
          "apo_fto_ctrl_bt_sum\tapo_fto_ctrl_ut_sum\t" \
          "apo_ripk1_ctrl_bt_1_sum\tapo_ripk1_ctrl_ut_1_sum\t" \
          "apo_ripk1_ctrl_bt_2_sum\tapo_ripk1_ctrl_ut_2_sum\t" \
          "tad_fto_ctrl_bt_sum\ttad_fto_ctrl_ut_sum\t" \
          "apo_fto_bt_diff\tapo_fto_ut_diff\t" \
          "apo_ripk1_bt_1_diff\tapo_ripk1_ut_1_diff\t" \
          "apo_ripk1_bt_2_diff\tapo_ripk1_ut_2_diff\t" \
          "tad_fto_bt_diff\ttad_fto_ut_diff\t" \
          "log2fc_apo_fto\tlog2fc_ripk1_1\tlog2fc_ripk1_2\tlog2fc_tad_fto\t" \
          "combined_score"
    next
}

{
    cluster = $1
    clusters[cluster] = 1

    # Track start and end positions
    if (!(cluster in cluster_start)) {
        cluster_start[cluster] = $2
    }
    cluster_end[cluster] = $3

    # Sum values for each condition
    apo_fto_bt_sum[cluster] += $4
    apo_fto_ut_sum[cluster] += $5
    apo_ripk1_bt_1_sum[cluster] += $6
    apo_ripk1_ut_1_sum[cluster] += $7
    apo_ripk1_bt_2_sum[cluster] += $8
    apo_ripk1_ut_2_sum[cluster] += $9
    tad_fto_bt_sum[cluster] += $10
    tad_fto_ut_sum[cluster] += $11
    apo_fto_ctrl_bt_sum[cluster] += $12
    apo_fto_ctrl_ut_sum[cluster] += $13
    apo_ripk1_ctrl_bt_1_sum[cluster] += $14
    apo_ripk1_ctrl_ut_1_sum[cluster] += $15
    apo_ripk1_ctrl_bt_2_sum[cluster] += $16
    apo_ripk1_ctrl_ut_2_sum[cluster] += $17
    tad_fto_ctrl_bt_sum[cluster] += $18
    tad_fto_ctrl_ut_sum[cluster] += $19
}

END {
    for (cluster in clusters) {
        # Calculate differences using summed values
        apo_fto_bt_diff = (apo_fto_bt_sum[cluster] - apo_fto_ctrl_bt_sum[cluster] > 0) ? \
            apo_fto_bt_sum[cluster] - apo_fto_ctrl_bt_sum[cluster] : 0
        apo_fto_ut_diff = (apo_fto_ut_sum[cluster] - apo_fto_ctrl_ut_sum[cluster] > 0) ? \
            apo_fto_ut_sum[cluster] - apo_fto_ctrl_ut_sum[cluster] : 0
        
        apo_ripk1_bt_1_diff = (apo_ripk1_bt_1_sum[cluster] - apo_ripk1_ctrl_bt_1_sum[cluster] > 0) ? \
            apo_ripk1_bt_1_sum[cluster] - apo_ripk1_ctrl_bt_1_sum[cluster] : 0
        apo_ripk1_ut_1_diff = (apo_ripk1_ut_1_sum[cluster] - apo_ripk1_ctrl_ut_1_sum[cluster] > 0) ? \
            apo_ripk1_ut_1_sum[cluster] - apo_ripk1_ctrl_ut_1_sum[cluster] : 0
        
        apo_ripk1_bt_2_diff = (apo_ripk1_bt_2_sum[cluster] - apo_ripk1_ctrl_bt_2_sum[cluster] > 0) ? \
            apo_ripk1_bt_2_sum[cluster] - apo_ripk1_ctrl_bt_2_sum[cluster] : 0
        apo_ripk1_ut_2_diff = (apo_ripk1_ut_2_sum[cluster] - apo_ripk1_ctrl_ut_2_sum[cluster] > 0) ? \
            apo_ripk1_ut_2_sum[cluster] - apo_ripk1_ctrl_ut_2_sum[cluster] : 0
        
        tad_fto_bt_diff = (tad_fto_bt_sum[cluster] - tad_fto_ctrl_bt_sum[cluster] > 0) ? \
            tad_fto_bt_sum[cluster] - tad_fto_ctrl_bt_sum[cluster] : 0
        tad_fto_ut_diff = (tad_fto_ut_sum[cluster] - tad_fto_ctrl_ut_sum[cluster] > 0) ? \
            tad_fto_ut_sum[cluster] - tad_fto_ctrl_ut_sum[cluster] : 0

        # Calculate log2 fold changes using summed values
        log2fc_apo_fto = log2((apo_fto_bt_diff + 1)/(apo_fto_ut_diff + 1))
        log2fc_ripk1_1 = log2((apo_ripk1_bt_1_diff + 1)/(apo_ripk1_ut_1_diff + 1))
        log2fc_ripk1_2 = log2((apo_ripk1_bt_2_diff + 1)/(apo_ripk1_ut_2_diff + 1))
        log2fc_tad_fto = log2((tad_fto_bt_diff + 1)/(tad_fto_ut_diff + 1))

        log2fc_ripk1_avg = (log2fc_ripk1_1 + log2fc_ripk1_2)/2

        # Calculate combined score with similarity weighting
        combined_score = 0
        if (log2fc_apo_fto > 0 && log2fc_ripk1_1 > 0 && log2fc_ripk1_2 > 0) {
            # Calculate fold change score
            fc_score = log2fc_apo_fto * log2fc_ripk1_avg
            
            # Calculate differences from RIPK1 average
            apo_fto_ripk1_diff = abs(log2fc_apo_fto - log2fc_ripk1_avg)
            
            # Calculate similarity score (inversely proportional to differences)
            sim_score = 1 / (1 + apo_fto_ripk1_diff)
            
            # Calculate final score incorporating both fold changes and similarity
            combined_score = fc_score * sim_score
        }

        printf "%s\t%d\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
            cluster, cluster_start[cluster], cluster_end[cluster],
            apo_fto_bt_sum[cluster], apo_fto_ut_sum[cluster],
            apo_ripk1_bt_1_sum[cluster], apo_ripk1_ut_1_sum[cluster],
            apo_ripk1_bt_2_sum[cluster], apo_ripk1_ut_2_sum[cluster],
            tad_fto_bt_sum[cluster], tad_fto_ut_sum[cluster],
            apo_fto_ctrl_bt_sum[cluster], apo_fto_ctrl_ut_sum[cluster],
            apo_ripk1_ctrl_bt_1_sum[cluster], apo_ripk1_ctrl_ut_1_sum[cluster],
            apo_ripk1_ctrl_bt_2_sum[cluster], apo_ripk1_ctrl_ut_2_sum[cluster],
            tad_fto_ctrl_bt_sum[cluster], tad_fto_ctrl_ut_sum[cluster],
            apo_fto_bt_diff, apo_fto_ut_diff,
            apo_ripk1_bt_1_diff, apo_ripk1_ut_1_diff,
            apo_ripk1_bt_2_diff, apo_ripk1_ut_2_diff,
            tad_fto_bt_diff, tad_fto_ut_diff,
            log2fc_apo_fto, log2fc_ripk1_1, log2fc_ripk1_2, log2fc_tad_fto,
            combined_score
    }
}

function log2(x) {
    return log(x)/log(2)
}

function abs(x) {
    return x < 0 ? -x : x
}' combined_edit_rates_with_controls.tsv | \
(head -n 1; tail -n +2 | sort -k32,32nr) > edit_rates_differences_apo.tsv

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
}' "$genes_file" edit_rates_differences_apo.tsv > summed_edit_rates_with_genes_trim_with_controls_apo.tsv

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
    }
    print $0
}' /Volumes/Lucky_me/9_STAMP_transcript/transcriptome_clean.fa.fai \
   summed_edit_rates_with_genes_trim_with_controls_apo.tsv > summed_edit_rates_with_genes_trim_with_controls_apo_cl.tsv 

R --vanilla --slave << 'EOF'
library(pheatmap)
library(RColorBrewer)

ht <- read.delim("summed_edit_rates_with_genes_trim_with_controls_apo_cl.tsv", 
                 sep="\t", header=TRUE, stringsAsFactors=FALSE)
ht = head(ht, 100)  # Take top 100 rows

# Calculate control sums and order
ht$ctrbtsum <- rowSums(ht[,grep("ctrl_bt",colnames(ht))[c(FALSE,TRUE)]])
ht = ht[order(ht$ctrbtsum),]
write.table(ht, 
            file = "processed_top100_ordered_apo.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)
# Create row labels with chromosome and gene information
row_labels <- ifelse(ht[, ncol(ht)-1] == "NA", 
                    ht$chromosome,
                    paste0(ht$chromosome, "_", ht[, ncol(ht)-1]))

# Prepare matrix data
newht <- ht[,grep("sum",colnames(ht))]
# Prepare matrix data with explicit column order
matrix_data <- as.matrix(ht[,c(
    # BT samples
    "apo_fto_bt_sum", "apo_ripk1_bt_1_sum", "apo_ripk1_bt_2_sum", "tad_fto_bt_sum",
    # UT samples
    "apo_fto_ut_sum", "apo_ripk1_ut_1_sum", "apo_ripk1_ut_2_sum", "tad_fto_ut_sum",
    # BT controls
    "apo_fto_ctrl_bt_sum", "apo_ripk1_ctrl_bt_1_sum", "apo_ripk1_ctrl_bt_2_sum", "tad_fto_ctrl_bt_sum",
    # UT controls
    "apo_fto_ctrl_ut_sum", "apo_ripk1_ctrl_ut_1_sum", "apo_ripk1_ctrl_ut_2_sum", "tad_fto_ctrl_ut_sum"
)])

rownames(matrix_data) <- row_labels

# Clean up column names
colnames(matrix_data) <- gsub("_sum", "", colnames(matrix_data))

# Create heatmap with columns in specified order
png("heatmap_output_apo.png", width=2400, height=3000, res=300)
pheatmap(matrix_data,
         scale="row",
         cluster_cols=FALSE,
         cluster_rows=FALSE,
         show_rownames=TRUE,
         color=colorRampPalette(rev(brewer.pal(n=7, name="RdYlBu")))(100),
         main="Edit Rate Heatmap",
         fontsize_row=7,
         fontsize_col=10)
dev.off()
EOF

#tad-fto and apo-ripk1
awk -v OFS='\t' '
NR == 1 {
    print "chromosome\tstart\tend\t" \
          "apo_fto_bt_sum\tapo_fto_ut_sum\t" \
          "apo_ripk1_bt_1_sum\tapo_ripk1_ut_1_sum\t" \
          "apo_ripk1_bt_2_sum\tapo_ripk1_ut_2_sum\t" \
          "tad_fto_bt_sum\ttad_fto_ut_sum\t" \
          "apo_fto_ctrl_bt_sum\tapo_fto_ctrl_ut_sum\t" \
          "apo_ripk1_ctrl_bt_1_sum\tapo_ripk1_ctrl_ut_1_sum\t" \
          "apo_ripk1_ctrl_bt_2_sum\tapo_ripk1_ctrl_ut_2_sum\t" \
          "tad_fto_ctrl_bt_sum\ttad_fto_ctrl_ut_sum\t" \
          "apo_fto_bt_diff\tapo_fto_ut_diff\t" \
          "apo_ripk1_bt_1_diff\tapo_ripk1_ut_1_diff\t" \
          "apo_ripk1_bt_2_diff\tapo_ripk1_ut_2_diff\t" \
          "tad_fto_bt_diff\ttad_fto_ut_diff\t" \
          "log2fc_apo_fto\tlog2fc_ripk1_1\tlog2fc_ripk1_2\tlog2fc_tad_fto\t" \
          "combined_score"
    next
}

{
    cluster = $1
    clusters[cluster] = 1

    # Track start and end positions
    if (!(cluster in cluster_start)) {
        cluster_start[cluster] = $2
    }
    cluster_end[cluster] = $3

    # Sum values for each condition
    apo_fto_bt_sum[cluster] += $4
    apo_fto_ut_sum[cluster] += $5
    apo_ripk1_bt_1_sum[cluster] += $6
    apo_ripk1_ut_1_sum[cluster] += $7
    apo_ripk1_bt_2_sum[cluster] += $8
    apo_ripk1_ut_2_sum[cluster] += $9
    tad_fto_bt_sum[cluster] += $10
    tad_fto_ut_sum[cluster] += $11
    apo_fto_ctrl_bt_sum[cluster] += $12
    apo_fto_ctrl_ut_sum[cluster] += $13
    apo_ripk1_ctrl_bt_1_sum[cluster] += $14
    apo_ripk1_ctrl_ut_1_sum[cluster] += $15
    apo_ripk1_ctrl_bt_2_sum[cluster] += $16
    apo_ripk1_ctrl_ut_2_sum[cluster] += $17
    tad_fto_ctrl_bt_sum[cluster] += $18
    tad_fto_ctrl_ut_sum[cluster] += $19
}

END {
    for (cluster in clusters) {
        # Calculate differences using summed values
        apo_fto_bt_diff = (apo_fto_bt_sum[cluster] - apo_fto_ctrl_bt_sum[cluster] > 0) ? \
            apo_fto_bt_sum[cluster] - apo_fto_ctrl_bt_sum[cluster] : 0
        apo_fto_ut_diff = (apo_fto_ut_sum[cluster] - apo_fto_ctrl_ut_sum[cluster] > 0) ? \
            apo_fto_ut_sum[cluster] - apo_fto_ctrl_ut_sum[cluster] : 0
        
        apo_ripk1_bt_1_diff = (apo_ripk1_bt_1_sum[cluster] - apo_ripk1_ctrl_bt_1_sum[cluster] > 0) ? \
            apo_ripk1_bt_1_sum[cluster] - apo_ripk1_ctrl_bt_1_sum[cluster] : 0
        apo_ripk1_ut_1_diff = (apo_ripk1_ut_1_sum[cluster] - apo_ripk1_ctrl_ut_1_sum[cluster] > 0) ? \
            apo_ripk1_ut_1_sum[cluster] - apo_ripk1_ctrl_ut_1_sum[cluster] : 0
        
        apo_ripk1_bt_2_diff = (apo_ripk1_bt_2_sum[cluster] - apo_ripk1_ctrl_bt_2_sum[cluster] > 0) ? \
            apo_ripk1_bt_2_sum[cluster] - apo_ripk1_ctrl_bt_2_sum[cluster] : 0
        apo_ripk1_ut_2_diff = (apo_ripk1_ut_2_sum[cluster] - apo_ripk1_ctrl_ut_2_sum[cluster] > 0) ? \
            apo_ripk1_ut_2_sum[cluster] - apo_ripk1_ctrl_ut_2_sum[cluster] : 0
        
        tad_fto_bt_diff = (tad_fto_bt_sum[cluster] - tad_fto_ctrl_bt_sum[cluster] > 0) ? \
            tad_fto_bt_sum[cluster] - tad_fto_ctrl_bt_sum[cluster] : 0
        tad_fto_ut_diff = (tad_fto_ut_sum[cluster] - tad_fto_ctrl_ut_sum[cluster] > 0) ? \
            tad_fto_ut_sum[cluster] - tad_fto_ctrl_ut_sum[cluster] : 0

        # Calculate log2 fold changes using summed values
        log2fc_apo_fto = log2((apo_fto_bt_diff + 1)/(apo_fto_ut_diff + 1))
        log2fc_ripk1_1 = log2((apo_ripk1_bt_1_diff + 1)/(apo_ripk1_ut_1_diff + 1))
        log2fc_ripk1_2 = log2((apo_ripk1_bt_2_diff + 1)/(apo_ripk1_ut_2_diff + 1))
        log2fc_tad_fto = log2((tad_fto_bt_diff + 1)/(tad_fto_ut_diff + 1))

        log2fc_ripk1_avg = (log2fc_ripk1_1 + log2fc_ripk1_2)/2

        # Calculate combined score with similarity weighting
        combined_score = 0
        if (log2fc_ripk1_1 > 0 && log2fc_ripk1_2 > 0 && log2fc_tad_fto > 0) {
            # Calculate fold change score
            fc_score = log2fc_tad_fto * log2fc_ripk1_avg
            
            # Calculate differences from RIPK1 average
            tad_fto_ripk1_diff = abs(log2fc_tad_fto - log2fc_ripk1_avg)
            
            # Calculate similarity score (inversely proportional to differences)
            sim_score = 1 / (1 + tad_fto_ripk1_diff)
            
            # Calculate final score incorporating both fold changes and similarity
            combined_score = fc_score * sim_score
        }

        printf "%s\t%d\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
            cluster, cluster_start[cluster], cluster_end[cluster],
            apo_fto_bt_sum[cluster], apo_fto_ut_sum[cluster],
            apo_ripk1_bt_1_sum[cluster], apo_ripk1_ut_1_sum[cluster],
            apo_ripk1_bt_2_sum[cluster], apo_ripk1_ut_2_sum[cluster],
            tad_fto_bt_sum[cluster], tad_fto_ut_sum[cluster],
            apo_fto_ctrl_bt_sum[cluster], apo_fto_ctrl_ut_sum[cluster],
            apo_ripk1_ctrl_bt_1_sum[cluster], apo_ripk1_ctrl_ut_1_sum[cluster],
            apo_ripk1_ctrl_bt_2_sum[cluster], apo_ripk1_ctrl_ut_2_sum[cluster],
            tad_fto_ctrl_bt_sum[cluster], tad_fto_ctrl_ut_sum[cluster],
            apo_fto_bt_diff, apo_fto_ut_diff,
            apo_ripk1_bt_1_diff, apo_ripk1_ut_1_diff,
            apo_ripk1_bt_2_diff, apo_ripk1_ut_2_diff,
            tad_fto_bt_diff, tad_fto_ut_diff,
            log2fc_apo_fto, log2fc_ripk1_1, log2fc_ripk1_2, log2fc_tad_fto,
            combined_score
    }
}

function log2(x) {
    return log(x)/log(2)
}

function abs(x) {
    return x < 0 ? -x : x
}' combined_edit_rates_with_controls.tsv | \
(head -n 1; tail -n +2 | sort -k32,32nr) > edit_rates_differences_tad.tsv

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
}' "$genes_file" edit_rates_differences_tad.tsv > summed_edit_rates_with_genes_trim_with_controls_tad.tsv

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
    }
    print $0
}' /Volumes/Lucky_me/9_STAMP_transcript/transcriptome_clean.fa.fai \
   summed_edit_rates_with_genes_trim_with_controls_tad.tsv > summed_edit_rates_with_genes_trim_with_controls_tad_cl.tsv 

R --vanilla --slave << 'EOF'
library(pheatmap)
library(RColorBrewer)

ht <- read.delim("summed_edit_rates_with_genes_trim_with_controls_tad_cl.tsv", 
                 sep="\t", header=TRUE, stringsAsFactors=FALSE)
ht = head(ht, 100)  # Take top 100 rows

# Calculate control sums and order
ht$ctrbtsum <- rowSums(ht[,grep("ctrl_bt",colnames(ht))[c(FALSE,TRUE)]])
ht = ht[order(ht$ctrbtsum),]
write.table(ht, 
            file = "processed_top100_ordered_tad.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)
# Create row labels with chromosome and gene information
row_labels <- ifelse(ht[, ncol(ht)-1] == "NA", 
                    ht$chromosome,
                    paste0(ht$chromosome, "_", ht[, ncol(ht)-1]))

# Prepare matrix data
newht <- ht[,grep("sum",colnames(ht))]
# Prepare matrix data with explicit column order
matrix_data <- as.matrix(ht[,c(
    # BT samples
    "apo_fto_bt_sum", "apo_ripk1_bt_1_sum", "apo_ripk1_bt_2_sum", "tad_fto_bt_sum",
    # UT samples
    "apo_fto_ut_sum", "apo_ripk1_ut_1_sum", "apo_ripk1_ut_2_sum", "tad_fto_ut_sum",
    # BT controls
    "apo_fto_ctrl_bt_sum", "apo_ripk1_ctrl_bt_1_sum", "apo_ripk1_ctrl_bt_2_sum", "tad_fto_ctrl_bt_sum",
    # UT controls
    "apo_fto_ctrl_ut_sum", "apo_ripk1_ctrl_ut_1_sum", "apo_ripk1_ctrl_ut_2_sum", "tad_fto_ctrl_ut_sum"
)])

rownames(matrix_data) <- row_labels

# Clean up column names
colnames(matrix_data) <- gsub("_sum", "", colnames(matrix_data))

# Create heatmap with columns in specified order
png("heatmap_output_tad.png", width=2400, height=3000, res=300)
pheatmap(matrix_data,
         scale="row",
         cluster_cols=FALSE,
         cluster_rows=FALSE,
         show_rownames=TRUE,
         color=colorRampPalette(rev(brewer.pal(n=7, name="RdYlBu")))(100),
         main="Edit Rate Heatmap",
         fontsize_row=7,
         fontsize_col=10)
dev.off()
EOF
