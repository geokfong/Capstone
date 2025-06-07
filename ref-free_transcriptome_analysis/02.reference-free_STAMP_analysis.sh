samples=(
"SKOV3_Apo-Cnt-FTO_BT1_ID3_R39"
"SKOV3_Apo-Cnt-FTO_BT1_ID4_R40"
"SKOV3_Apo-Cnt-FTO_UT_ID1_R37"
"SKOV3_Apo-Cnt-FTO_UT_ID2_R38"
"SKOV3_Apo-Cnt-RIPK1_BT1_ID10_R10"
"SKOV3_Apo-Cnt-RIPK1_BT1_ID10_R9_2"
"SKOV3_Apo-Cnt-RIPK1_BT1_ID11_R10_2"
"SKOV3_Apo-Cnt-RIPK1_BT1_ID9_R9"
"SKOV3_Apo-Cnt-RIPK1_UT_ID1_R1"
"SKOV3_Apo-Cnt-RIPK1_UT_ID2_R1_2"
"SKOV3_Apo-Cnt-RIPK1_UT_ID2_R2"
"SKOV3_Apo-Cnt-RIPK1_UT_ID3_R2_2"
"SKOV3_Apo-FTO_BT1_ID10_R45"
"SKOV3_Apo-FTO_BT1_ID11_R46"
"SKOV3_Apo-FTO_UT_ID7_R43"
"SKOV3_Apo-FTO_UT_ID9_R44"
"SKOV3_Apo-RIPK1_BT1_ID11_R11"
"SKOV3_Apo-RIPK1_BT1_ID12_R11_2"
"SKOV3_Apo-RIPK1_BT1_ID12_R12"
"SKOV3_Apo-RIPK1_BT1_ID13_R12_2"
"SKOV3_Apo-RIPK1_UT_ID5_R5"
"SKOV3_Apo-RIPK1_UT_ID6_R5_2"
"SKOV3_Apo-RIPK1_UT_ID6_R6"
"SKOV3_Apo-RIPK1_UT_ID7_R6_2"
"SKOV3_Tad-Cnt-FTO_BT1_ID3_R57"
"SKOV3_Tad-Cnt-FTO_BT1_ID4_R58"
"SKOV3_Tad-Cnt-FTO_UT_ID1_R55"
"SKOV3_Tad-Cnt-FTO_UT_ID2_R56"
"SKOV3_Tad-FTO_BT1_ID10_R64"
"SKOV3_Tad-FTO_BT1_ID9_R63"
"SKOV3_Tad-FTO_UT_ID7_R61"
"SKOV3_Tad-FTO_UT_ID8_R62"
)

genome_fa="/home/onggf/data/7_STAMP/transcriptome/output_4/transcriptome_clean.fa"
ncpus=48
cov=1
strandness="R"
parent_folder="/home/onggf/data/7_STAMP/data/3Mar_data/01.final_fastq"

hisat3n_index_name_tad="Transcript_AG"
hisat3n_index_tad="/home/onggf/data/7_STAMP/transcriptome/output_4/$hisat3n_index_name_tad/$hisat3n_index_name_tad"
hisat3n_out_dir_tad="/home/onggf/data/7_STAMP/data/3Mar_data/SKOV3/03.transcript_analysis"

hisat3n_index_name_apo="Transcript_CT"
hisat3n_index_apo="/home/onggf/data/7_STAMP/transcriptome/output_4/$hisat3n_index_name_apo/$hisat3n_index_name_apo"
hisat3n_out_dir_apo="/home/onggf/data/7_STAMP/data/3Mar_data/SKOV3/03.transcript_analysis"

process_sample() {
  local sample_type=$1
  local hisat3n_index=$2
  local hisat3n_out_dir=$3
  local log_file="$hisat3n_out_dir/$sample_type.processing.log"

  for sample_name in "${samples[@]}"; do
    echo "Processing $sample_name for $sample_type"
    
    sample_dir="$parent_folder/$sample_name"
    
    sample_r1=$(find "$sample_dir" -type f -name "*1.fq.gz" | head -n 1)
    sample_r2=$(find "$sample_dir" -type f -name "*2.fq.gz" | head -n 1)

    if [ ! -f "$sample_r1" ] || [ ! -f "$sample_r2" ]; then
      echo "Error: Paired-end files for $sample_name not found!" | tee -a "$log_file"
      continue
    fi

    if [[ "$sample_type" == "Tad" ]]; then
      base_change="A,G"
      hisat3n_index_name=$hisat3n_index_name_tad
    elif [[ "$sample_type" == "Apo" ]]; then
      base_change="C,T"
      hisat3n_index_name=$hisat3n_index_name_apo
    else
      echo "Unknown sample type: $sample_type" | tee -a "$log_file"
      continue
    fi

    mkdir -p "$hisat3n_out_dir/$sample_name"
    read_len_r1=$(zcat "$sample_r1" | head -2 | tail -1 | wc -c)
    read_len_r2=$(zcat "$sample_r2" | head -2 | tail -1 | wc -c)
    crop_r1=$(expr $read_len_r1 - 16)
    crop_r2=$(expr $read_len_r2 - 16)

    echo "Running Trimmomatic for $sample_name"
    trimmomatic PE -threads $ncpus -phred33 \
    "$sample_r1" "$sample_r2" \
    "$hisat3n_out_dir/$sample_name/${sample_name}_trimmed_R1.fastq.gz" \
    "$hisat3n_out_dir/$sample_name/${sample_name}_unpaired_R1.fastq.gz" \
    "$hisat3n_out_dir/$sample_name/${sample_name}_trimmed_R2.fastq.gz" \
    "$hisat3n_out_dir/$sample_name/${sample_name}_unpaired_R2.fastq.gz" \
    HEADCROP:15 CROP:$crop_r1 \
    HEADCROP:15 CROP:$crop_r2 \
    2>&1
    if [ $? -ne 0 ]; then
      echo "Error during Trimmomatic for $sample_name" | tee -a "$log_file"
      continue
    fi

    echo "Running HISAT-3n for $sample_name"
    /home/onggf/software/hisat/hisat-3n/hisat-3n -p $ncpus --time --base-change "$base_change" \
      --repeat --repeat-limit 1000 --bowtie2-dp 0 --no-unal --rna-strandness $strandness \
      -x "$hisat3n_index" -U "$hisat3n_out_dir/$sample_name/${sample_name}_trimmed_R2.fastq.gz" | \
      samtools view -@ $ncpus -Shb - -o "$hisat3n_out_dir/$sample_name/${sample_name}.$hisat3n_index_name.align.raw.bam" 2>&1
    if [ $? -ne 0 ]; then
      echo "Error during HISAT-3n for $sample_name" | tee -a "$log_file"
      continue
    fi

    echo "Sorting and indexing BAM for $sample_name"
    samtools view -@ $ncpus -q 1 -hb "$hisat3n_out_dir/$sample_name/${sample_name}.$hisat3n_index_name.align.raw.bam" | \
    samtools sort -T $hisat3n_out_dir/$sample_name -@ $ncpus -o "$hisat3n_out_dir/$sample_name/${sample_name}.$hisat3n_index_name.align.sorted.bam" -
    samtools index "$hisat3n_out_dir/$sample_name/${sample_name}.$hisat3n_index_name.align.sorted.bam" 2>&1
    if [ $? -ne 0 ]; then
      echo "Error during BAM indexing for $sample_name" | tee -a "$log_file"
      continue
    fi

    echo "Running pileup2var for $sample_name"
    /home/onggf/software/pileup2var/pileup2var -t $ncpus -f 524 -c $cov -s $strandness \
      -g "$genome_fa" -b "$hisat3n_out_dir/$sample_name/${sample_name}.$hisat3n_index_name.align.sorted.bam" \
      -o "$hisat3n_out_dir/$sample_name/${sample_name}.$hisat3n_index_name.pileup2var.flt.txt" 2>&1
    if [ $? -ne 0 ]; then
      echo "Error during pileup2var for $sample_name" | tee -a "$log_file"
      continue
    fi

    echo "Completed processing $sample_name for $sample_type"
  done
}

echo "Starting analysis..."

# For Tad samples
process_sample "Tad" "$hisat3n_index_tad" "$hisat3n_out_dir_tad"

# For Apo samples
process_sample "Apo" "$hisat3n_index_apo" "$hisat3n_out_dir_apo"

# Loop over all directories containing the files
for dir in /home/onggf/data/7_STAMP/data/3Mar_data/SKOV3/03.transcript_analysis/*; do
    if [ -d "$dir" ]; then  
        pile_file=$(find "$dir" -type f -name "*.Transcript_CT.pileup2var.flt.txt" -o -name "*.Transcript_AG.pileup2var.flt.txt" | head -n 1)
        
        # Proceed only if a pileup file exists
        if [ -f "$pile_file" ]; then
            # Check if the file is smaller than 4KB
            if [ $(stat -c %s "$pile_file") -lt 4096 ]; then
                echo "$pile_file is less than 4KB, concatenating cluster files..."
                
                # Find and concatenate *cluster* files in the directory
                find "$dir" -type f -name "*cluster*" -print0 | xargs -0 cat >> "$pile_file"
                echo "Appended cluster files to $pile_file"
            else
                echo "$pile_file is not smaller than 4KB, skipping concatenation."
            fi
        else
            echo "No pileup file found in $dir, skipping."
        fi
    fi
done
