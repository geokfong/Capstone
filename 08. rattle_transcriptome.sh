# Step 1: Adapter trimming with Porechop
porechop -i /home/onggf/data/7_STAMP/transcriptome/ONT_BT/fastq/concatenated_SKOV3_BT6_MDA_BT2.fastq \
         -o /home/onggf/data/7_STAMP/transcriptome/output_porechop/concatenated_SKOV3_BT6_MDA_BT2_trimmed.fastq \
         -t 48  # Use 48 threads for faster processing

# Step 2: Environment setup for RATTLE
mamba deactivate
mamba activate RATTLE  # Activate the RATTLE environment

INPUT="/home/onggf/data/7_STAMP/transcriptome/output_porechop/concatenated_SKOV3_BT6_MDA_BT2_trimmed.fastq"
BASE_OUTPUT="/home/onggf/data/7_STAMP/transcriptome/output_4"
RATTLE="/home/onggf/software/RATTLE/rattle"

mkdir -p "${BASE_OUTPUT}/split_fastq"
mkdir -p "${BASE_OUTPUT}/split_output"

# Step 3: Split input FASTQ into 4 parts for parallel processing
echo "Splitting input file into 4 parts..."
seqkit seq -m 150 "${INPUT}" | \  # Keep reads longer than 150bp
seqkit split2 --by-part 4 -O "${BASE_OUTPUT}/split_fastq" -f  # Split into 4 equal parts

if ! ls "${BASE_OUTPUT}/split_fastq"/*.fastq 1> /dev/null 2>&1; then
    echo "Error: No split FASTQ files found. Exiting."
    exit 1
fi

# Step 4: Process each split FASTQ file with RATTLE
for split_file in "${BASE_OUTPUT}/split_fastq"/*.fastq; do
    base_name=$(basename "${split_file}" .fastq)
    output_dir="${BASE_OUTPUT}/split_output/${base_name}"
    
    echo "Processing ${base_name}..."
    mkdir -p "${output_dir}/clusters"

    # RATTLE clustering step
    "${RATTLE}" cluster \
        -i "${split_file}" \
        -t 32 \
        -o "${output_dir}" \
        --rna

    # Check if clustering was successful
    if [ ! -f "${output_dir}/clusters.out" ]; then
        echo "Error: Clustering failed for ${base_name}. Skipping..."
        continue
    fi

    # Generate cluster summary
    "${RATTLE}" cluster_summary \
        -i "${split_file}" \
        -c "${output_dir}/clusters.out" \
        > "${output_dir}/cluster_summary.tsv"

    # Extract clusters into individual files
    "${RATTLE}" extract_clusters \
        -i "${split_file}" \
        -c "${output_dir}/clusters.out" \
        -o "${output_dir}/clusters" \
        --fastq

    # Correct reads within clusters
    "${RATTLE}" correct \
        -i "${split_file}" \
        -c "${output_dir}/clusters.out" \
        -o "${output_dir}" \
        -t 32

    echo "Completed processing ${base_name}"
done

echo "All processing complete!"

# Step 5: Merge consensus sequences if available
if ls ${BASE_OUTPUT}/split_output/*/consensi.fq 1> /dev/null 2>&1; then
    echo "Merging consensus FASTQ files..."
    cat ${BASE_OUTPUT}/split_output/*/consensi.fq > ${BASE_OUTPUT}/combined_consensi.fq
else
    echo "Error: No consensi.fq files found. Skipping polishing step."
    exit 1
fi

# Step 6: RATTLE polishing to generate transcriptome
/home/onggf/software/RATTLE/rattle polish \
    -i ${BASE_OUTPUT}/combined_consensi.fq \
    -o ${BASE_OUTPUT} \
    -t 32 \
    --rna

# Step 7: Convert polished transcriptome from FASTQ to FASTA
seqtk seq -a /home/onggf/data/7_STAMP/transcriptome/output_4/transcriptome.fq > /home/onggf/data/7_STAMP/transcriptome/output_4/transcriptome.fa

# Step 8: Replace uracil (U) with thymine (T) for DNA-based indexing
sed 's/U/T/g' /home/onggf/data/7_STAMP/transcriptome/output_4/transcriptome.fa > /home/onggf/data/7_STAMP/transcriptome/output_4/transcriptome_clean.fa

# Step 9: Build HISAT-3N indices with specific base conversions
# AG conversion: for detecting A-to-G RNA (for TadA) editing events
/home/onggf/software/hisat/hisat-3n/hisat-3n-build --base-change A,G \
    /home/onggf/data/7_STAMP/transcriptome/output_4/transcriptome_clean.fa \
    /home/onggf/data/7_STAMP/transcriptome/output_4/Transcript_AG/Transcript_AG

# CT conversion: for detecting C-to-T (for APOBEC) editing events
/home/onggf/software/hisat/hisat-3n/hisat-3n-build --base-change C,T \
    /home/onggf/data/7_STAMP/transcriptome/output_4/transcriptome_clean.fa \
    /home/onggf/data/7_STAMP/transcriptome/output_4/Transcript_CT/Transcript_CT
