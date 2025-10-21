# Create a file with the list of filenames
cat <<EOL > file_list.txt
GTEX-1C475-0726-SM-73KVL.Esophagus_Muscularis.RNAseq.bw
GTEX-1C475-1826-SM-73KWA.Skin_Sun_Exposed_Lower_leg.RNAseq.bw
GTEX-1GMR3-0626-SM-9WYT3.Artery_Coronary.RNAseq.bw
GTEX-1H1E6-0826-SM-9WG83.Pancreas.RNAseq.bw
GTEX-1HFI6-0011-R7b-SM-CM2SS.Brain_Putamen_basal_ganglia.RNAseq.bw
GTEX-1HGF4-0011-R5b-SM-CM2ST.Brain_Caudate_basal_ganglia.RNAseq.bw
GTEX-1HSGN-0726-SM-A9G2F.Thyroid.RNAseq.bw
GTEX-1HSKV-0011-R1b-SM-CMKH7.Brain_Hippocampus.RNAseq.bw
GTEX-1I1GU-1226-SM-A9SKT.Esophagus_Gastroesophageal_Junction.RNAseq.bw
GTEX-1IDJC-1326-SM-CL53H.Colon_Transverse.RNAseq.bw
GTEX-1IDJU-1026-SM-AHZ2U.Vagina.RNAseq.bw
GTEX-1JKYN-1026-SM-CGQG4.Testis.RNAseq.bw
GTEX-1JN76-0626-SM-CKZOQ.Skin_Not_Sun_Exposed_Suprapubic.RNAseq.bw
GTEX-1KXAM-1926-SM-D3LAG.Colon_Sigmoid.RNAseq.bw
GTEX-1LG7Z-0005-SM-DKPQ6.Whole_Blood.RNAseq.bw
GTEX-1MA7W-1526-SM-DHXKS.Uterus.RNAseq.bw
GTEX-1PIEJ-1526-SM-E6CP8.Small_Intestine_Terminal_Ileum.RNAseq.bw
GTEX-11NSD-1126-SM-5N9BQ.Esophagus_Mucosa.RNAseq.bw
GTEX-13OVI-1126-SM-5KLZF.Kidney_Cortex.RNAseq.bw
GTEX-13S86-0326-SM-5SI6K.Heart_Atrial_Appendage.RNAseq.bw
GTEX-13X6J-0011-R11a-SM-5P9HE.Brain_Cerebellar_Hemisphere.RNAseq.bw
GTEX-14BIN-0011-R6a-SM-5S2RH.Brain_Nucleus_accumbens_basal_ganglia.RNAseq.bw
GTEX-14BMU-0626-SM-73KZ6.Adipose_Visceral_Omentum.RNAseq.bw
GTEX-14DAR-1026-SM-73KV3.Prostate.RNAseq.bw
GTEX-14PKU-0526-SM-6871A.Spleen.RNAseq.bw
GTEX-14PN4-0011-R3b-SM-686ZU.Brain_Anterior_cingulate_cortex_BA24.RNAseq.bw
GTEX-117XS-0008-SM-5Q5DQ.Cells_Cultured_fibroblasts.RNAseq.bw
GTEX-145MH-2926-SM-5Q5D2.Brain_Cerebellum.RNAseq.bw
GTEX-1122O-0003-SM-5Q5DL.Cells_EBV-transformed_lymphocytes.RNAseq.bw
GTEX-NFK9-0326-SM-3MJGV.Adipose_Subcutaneous.RNAseq.bw
GTEX-NFK9-0626-SM-2HMIV.Muscle_Skeletal.RNAseq.bw
GTEX-NFK9-0926-SM-2HMJU.Heart_Left_Ventricle.RNAseq.bw
GTEX-NFK9-1526-SM-3LK7B.Stomach.RNAseq.bw
GTEX-OHPK-2326-SM-3MJH2.Fallopian_Tube.RNAseq.bw
GTEX-S3XE-1226-SM-4AD4L.Bladder.RNAseq.bw
GTEX-S341-1126-SM-4AD6T.Cervix_Ectocervix.RNAseq.bw
GTEX-T5JC-0011-R4A-SM-32PLT.Brain_Amygdala.RNAseq.bw
GTEX-T5JC-0011-R8A-SM-32PLM.Brain_Hypothalamus.RNAseq.bw
GTEX-T5JC-0011-R10A-SM-32PM2.Brain_Frontal_Cortex_BA9.RNAseq.bw
GTEX-T5JC-1626-SM-EZ6KW.Kidney_Medulla.RNAseq.bw
GTEX-TML8-1626-SM-32QOO.Nerve_Tibial.RNAseq.bw
GTEX-UTHO-3026-SM-3GAFB.Brain_Cortex.RNAseq.bw
GTEX-WYVS-0426-SM-4ONDL.Artery_Aorta.RNAseq.bw
GTEX-XPT6-2226-SM-4B66R.Artery_Tibial.RNAseq.bw
GTEX-Y5LM-0126-SM-4VBRL.Adrenal_Gland.RNAseq.bw
GTEX-Y5LM-0426-SM-4VBRO.Liver.RNAseq.bw
GTEX-Y5LM-1826-SM-4VDT9.Minor_Salivary_Gland.RNAseq.bw
GTEX-Y5V5-0826-SM-4VBQD.Lung.RNAseq.bw
GTEX-Y111-2926-SM-4TT25.Pituitary.RNAseq.bw
GTEX-YFC4-0011-R9a-SM-4SOK4.Brain_Spinal_cord_cervical_c-1.RNAseq.bw
GTEX-Z93S-0011-R2a-SM-4RGNG.Brain_Substantia_nigra.RNAseq.bw
GTEX-ZPIC-1326-SM-DO91Y.Cervix_Endocervix.RNAseq.bw
GTEX-ZT9W-2026-SM-51MRA.Breast_Mammary_Tissue.RNAseq.bw
GTEX-ZVT2-0326-SM-5E44G.Ovary.RNAseq.bw
EOL
mkdir -p tissue_expression

while read -r file; do
  wget -P tissue_expression "https://hgdownload.soe.ucsc.edu/gbdb/hg38/gtex/cov/$file"
done < file_list.txt

cd ~/Downloads
curl -O http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.arm64/bigWigToBedGraph
chmod +x bigWigToBedGraph
sudo mv bigWigToBedGraph /usr/local/bin/

which bigWigToBedGraph
bigWigToBedGraph -h

mkdir -p bedgraph
for bw in *.bw; do
    base=$(basename "$bw" .bw)
    bigWigToBedGraph "$bw" "bedgraph/${base}.bedgraph"
done

for bedgraph in bedgraph/*.bedgraph; do
    tissue=$(basename "$bedgraph" | cut -d'.' -f2)

    HEADER="chrom\tstart\tend\tfeature_info\texpression_count\texpression_mean\texpression_sum"
    echo -e "$HEADER" > "/Volumes/Lucky me/8_STAMP/bin/04.sample_metrics_final_v3/expression_level/tissue_expression/output/${tissue}_output.tsv"

    bedtools map \
        -a "/Volumes/Lucky me/8_STAMP/bin/04.sample_metrics_final_v3/expression_level/region_data.bed" \
        -b "$bedgraph" \
        -c 4 \
        -o count,mean,sum >> "/Volumes/Lucky me/8_STAMP/bin/04.sample_metrics_final_v3/expression_level/tissue_expression/output/${tissue}_output.tsv"
done

output_dir="/Volumes/Lucky me/8_STAMP/bin/04.sample_metrics_final_v3/expression_level/tissue_expression/output"
merged_file="${output_dir}/combined_expression_matrix.tsv"

first_file=$(ls "${output_dir}"/*_output.tsv | head -n 1)
cut -f1-4 "$first_file" | tail -n +2 > "$merged_file"
headers=("chrom" "start" "end" "feature_info")

for file in "${output_dir}"/*_output.tsv; do
    tissue=$(basename "$file" | sed 's/_output\.tsv//')

    # Skip the header line and cut metrics columns (5–7), save them to a temporary file
    temp_file=$(mktemp)
    cut -f5-7 "$file" | tail -n +2 > "$temp_file"

    # Append tissue-specific columns to headers
    headers+=("${tissue}_count" "${tissue}_mean" "${tissue}_sum")

    paste "$merged_file" "$temp_file" > "${merged_file}.tmp"
    mv "${merged_file}.tmp" "$merged_file"
done

header_line=$(IFS=$'\t'; echo "${headers[*]}")
echo -e "$header_line" | cat - "$merged_file" > "${merged_file}.tmp" && mv "${merged_file}.tmp" "$merged_file"
