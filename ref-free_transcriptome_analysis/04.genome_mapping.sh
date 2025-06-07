minimap2 -ax splice -k14 \
    /home/onggf/data/hg38_etam/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
    /home/onggf/data/7_STAMP/transcriptome/selected_clusters.fa \
    > minimap2_results/mapping.sam

samtools view -bS minimap2_results/mapping.sam | samtools sort > minimap2_results/mapping_sorted.bam
samtools index minimap2_results/mapping_sorted.bam
samtools view -b -F 256 -F 2048 minimap2_results/mapping_sorted.bam > minimap2_results/mapping_primary_only.bam
samtools sort minimap2_results/mapping_primary_only.bam > minimap2_results/mapping_primary_sorted.bam
samtools index minimap2_results/mapping_primary_sorted.bam
bedtools bamtobed -splitD -i minimap2_results/mapping_primary_sorted.bam > minimap2_results/transcript_locations_primary.bed

mkdir -p bt_sample_edits

find /Volumes/Lucky_me/9_STAMP_transcript/07.transcript_analysis_v3 -name "*_edits_trim.bed" | grep -i "bt" | grep -v -i "ut" > bt_files.txt
echo "#chr	start	end	name	score	strand" > bt_sample_edits/all_cluster_edits.bed

while IFS= read -r cluster; do
    echo "Searching for $cluster..."
    
    cluster_id=$(echo "$cluster" | cut -d'_' -f1-2)
    
    while IFS= read -r bed_file; do
        if [ -f "$bed_file" ]; then
            matches=$(grep "$cluster_id" "$bed_file" 2>/dev/null)
            
            if [ ! -z "$matches" ]; then
                sample_name=$(basename "$bed_file" _edits_trim.bed)
                echo "Found matches in $sample_name for $cluster_id"
                
                echo "$matches" >> bt_sample_edits/all_cluster_edits.bed
            fi
        fi
    done < bt_files.txt
    
done < cluster_list.txt

sort -k1,1 -k2,2n bt_sample_edits/all_cluster_edits.bed > bt_sample_edits/all_cluster_edits_sorted.bed
sort -u bt_sample_edits/all_cluster_edits_sorted.bed > bt_sample_edits/all_cluster_edits_unique.bed

awk 'BEGIN{OFS="\t"} 
NR>1 && NF>=6 {
    printf "%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $5, $8
}' bt_sample_edits/all_cluster_edits_unique.bed > bt_sample_edits/all_cluster_edits_clean.bed

#Extract only the transcript mappings for the clusters of interest
grep -f <(cut -d'_' -f1-2 cluster_list.txt) minimap2_results/transcript_locations_primary.bed > bt_sample_edits/relevant_transcript_mappings.bed

#For each edit, find the corresponding genomic location with proper exon mapping
awk 'BEGIN{OFS="\t"; print "#genomic_chr", "genomic_start", "genomic_end", "edit_name", "edit_score", "edit_strand", "transcript_chr", "transcript_start", "transcript_end", "cluster_name"}
# First pass: read transcript mappings and build cumulative exon positions
FNR==NR {
    cluster = $4
    genomic_chr = $1
    genomic_start = $2
    genomic_end = $3
    strand = $6
    
    # Store exon information for each cluster
    if (cluster in exon_count) {
        exon_count[cluster]++
    } else {
        exon_count[cluster] = 1
        transcript_length[cluster] = 0
    }
    
    exon_length = genomic_end - genomic_start
    exon_genomic_chr[cluster, exon_count[cluster]] = genomic_chr
    exon_genomic_start[cluster, exon_count[cluster]] = genomic_start
    exon_genomic_end[cluster, exon_count[cluster]] = genomic_end
    exon_strand[cluster, exon_count[cluster]] = strand
    exon_length_array[cluster, exon_count[cluster]] = exon_length
    
    exon_transcript_start[cluster, exon_count[cluster]] = transcript_length[cluster]
    exon_transcript_end[cluster, exon_count[cluster]] = transcript_length[cluster] + exon_length
    transcript_length[cluster] += exon_length
    
    next
}

# Function to get the correct exon index for transcript coordinate
function get_transcript_exon(cluster, transcript_pos) {
    if (exon_strand[cluster, 1] == "-") {
        # For negative strand: reverse mapping
        cumulative = 0
        for (i = exon_count[cluster]; i >= 1; i--) {
            exon_start = cumulative
            exon_end = cumulative + exon_length_array[cluster, i]
            if (transcript_pos >= exon_start && transcript_pos < exon_end) {
                return i
            }
            cumulative += exon_length_array[cluster, i]
        }
    } else {
        # For positive strand: direct mapping
        for (i = 1; i <= exon_count[cluster]; i++) {
            if (transcript_pos >= exon_transcript_start[cluster, i] && 
                transcript_pos < exon_transcript_end[cluster, i]) {
                return i
            }
        }
    }
    return 0
}

# Second pass: process edits and find which exon they belong to
{
    edit_cluster = $1
    edit_start = $2
    edit_end = $3
    edit_name = $4
    edit_score = $5
    edit_strand = $6
    
    if (edit_cluster in exon_count) {
        found_exon = 0
        
        if (exon_strand[edit_cluster, 1] == "-") {
            cumulative = 0
            for (i = exon_count[edit_cluster]; i >= 1; i--) {
                exon_start = cumulative
                exon_end = cumulative + exon_length_array[edit_cluster, i]
                
                if (edit_start >= exon_start && edit_start < exon_end) {
                    # Calculate position within this exon
                    position_in_exon = edit_start - exon_start
                    
                    # For negative strand: count backwards from exon end
                    genomic_edit_end = exon_genomic_end[edit_cluster, i] - position_in_exon
                    genomic_edit_start = genomic_edit_end - (edit_end - edit_start)
                    
                    print exon_genomic_chr[edit_cluster, i], genomic_edit_start, genomic_edit_end, edit_name, edit_score, edit_strand, edit_cluster, edit_start, edit_end, edit_cluster
                    found_exon = 1
                    break
                }
                cumulative += exon_length_array[edit_cluster, i]
            }
        } else {
            # For positive strand: direct mapping
            for (i = 1; i <= exon_count[edit_cluster]; i++) {
                if (edit_start >= exon_transcript_start[edit_cluster, i] && 
                    edit_start < exon_transcript_end[edit_cluster, i]) {
                    
                    position_in_exon = edit_start - exon_transcript_start[edit_cluster, i]
                    
                    # For positive strand: add to exon start
                    genomic_edit_start = exon_genomic_start[edit_cluster, i] + position_in_exon
                    genomic_edit_end = genomic_edit_start + (edit_end - edit_start)
                    
                    print exon_genomic_chr[edit_cluster, i], genomic_edit_start, genomic_edit_end, edit_name, edit_score, edit_strand, edit_cluster, edit_start, edit_end, edit_cluster
                    found_exon = 1
                    break
                }
            }
        }
        
        if (!found_exon) {
            print "Edit not found in any exon:", edit_cluster, edit_start, edit_end > "bt_sample_edits/edits_not_found_in_exons.txt"
        }
    } else {
        print "Cluster not found in mapping:", edit_cluster > "bt_sample_edits/edits_not_found_in_exons.txt"
    }
}' bt_sample_edits/relevant_transcript_mappings.bed bt_sample_edits/all_cluster_edits_clean.bed > bt_sample_edits/cluster_edits_genomic_coordinates.bed

# Check for unmapped edits
if [ -f bt_sample_edits/edits_not_found_in_exons.txt ]; then
    echo "WARNING: Some edits could not be mapped to genomic coordinates"
    echo "Unmapped edits: $(wc -l < bt_sample_edits/edits_not_found_in_exons.txt)"
    echo "Check bt_sample_edits/edits_not_found_in_exons.txt for details"
else
    echo "All edits successfully mapped to genomic coordinates"
fi
