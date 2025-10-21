#identify consensus RNA editing targets across multiple analysis
library(pheatmap)
library(RColorBrewer)
library(dplyr)        # For data manipulation
library(tidyr)        # For data reshaping
library(VennDiagram)  # For venn diagrams (if you plan to visualize overlaps)
library(UpSetR) 

input_files <- c(
    # v1
    "/Volumes/Lucky_me/9_STAMP_transcript/05.transcript_analysis/x/summed_edit_rates_with_genes_v2_top200.tsv",
    # v2 group
    "/Volumes/Lucky_me/9_STAMP_transcript/06.transcript_analysis_v2/processed_top200_ordered.tsv",
    "/Volumes/Lucky_me/9_STAMP_transcript/06.transcript_analysis_v2/processed_top200_ordered_apo.tsv",
    "/Volumes/Lucky_me/9_STAMP_transcript/06.transcript_analysis_v2/processed_top200_ordered_tad.tsv",
    # v3 group
    "/Volumes/Lucky_me/9_STAMP_transcript/07.transcript_analysis_v3/1_all/processed_top200_ordered.tsv",
    "/Volumes/Lucky_me/9_STAMP_transcript/07.transcript_analysis_v3/2_apo/processed_top200_ordered.tsv",
    "/Volumes/Lucky_me/9_STAMP_transcript/07.transcript_analysis_v3/3_tad/processed_top200_ordered.tsv"
)

# Extract transcripts from each file
process_file <- function(file_path) {
    ht <- read.delim(file_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    ht <- head(ht, 200)  # Get top 200 rows
    
    # Extract clusters and genes based on file type
    clusters <- ht[, 1]
    if (grepl("processed_top200", file_path)) {
        genes <- ht[, ncol(ht)-1]  # Second to last column for processed files
    } else {
        genes <- ht[, ncol(ht)]    # Last column for summed_edit_rates files
    }
    
    # Combine cluster and gene information
    transcripts <- paste(clusters, genes, sep="_")
    return(transcripts)
}

transcript_lists <- lapply(input_files, process_file)
names(transcript_lists) <- c("v1_all",                # v1
                           "v2_all", "v2_apo", "v2_tad",  # v2 group
                           "v3_all", "v3_apo", "v3_tad")# Find common transcripts
common_transcripts <- Reduce(intersect, transcript_lists)

write.table(data.frame(transcript=common_transcripts),
            "common_transcripts_all_versions.tsv",
            sep="\t", row.names=FALSE, quote=FALSE)

cat("Number of common transcripts:", length(common_transcripts), "\n")

# Calculate pairwise overlaps
n_files <- length(transcript_lists)
versions <- c("v1_all", "v2_all", "v2_apo", "v2_tad", "v3_all", "v3_apo", "v3_tad")
overlap_matrix <- matrix(0, n_files, n_files)
rownames(overlap_matrix) <- versions
colnames(overlap_matrix) <- versions

for(i in 1:n_files) {
    overlap_matrix[i,i] <- length(transcript_lists[[versions[i]]])
}

# Fill pairwise overlaps
for(i in 1:(n_files-1)) {
    for(j in (i+1):n_files) {
        overlap <- length(intersect(transcript_lists[[versions[i]]], transcript_lists[[versions[j]]]))
        overlap_matrix[i,j] <- overlap
        overlap_matrix[j,i] <- overlap
    }
}

write.table(overlap_matrix,
            "pairwise_overlaps_top200.tsv",
            sep="\t", 
            quote=FALSE,
            col.names=NA)

# Visualize pairwise overlaps as heatmap
png("overlap_heatmap_top200.png", width=2000, height=2000, res=300)
pheatmap(overlap_matrix,
         display_numbers=TRUE,  # Show numbers in cells
         number_format="%.0f",  # No decimal places
         cluster_rows=FALSE,    # Don't cluster
         cluster_cols=FALSE,
         color=colorRampPalette(c("white", "red"))(100),
         main="Pairwise Overlaps Between Versions",
         fontsize=14,
         fontsize_number=12)
dev.off()

library(UpSetR)

# Convert transcript lists to binary matrix format
all_transcripts <- unique(unlist(transcript_lists))
upset_matrix <- matrix(0, nrow=length(all_transcripts), ncol=length(transcript_lists))
colnames(upset_matrix) <- names(transcript_lists)
rownames(upset_matrix) <- all_transcripts

for(set_name in names(transcript_lists)) {
    upset_matrix[, set_name] <- all_transcripts %in% transcript_lists[[set_name]]
}
upset_matrix <- as.data.frame(upset_matrix)

png("upset_plot_top200.png", width=3000, height=2000, res=300)
upset(upset_matrix, 
      sets = colnames(upset_matrix),
      sets.bar.color = "darkred",
      main.bar.color = "darkblue",
      text.scale = 1.5,
      point.size = 3,
      line.size = 1,
      order.by = "freq")
dev.off()

# Create a presence/absence matrix for all transcripts
all_transcripts <- unique(unlist(transcript_lists))
presence_matrix <- matrix(FALSE, 
                         nrow = length(all_transcripts), 
                         ncol = length(transcript_lists))
colnames(presence_matrix) <- names(transcript_lists)
rownames(presence_matrix) <- all_transcripts

for(set_name in names(transcript_lists)) {
    presence_matrix[, set_name] <- all_transcripts %in% transcript_lists[[set_name]]
}

overlap_df <- as.data.frame(presence_matrix)
overlap_df$transcript <- rownames(overlap_df)
overlap_df$present_in_n_sets <- rowSums(presence_matrix)
overlap_df$present_in_sets <- apply(presence_matrix, 1, function(x) {
    paste(names(transcript_lists)[x], collapse=";")
})

overlap_df <- overlap_df[order(-overlap_df$present_in_n_sets), ]
write.table(head(overlap_df, 2000),
            "transcript_overlaps_top200.tsv",
            sep="\t",
            row.names=FALSE,
            quote=FALSE)

cat("\nNumber of transcripts by overlap count:\n")
print(table(overlap_df$present_in_n_sets))
