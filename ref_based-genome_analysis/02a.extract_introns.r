library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)

gtf_file <- "Transcripts_to_work_on.gtf"
gtf <- import(gtf_file)
txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")
introns <- intronsByTranscript(txdb, use.names=TRUE)
introns_gr <- unlist(introns)
mcols(introns_gr)$type <- "intron"
mcols(introns_gr)$transcript_id <- names(introns_gr)

export(introns_gr, "introns.gtf", format="gtf")
