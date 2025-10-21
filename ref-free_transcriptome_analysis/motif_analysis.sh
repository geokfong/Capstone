bedtools getfasta -fi  /Volumes/Lucky_me/10_genome_transcript_analysis/genome.fa -bed /Volumes/Lucky_me/10_genome_transcript_analysis/bt_sample_edits/bt_intron_3utr_only.bed -fo bt_intron_3utr_only.fa -s -name
bedtools getfasta -fi  /Volumes/Lucky_me/10_genome_transcript_analysis/genome.fa -bed /Volumes/Lucky_me/8_STAMP/genemodels/temporary/selected_regions.bed -fo introns_and_utrs.fa -s -name

fasta-get-markov -m 2 introns_and_utrs.fa bg_model.txt

meme bt_intron_3utr_only.fa \
    -dna -mod zoops -nmotifs 3 -minw 6 -maxw 30 \
    -bfile bg_model.txt -oc meme_out_3

fimo --oc fimo_high_v3 /Volumes/Lucky_me/10_genome_transcript_analysis/test/playground_2/v3/meme_out_3/meme.txt /Volumes/Lucky_me/10_genome_transcript_analysis/test/playground_2/bt_intron_3utr_only.fa
fimo --oc fimo_low_v3 /Volumes/Lucky_me/10_genome_transcript_analysis/test/playground_2/v3/meme_out_3/meme.txt /Volumes/Lucky_me/10_genome_transcript_analysis/test/playground_2/introns_and_utrs.fa

meme high.fa -oc meme_out -dna -mod zoops -nmotifs 5 -minw 6 -maxw 30

awk 'BEGIN{OFS="\t"} NR>1 {print $3, $4, $5, $1, $7, $6, $9}' fimo.tsv > fimo.bed

meme introns_and_utrs.fa \
    -dna -mod zoops -nmotifs 3 -minw 6 -maxw 30 \
    -oc meme_out_neg
