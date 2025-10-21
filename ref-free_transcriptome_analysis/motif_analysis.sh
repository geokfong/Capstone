bedtools getfasta -fi  /Volumes/Lucky_me/10_genome_transcript_analysis/genome.fa -bed /Volumes/Lucky_me/10_genome_transcript_analysis/bt_sample_edits/bt_intron_3utr_only.bed -fo bt_intron_3utr_only.fa -s -name
bedtools getfasta -fi  /Volumes/Lucky_me/10_genome_transcript_analysis/genome.fa -bed /Volumes/Lucky_me/8_STAMP/genemodels/temporary/selected_regions.bed -fo introns_and_utrs.fa -s -name
bedtools getfasta -fi /Volumes/Lucky_me/10_genome_transcript_analysis/genome.fa -bed /Volumes/Lucky_me/ori_regions.bed -fo ori_regions.fa -s -name

fasta-get-markov -m 2 introns_and_utrs.fa bg_model.txt

meme ori_regions.fa \
    -dna -mod zoops -nmotifs 5 -minw 6 -maxw 30 \
    -bfile /Volumes/Lucky_me/10_genome_transcript_analysis/test/playground_2/bg_model.txt -oc meme_out_ori

meme introns_and_utrs.fa \
    -dna -mod zoops -nmotifs 3 -minw 6 -maxw 30 \
    -oc meme_out_neg
 
fimo --oc fimo_high_ori /Volumes/Lucky_me/meme_out_ori/meme.txt /Volumes/Lucky_me/10_genome_transcript_analysis/test/playground_2/bt_intron_3utr_only.fa
fimo --oc fimo_low_ori /Volumes/Lucky_me/meme_out_ori/meme.txt /Volumes/Lucky_me/10_genome_transcript_analysis/test/playground_2/introns_and_utrs.fa

awk 'BEGIN{OFS="\t"} NR>1 {print $3, $4, $5, $1, $7, $6, $9}' fimo.tsv > fimo.bed
