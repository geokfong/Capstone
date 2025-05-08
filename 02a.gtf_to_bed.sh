awk -F"\t" '{if(FNR == NR){hash[$2]=$1}else{if($0 !~ /#/){split($9, info, ";"); OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8,hash[" "info[1]]"; "$9}}}' ENSG_TRANSID.txt introns.gtf > Introns_wnames.gtf

cat Introns_wnames.gtf Transcripts_to_work_on.gtf > combined.gtf
