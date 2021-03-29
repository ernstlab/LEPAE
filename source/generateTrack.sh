#!/usr/bin/env bash


a_data_filename=$1
b_data_filename=$2
score_filename=$3
score_threshold=$4
score_name=$5
output_filename=$6

# echo $output_filename

printf "track type=interact name="$score_name" " > $output_filename
printf "useScore=on maxHeightPixels=200:100:50 visibility=full scoreFilter="$score_threshold"\n" >> $output_filename
printf "browser position chr1:0-10,000,000\n" >> $output_filename

paste -d "\t" \
<(gzip -cd $a_data_filename | cut -f 1-3) \
<(gzip -cd $b_data_filename | cut -f 1-3) \
<(gzip -cd $score_filename | paste -d"\t" - -) |\
bedtools sort |\
awk -v OFS="\t" '{ printf "%s\t%d\t%d\t.\t%d\t%.5f\t.\t0\t%s\t%d\t%d\t.\t.\t%s\t%d\t%d\t.\t.\n", $1, $2, $6, ($7+$8)*500, ($7+$8)/2, $1, $2, $3, $4, $5, $6 }' >> $output_filename

# ./generateTrack.sh ../data/chr1_dist10000.base.a.gz ../data/chr1_dist10000.base.b.gz ../prediction/chr1_dist10000_score.txt.gz 950 pairwise_dist10kb ../prediction/chr1_dist10000.interact; ./generateTrack.sh ../data/chr1_dist20000.base.a.gz ../data/chr1_dist20000.base.b.gz ../prediction/chr1_dist20000_score.txt.gz 950 pairwise_dist20kb ../prediction/chr1_dist20000.interact; ./generateTrack.sh ../data/chr1_dist30000.base.a.gz ../data/chr1_dist30000.base.b.gz ../prediction/chr1_dist30000_score.txt.gz 950 pairwise_dist30kb ../prediction/chr1_dist30000.interact