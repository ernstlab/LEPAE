#!/usr/bin/env bash

chrom_size_filename=$1
distance=$2

tuning_chrom=$3
test_chrom=$4
test_start=$5
test_end=$6
test_window_size=$7

training_data_size=$8
tuning_data_size=$9

seed=${10}

training_position_output_filename_prefix=${11}
tuning_position_output_filename_prefix=${12}
test_position_output_filename_prefix=${13}
prediction_position_output_filename_prefix=${14}



# Training pairs
cat $chrom_size_filename | awk -v tuning_chrom=$tuning_chrom -v test_chrom=$test_chrom '$1!=tuning_chrom && $1!=test_chrom' \
> "$training_position_output_filename_prefix".tmp

bedtools random -l $distance -n $training_data_size -g "$training_position_output_filename_prefix".tmp -seed $seed |\
cut -f 1-3 | bedtools sort | awk -v OFS="\t" '{print $1,$2,$3,NR}' |\
gzip > "$training_position_output_filename_prefix".pair.gz 

gzip -cd "$training_position_output_filename_prefix".pair.gz  | awk -v OFS="\t" '{print $1,$2,$2+1,$4"a","\n"$1,$3,$3+1,$4"b"}' |\
bedtools sort | gzip > "$training_position_output_filename_prefix".base.gz 

rm "$training_position_output_filename_prefix".tmp



# Tuning pairs
cat $chrom_size_filename | awk -v tuning_chrom=$tuning_chrom '$1==tuning_chrom' \
> "$tuning_position_output_filename_prefix".tmp

bedtools random -l $distance -n $tuning_data_size -g "$tuning_position_output_filename_prefix".tmp -seed $seed |\
cut -f 1-3 | bedtools sort | awk -v OFS="\t" '{print $1,$2,$3,NR}' |\
gzip > "$tuning_position_output_filename_prefix".pair.gz 

gzip -cd "$tuning_position_output_filename_prefix".pair.gz   | awk -v OFS="\t" '{print $1,$2,$2+1,$4"a","\n"$1,$3,$3+1,$4"b"}' |\
bedtools sort | gzip > "$tuning_position_output_filename_prefix".base.gz 

rm "$tuning_position_output_filename_prefix".tmp



# Test pairs
cat $chrom_size_filename | awk -v test_chrom=$test_chrom '$1==test_chrom' \
> "$test_position_output_filename_prefix".tmp

bedtools random -l $distance -n $tuning_data_size -g "$test_position_output_filename_prefix".tmp -seed $seed |\
cut -f 1-3 | bedtools sort | awk -v OFS="\t" '{print $1,$2,$3,NR}' |\
gzip > "$test_position_output_filename_prefix".pair.gz 

gzip -cd "$test_position_output_filename_prefix".pair.gz   | awk -v OFS="\t" '{print $1,$2,$2+1,$4"a","\n"$1,$3,$3+1,$4"b"}' |\
bedtools sort | gzip > "$test_position_output_filename_prefix".base.gz 

rm "$test_position_output_filename_prefix".tmp



# Prediction pairs
awk -v start=$test_start -v end=$test_end -v test_window_size=$test_window_size -v chrom=$test_chrom 'BEGIN {for (i=start; i<end; i+=test_window_size) print chrom,i,i+test_window_size}' |\
awk -v OFS="\t" -v end=$test_end '$3<end {print $1,$2,$3,NR}' | gzip > "$prediction_position_output_filename_prefix".pair.gz

gzip -cd "$prediction_position_output_filename_prefix".pair.gz | awk -v OFS="\t" '{print $1,$2,$2+1,$4"a","\n"$1,$3,$3+1,$4"b"}' |\
bedtools sort | gzip > "$prediction_position_output_filename_prefix".base.gz