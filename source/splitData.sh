#!/usr/bin/env bash

input_filename=$1
output_filename_prefix=$2
seed=$3

get_seeded_random()
{
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt \
    </dev/zero 2>/dev/null;
}

gzip -cd $input_filename | awk '$4 ~ /a/' | sort -k4 | shuf --random-source=<(get_seeded_random) | gzip > "$output_filename_prefix".a.gz
gzip -cd $input_filename | awk '$4 ~ /b/' | sort -k4 | shuf --random-source=<(get_seeded_random) | gzip > "$output_filename_prefix".b.gz