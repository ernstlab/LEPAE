#!/usr/bin/env bash

input_filename=$1
seed=$2
output_filename=$3

get_seeded_random()
{
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt \
    </dev/zero 2>/dev/null;
}

gzip -cd $input_filename | shuf --random-source=<(get_seeded_random) | gzip > $output_filename