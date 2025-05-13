#!/bin/bash

# convert datset of h5 to ftr 
./autoexec_convert_dataset.sh

# convert otest to ftr
./autoexec_convert_otest.sh

# convert ground truth to answer file
./autoexec_gt_otest_to_answer.sh

# make bucket file for sketch
./autoexec_make_bucket.sh sketch_w16
./autoexec_make_bucket.sh sketch_w18
./autoexec_make_bucket.sh sketch_w20

# search by double filtering
TEMP_DIR="temp"
if [ -d "$TEMP_DIR" ] && [ -n "$(ls -A "$TEMP_DIR")" ]; then
  rm -rf "$TEMP_DIR"/*
fi
./autoexec_eval.sh

# copy result and summary
cp temp/result*.csv result/
cat temp/summary*.csv > result/summary_all.csv
