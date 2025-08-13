#!/bin/bash

# prepare for search
./autoexec_prepare.sh

# search by double filtering
TEMP_DIR="temp"
if [ -d "$TEMP_DIR" ] && [ -n "$(ls -A "$TEMP_DIR")" ]; then
  rm -rf "$TEMP_DIR"/*
fi

./autoexec_eval.sh temp NONE v5_2 8

# copy result and summary
cp temp/result*.csv result/
cat temp/summary*.csv > result/summary_all.csv
