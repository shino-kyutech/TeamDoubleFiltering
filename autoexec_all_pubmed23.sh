#!/bin/bash

# convert datset of h5 to ftr 
./autoexec_convert_dataset_pubmed23.sh

# convert otest to ftr
./autoexec_convert_otest_pubmed23.sh

# convert ground truth to answer file
./autoexec_gt_otest_to_answer_pubmed23.sh

# make bucket file for sketch
./autoexec_make_bucket.sh sketch_w21

# search by double filtering
./search_by_double_filtering_pubmed23.sh sketch_w21 QPSMAP_p192 00_23 otest NONE NONE temp 30 0 1 2.0 0.5 10 25 8 7 13 3 1.05 100 0 8 1 1 df_parameters_w21_p192_nc2_1000.txt summary.csv NONE
#./search_by_double_filtering_pubmed23.sh sketch_w21 QPSMAP_p192 00_23 otest NONE NONE result.csv 30 0 2 2.0 0.5 10 25 8 7 13 3 1.05 100 0 8 1 1 df_parameters_w21_p192_2bit.txt result/summary.csv result/search_cost.csv
cp temp/*.csv result/
