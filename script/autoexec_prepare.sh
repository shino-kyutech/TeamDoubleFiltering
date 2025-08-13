#!/bin/bash

# convert datset of h5 to ftr 
time ./autoexec_convert_dataset.sh

# convert otest to ftr
time ./autoexec_convert_otest.sh

# convert ground truth to answer file
time ./autoexec_gt_otest_to_answer.sh

# make bucket file for sketch
time ./autoexec_make_bucket.sh sketch_w16
time ./autoexec_make_bucket.sh sketch_w18
time ./autoexec_make_bucket.sh sketch_w20
