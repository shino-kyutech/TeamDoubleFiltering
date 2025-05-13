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
