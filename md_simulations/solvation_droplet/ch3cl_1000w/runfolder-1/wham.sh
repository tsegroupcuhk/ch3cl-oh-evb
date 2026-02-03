#!/bin/bash

python python get_CV_corr_time_edit_metafile.py metafile metafile-trun
wham 0 40 41 0.0000001 300.0 0 metafile-trun WHAM_output 200 1234

