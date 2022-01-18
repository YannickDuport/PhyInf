#!/bin/bash

path_to_repo=$1

`$path_to_repo/org_script.py find -q "NUM_TAXA > 100 and DATA_TYPE = 'DNA' and OVERALL_NUM_PARTITIONS = 1"`

