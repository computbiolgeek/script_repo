#!/bin/bash

grep 'completeness_mem' $1/*.pdb | awk '$5>=0.9 && $5 != "nan" {print}' | awk -F: '{print $1}' > $1/$1_complete_mem.ls
