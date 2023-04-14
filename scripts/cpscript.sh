#!/bin/bash

for i in ls /mnt/data/research_data/2021-03-29_rrbs61-65_peace_hg_himani_mm/usftp21.novogene.com/raw_data/*
    do                 # Line breaks are important
        if [ -d $i ]   # Spaces are important
            then
                ln -s $i/*.gz .
        fi
    done