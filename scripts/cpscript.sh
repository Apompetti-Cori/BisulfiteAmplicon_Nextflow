#!/bin/bash

for i in ls /mnt/data/research_data/2021-12-07_mouse_rrbs_peace/usftp21.novogene.com/raw_data/*
    do                 # Line breaks are important
        if [ -d $i ]   # Spaces are important
            then
                ln -s $i/*.gz .
        fi
    done