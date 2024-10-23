#!/usr/bin/env bash

DEF_DIR=$(find . -name "*.def" -exec dirname {} \;)
singularity build --fakeroot "${DEF_DIR}/singularity.sif" "${DEF_DIR}/singularity_image.def"