#!/bin/bash

# usage:
#  ./run.sh <filename> <size>
# 
# example:
#  ./run.sh lena.jpg 256x256
#
# decompressed output file is saved as 'output.jpg'

# convert to raw format
convert -depth 8 -size "$2" $1 in.rgb
# Compress the image
./fractal -r -o 1 -t 100 -p 5 in.rgb
# change extension to .rgb for convert utility
mv output.raw output.rgb
# convert raw image back to .jpg
convert -depth 8 -size "$2" output.rgb  output.jpg
# remove temp files
rm in.rgb output.rgb
