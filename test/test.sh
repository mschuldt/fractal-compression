#!/bin/bash

run () {
    # convert to raw format
    convert -depth 8 -size "$2" $1 in.rgb
    # Compress the image
    time ./../fractal -r -o 1 -t 100 -p 5 in.rgb
    # change extension to .rgb for convert utility
    mv output.raw output.rgb
    # convert raw image back to .jpg
    convert -depth 8 -size "$2" output.rgb $1_out.jpg
    # remove temp files
    rm in.rgb output.rgb
}

check (){
    ref_file=$1_compressed_ref.jpg
    out_file=$1_out.jpg
    cmp -s $ref_file $out_file > /dev/null
    if [ $? -eq 1 ]; then
        compare $ref_file $out_file $1_diff.jpg
        echo "FAIL"
    else
        echo "OK"
    fi
}

run lena256.jpg "256x256"
check lena256.jpg
