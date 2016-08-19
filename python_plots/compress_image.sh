#!/bin/bash

i=$1
n=`basename $i .tif`
convert $i -set colorspace RGB -layers flatten -alpha off -compress lzw -depth 8 \
-density 600 -adaptive-resize 4500x2400  images/$n"_compressed".tif

convert images/$n"_compressed".tif images/$n"_compressed".pdf
