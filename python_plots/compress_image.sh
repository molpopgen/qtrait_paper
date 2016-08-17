#!/bin/bash

for i in *.tif
do
    n=`basename $i .tif`
    convert $i -set colorspace RGB -layers flatten -alpha off -compress lzw -depth 8 \
    -density 600 -adaptive-resize 4500x2400  $n"_compressed".tif

    convert $n"_compressed".tif $n"_compressed".pdf
done
