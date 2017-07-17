#!/bin/bash
for i in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
do 
  wget http://r0k.us/graphics/kodak/kodak/kodim${i}.png
  convert kodim${i}.png kodim${i}.ppm
done
