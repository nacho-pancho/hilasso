#!/bin/bash
for file in mnist_mat.7z usps_mat.7z
do
  wget -c http://iie.fing.edu.uy/~nacho/data/${file}
  7zr x ${file}
done
