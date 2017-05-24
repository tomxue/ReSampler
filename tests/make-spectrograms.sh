#!/usr/bin/env bash

echo $(tput setaf 2)cleaning ./spectrograms folder ...$(tput setaf 7)
rm ./spectrograms/*.*

echo $(tput setaf 2)generating spectrograms ...$(tput setaf 1)
find ./outputs -type f ! -name '*.png' -exec spectrogram --dyn-range=190 '{}' 1200 960 '{}'.png \;

echo $(tput setaf 2)moving spectrograms to ./spectrograms folder$(tput setaf 7)
mv ./outputs/*.png ./spectrograms