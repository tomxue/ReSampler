#!/usr/bin/env bash

# notes on obtaining spectrogram utility:
# spectrogram tool available from http://www.mega-nerd.com/libsndfile/tools/#spectrogram
# ... and https://github.com/erikd/sndfile-tools/tree/master/src

# Windows binary for spectrogram.exe available here: https://hydrogenaud.io/index.php/topic,102495.0.html
# Windows dependencies: libcairo-2.dll libpng12-0.dll libsndfile-1.dll libfftw3-3.dll

spectrogram_tool="spectrogram"

echo $(tput setaf 2)cleaning ./spectrograms folder ...$(tput setaf 7)
rm ./spectrograms/*.*

echo $(tput setaf 2)generating spectrograms ...$(tput setaf 1)
find ./outputs -type f ! -name '*.png' -exec $spectrogram_tool --dyn-range=190 '{}' 1200 960 '{}'.png \;

echo $(tput setaf 2)moving spectrograms to ./spectrograms folder$(tput setaf 7)
mv ./outputs/*.png ./spectrograms