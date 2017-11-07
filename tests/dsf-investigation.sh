#!/usr/bin/env bash

input_path=~
output_path=./outputs

function tolower(){
    echo $1 | sed "y/ABCDEFGHIJKLMNOPQRSTUVWXYZ/abcdefghijklmnopqrstuvwxyz/"
}

os=`tolower $OSTYPE`

# set converter path according to OS:
if [ $os == 'cygwin' ] || [ $os == 'msys' ]
then 
    #Windows ...
    resampler_path=../x64/Release/ReSampler.exe
    #resampler_path="E:\Temp\ReSampler1.3.6\ReSampler.exe"
    #resampler_path=../x64/minGW-W64/ReSampler.exe
else
    resampler_path=../ReSampler
fi

# clear old outputs:
rm $output_path/*.*
rm $output_path/._* 

$resampler_path --single-stage -b 32f -i /Users/newuser/sweep-176400hz-0-22050hz-20s-D64-2.8mhz.dsf -o ./outputs/96khz_sweep-3dBFS_32f-to44k-dp-stage1.wav -r 705600 --doubleprecision --lpf-cutoff 5.681818 --lpf-transition 94.318182 --maxStages 1
$resampler_path --single-stage -b 32f -i ./outputs/96khz_sweep-3dBFS_32f-to44k-dp-stage1.wav -o ./outputs/96khz_sweep-3dBFS_32f-to44k-dp-stage2.wav -r 176400 --doubleprecision --lpf-cutoff 22.727273 --lpf-transition 77.272727 --maxStages 1
$resampler_path --single-stage -b 32f -i ./outputs/96khz_sweep-3dBFS_32f-to44k-dp-stage2.wav -o ./outputs/96khz_sweep-3dBFS_32f-to44k-dp.wav -r 44100 --doubleprecision --lpf-cutoff 90.909091 --lpf-transition 9.090909 --maxStages 1

$resampler_path --single-stage -b 32f -i /Users/newuser/sweep-176400hz-0-22050hz-20s-D64-2.8mhz.dsf -o ./outputs/96khz_sweep-3dBFS_32f-to44k-dp-ss.wav -r 44100 --doubleprecision --steepLPF --maxStages 1

# 32-bit float, double-precision sweeps
#$resampler_path -i $input_path/sweep-176400hz-0-22050hz-20s-D64-2.8mhz.dsf --showStages -o $output_path/96khz_sweep-3dBFS_32f-to22k-dp.wav -b 32f -r 22050 --noPeakChunk --doubleprecision
#$resampler_path -i $input_path/sweep-176400hz-0-22050hz-20s-D64-2.8mhz.dsf --showStages -o $output_path/96khz_sweep-3dBFS_32f-to32k-dp.wav -b 32f -r 32000 --noPeakChunk --doubleprecision --mt
#$resampler_path -i $input_path/sweep-176400hz-0-22050hz-20s-D64-2.8mhz.dsf --showStages -o $output_path/96khz_sweep-3dBFS_32f-to32k-dp.wav -b 32f -r 32000 --noPeakChunk --doubleprecision --mt
#$resampler_path -i $input_path/sweep-176400hz-0-22050hz-20s-D64-2.8mhz.dsf --showStages -o $output_path/96khz_sweep-3dBFS_32f-to44k-dp.wav -b 32f -r 44100 --noPeakChunk --doubleprecision
#$resampler_path -i $input_path/sweep-176400hz-0-22050hz-20s-D64-2.8mhz.dsf --showStages -o $output_path/96khz_sweep-3dBFS_32f-to48k-dp.wav -b 32f -r 48000 --noPeakChunk --doubleprecision --mt
#$resampler_path -i $input_path/sweep-176400hz-0-22050hz-20s-D64-2.8mhz.dsf --showStages -o $output_path/96khz_sweep-3dBFS_32f-to88k-dp.wav -b 32f -r 88200 --noPeakChunk --doubleprecision
#$resampler_path -i $input_path/sweep-176400hz-0-22050hz-20s-D64-2.8mhz.dsf --showStages -o $output_path/96khz_sweep-3dBFS_32f-to96k-dp.wav -b 32f -r 96000 --noPeakChunk --doubleprecision --mt
#$resampler_path -i $input_path/sweep-176400hz-0-22050hz-20s-D64-2.8mhz.dsf --showStages -o $output_path/96khz_sweep-3dBFS_32f-to176k-dp.wav -b 32f -r 176400 --noPeakChunk --doubleprecision
#$resampler_path -i $input_path/sweep-176400hz-0-22050hz-20s-D64-2.8mhz.dsf --showStages -o $output_path/96khz_sweep-3dBFS_32f-to192k-dp.wav -b 32f -r 192000 --noPeakChunk --doubleprecision --mt

# make spectrograms
./make-spectrograms.sh