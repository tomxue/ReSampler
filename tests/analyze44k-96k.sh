#!/usr/bin/env bash

input_path=./inputs
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

$resampler_path -i ./inputs/44khz_sweep-3dBFS_32f.wav -o ./outputs/44khz_sweep-3dBFS_32f-to96k-dp-final.wav -r 96000 --doubleprecision --showstages

#$resampler_path -i ./inputs/44khz_sweep-3dBFS_32f.wav -o ./outputs/44khz_sweep-3dBFS_32f-to96k-dp-stage1.wav -r 73500 --doubleprecision --lpf-cutoff 90.909091 --lpf-transition 5.454545 --maxStages 1
#$resampler_path -i ./outputs/44khz_sweep-3dBFS_32f-to96k-dp-stage1.wav -o ./outputs/44khz_sweep-3dBFS_32f-to96k-dp-stage2.wav -r 84000 --doubleprecision --lpf-cutoff 90.909091 --lpf-transition 74.772727 --maxStages 1
#$resampler_path -i ./outputs/44khz_sweep-3dBFS_32f-to96k-dp-stage2.wav -o ./outputs/44khz_sweep-3dBFS_32f-to96k-dp.wav -r 96000 --doubleprecision --lpf-cutoff 90.909091 --lpf-transition 45.738636 --maxStages 1

$resampler_path -i ./inputs/44khz_sweep-3dBFS_32f.wav -o ./outputs/44khz_sweep-3dBFS_32f-to96k-dp-stage1.wav -r 73500 --doubleprecision --lpf-cutoff 90.909091 --lpf-transition 5.454545 --maxStages 1
$resampler_path -i ./outputs/44khz_sweep-3dBFS_32f-to96k-dp-stage1.wav -o ./outputs/44khz_sweep-3dBFS_32f-to96k-dp-stage2.wav -r 84000 --doubleprecision --lpf-cutoff 54.5442 --lpf-transition 74.772727 --maxStages 1
$resampler_path -i ./outputs/44khz_sweep-3dBFS_32f-to96k-dp-stage2.wav -o ./outputs/44khz_sweep-3dBFS_32f-to96k-dp.wav -r 96000 --doubleprecision --lpf-cutoff 41.7604 --lpf-transition 45.738636 --maxStages 1