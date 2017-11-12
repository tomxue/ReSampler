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

# 32-bit float, double-precision sweeps
#$resampler_path -i $input_path/sweep-176400hz-0-22050hz-20s-D64-2.8mhz.dsf -o $output_path/96khz_sweep-3dBFS_32f-to22k-dp.wav -b 32f -r 22050 --noPeakChunk --doubleprecision
#$resampler_path -i $input_path/sweep-176400hz-0-22050hz-20s-D64-2.8mhz.dsf -o $output_path/96khz_sweep-3dBFS_32f-to32k-dp.wav -b 32f -r 32000 --noPeakChunk --doubleprecision --mt
#$resampler_path -i $input_path/sweep-176400hz-0-22050hz-20s-D64-2.8mhz.dsf -o $output_path/96khz_sweep-3dBFS_32f-to32k-dp.wav -b 32f -r 32000 --noPeakChunk --doubleprecision --mt
#$resampler_path -i $input_path/sweep-176400hz-0-22050hz-20s-D64-2.8mhz.dsf -o $output_path/96khz_sweep-3dBFS_32f-to44k-dp.wav -b 32f -r 44100 --noPeakChunk --doubleprecision
#$resampler_path -i $input_path/sweep-176400hz-0-22050hz-20s-D64-2.8mhz.dsf -o $output_path/96khz_sweep-3dBFS_32f-to48k-dp.wav -b 32f -r 48000 --noPeakChunk --doubleprecision --mt
#$resampler_path -i $input_path/sweep-176400hz-0-22050hz-20s-D64-2.8mhz.dsf -o $output_path/96khz_sweep-3dBFS_32f-to88k-dp.wav -b 32f -r 88200 --noPeakChunk --doubleprecision
#$resampler_path -i $input_path/sweep-176400hz-0-22050hz-20s-D64-2.8mhz.dsf -o $output_path/96khz_sweep-3dBFS_32f-to96k-dp.wav -b 32f -r 96000 --noPeakChunk --doubleprecision --mt
#$resampler_path -i $input_path/sweep-176400hz-0-22050hz-20s-D64-2.8mhz.dsf -o $output_path/96khz_sweep-3dBFS_32f-to176k-dp.wav -b 32f -r 176400 --noPeakChunk --doubleprecision
#$resampler_path -i $input_path/sweep-176400hz-0-22050hz-20s-D64-2.8mhz.dsf -o $output_path/96khz_sweep-3dBFS_32f-to192k-dp.wav -b 32f -r 192000 --noPeakChunk --doubleprecision --mt

# 32-bit float, double-precision sweeps
$resampler_path --singlestage -i $input_path/sweep-176400hz-0-22050hz-20s-D64-2.8mhz.dsf -o $output_path/96khz_sweep-3dBFS_32f-to22k-SSdp.wav -b 32f -r 22050 --noPeakChunk --doubleprecision
$resampler_path --singlestage -i $input_path/sweep-176400hz-0-22050hz-20s-D64-2.8mhz.dsf -o $output_path/96khz_sweep-3dBFS_32f-to32k-SSdp.wav -b 32f -r 32000 --noPeakChunk --doubleprecision --mt
$resampler_path --singlestage -i $input_path/sweep-176400hz-0-22050hz-20s-D64-2.8mhz.dsf -o $output_path/96khz_sweep-3dBFS_32f-to32k-SSdp.wav -b 32f -r 32000 --noPeakChunk --doubleprecision --mt
$resampler_path --singlestage -i $input_path/sweep-176400hz-0-22050hz-20s-D64-2.8mhz.dsf -o $output_path/96khz_sweep-3dBFS_32f-to44k-SSdp.wav -b 32f -r 44100 --noPeakChunk --doubleprecision
$resampler_path --singlestage -i $input_path/sweep-176400hz-0-22050hz-20s-D64-2.8mhz.dsf -o $output_path/96khz_sweep-3dBFS_32f-to48k-SSdp.wav -b 32f -r 48000 --noPeakChunk --doubleprecision --mt
$resampler_path --singlestage -i $input_path/sweep-176400hz-0-22050hz-20s-D64-2.8mhz.dsf -o $output_path/96khz_sweep-3dBFS_32f-to88k-SSdp.wav -b 32f -r 88200 --noPeakChunk --doubleprecision
$resampler_path --singlestage -i $input_path/sweep-176400hz-0-22050hz-20s-D64-2.8mhz.dsf -o $output_path/96khz_sweep-3dBFS_32f-to96k-SSdp.wav -b 32f -r 96000 --noPeakChunk --doubleprecision --mt
$resampler_path --singlestage -i $input_path/sweep-176400hz-0-22050hz-20s-D64-2.8mhz.dsf -o $output_path/96khz_sweep-3dBFS_32f-to176k-SSdp.wav -b 32f -r 176400 --noPeakChunk --doubleprecision
$resampler_path --singlestage -i $input_path/sweep-176400hz-0-22050hz-20s-D64-2.8mhz.dsf -o $output_path/96khz_sweep-3dBFS_32f-to192k-SSdp.wav -b 32f -r 192000 --noPeakChunk --doubleprecision --mt

# 32-bit float, double-precision sweeps
#$resampler_path --maxStages 1 -i $input_path/sweep-176400hz-0-22050hz-20s-D64-2.8mhz.dsf -o $output_path/96khz_sweep-3dBFS_32f-to22k-SSdp.wav -b 32f -r 22050 --noPeakChunk --doubleprecision
#$resampler_path --maxStages 1 -i $input_path/sweep-176400hz-0-22050hz-20s-D64-2.8mhz.dsf -o $output_path/96khz_sweep-3dBFS_32f-to32k-SSdp.wav -b 32f -r 32000 --noPeakChunk --doubleprecision --mt
#$resampler_path --maxStages 1 -i $input_path/sweep-176400hz-0-22050hz-20s-D64-2.8mhz.dsf -o $output_path/96khz_sweep-3dBFS_32f-to32k-SSdp.wav -b 32f -r 32000 --noPeakChunk --doubleprecision --mt
#$resampler_path --maxStages 1 -i $input_path/sweep-176400hz-0-22050hz-20s-D64-2.8mhz.dsf -o $output_path/96khz_sweep-3dBFS_32f-to44k-SSdp.wav -b 32f -r 44100 --noPeakChunk --doubleprecision
#$resampler_path --maxStages 1 -i $input_path/sweep-176400hz-0-22050hz-20s-D64-2.8mhz.dsf -o $output_path/96khz_sweep-3dBFS_32f-to48k-SSdp.wav -b 32f -r 48000 --noPeakChunk --doubleprecision --mt
#$resampler_path --maxStages 1 -i $input_path/sweep-176400hz-0-22050hz-20s-D64-2.8mhz.dsf -o $output_path/96khz_sweep-3dBFS_32f-to88k-SSdp.wav -b 32f -r 88200 --noPeakChunk --doubleprecision
#$resampler_path --maxStages 1 -i $input_path/sweep-176400hz-0-22050hz-20s-D64-2.8mhz.dsf -o $output_path/96khz_sweep-3dBFS_32f-to96k-SSdp.wav -b 32f -r 96000 --noPeakChunk --doubleprecision --mt
#$resampler_path --maxStages 1 -i $input_path/sweep-176400hz-0-22050hz-20s-D64-2.8mhz.dsf -o $output_path/96khz_sweep-3dBFS_32f-to176k-SSdp.wav -b 32f -r 176400 --noPeakChunk --doubleprecision
#$resampler_path --maxStages 1 -i $input_path/sweep-176400hz-0-22050hz-20s-D64-2.8mhz.dsf -o $output_path/96khz_sweep-3dBFS_32f-to192k-SSdp.wav -b 32f -r 192000 --noPeakChunk --doubleprecision --mt



# make spectrogr--ms
./make-spectrograms.sh