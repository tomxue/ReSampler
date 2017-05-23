#!/usr/bin/env bash

resampler_path=../x64/Release/ReSampler.exe
input_path=./inputs
output_path=./outputs

# clear old outputs:
rm $output_path/*

# 32-bit float sweeps:
$resampler_path -i $input_path/96khz_sweep-3dBFS_32f.wav -o $output_path/96khz_sweep-3dBFS_32f-to22k.wav -r 22050 --noPeakChunk
$resampler_path -i $input_path/96khz_sweep-3dBFS_32f.wav -o $output_path/96khz_sweep-3dBFS_32f-to32k.wav -r 32000 --noPeakChunk --mt
$resampler_path -i $input_path/96khz_sweep-3dBFS_32f.wav -o $output_path/96khz_sweep-3dBFS_32f-to44k.wav -r 44100 --noPeakChunk
$resampler_path -i $input_path/96khz_sweep-3dBFS_32f.wav -o $output_path/96khz_sweep-3dBFS_32f-to48k.wav -r 48000 --noPeakChunk --mt
$resampler_path -i $input_path/96khz_sweep-3dBFS_32f.wav -o $output_path/96khz_sweep-3dBFS_32f-to88k.wav -r 88200 --noPeakChunk
$resampler_path -i $input_path/96khz_sweep-3dBFS_32f.wav -o $output_path/96khz_sweep-3dBFS_32f-to96k.wav -r 96000 --noPeakChunk --mt
$resampler_path -i $input_path/96khz_sweep-3dBFS_32f.wav -o $output_path/96khz_sweep-3dBFS_32f-to176k.wav -r 176400 --noPeakChunk
$resampler_path -i $input_path/96khz_sweep-3dBFS_32f.wav -o $output_path/96khz_sweep-3dBFS_32f-to192k.wav -r 192000 --noPeakChunk --mt
# note: with floating-point formats, the PEAK chunk has a unique timestamp in it, which throws the MD5 checksums, hence --noPeakChunk 

# 32-bit float, double-precision sweeps:
$resampler_path -i $input_path/96khz_sweep-3dBFS_32f.wav -o $output_path/96khz_sweep-3dBFS_32f-to22k-dp.wav -r 22050 --noPeakChunk --doubleprecision
$resampler_path -i $input_path/96khz_sweep-3dBFS_32f.wav -o $output_path/96khz_sweep-3dBFS_32f-to32k-dp.wav -r 32000 --noPeakChunk --doubleprecision --mt
$resampler_path -i $input_path/96khz_sweep-3dBFS_32f.wav -o $output_path/96khz_sweep-3dBFS_32f-to44k-dp.wav -r 44100 --noPeakChunk --doubleprecision
$resampler_path -i $input_path/96khz_sweep-3dBFS_32f.wav -o $output_path/96khz_sweep-3dBFS_32f-to48k-dp.wav -r 48000 --noPeakChunk --doubleprecision --mt
$resampler_path -i $input_path/96khz_sweep-3dBFS_32f.wav -o $output_path/96khz_sweep-3dBFS_32f-to88k-dp.wav -r 88200 --noPeakChunk --doubleprecision
$resampler_path -i $input_path/96khz_sweep-3dBFS_32f.wav -o $output_path/96khz_sweep-3dBFS_32f-to96k-dp.wav -r 96000 --noPeakChunk --doubleprecision --mt
$resampler_path -i $input_path/96khz_sweep-3dBFS_32f.wav -o $output_path/96khz_sweep-3dBFS_32f-to176k-dp.wav -r 176400 --noPeakChunk --doubleprecision
$resampler_path -i $input_path/96khz_sweep-3dBFS_32f.wav -o $output_path/96khz_sweep-3dBFS_32f-to192k-dp.wav -r 192000 --noPeakChunk --doubleprecision --mt

# 32-bit float, double-precision sweeps, steep LPF:
$resampler_path -i $input_path/96khz_sweep-3dBFS_32f.wav -o $output_path/96khz_sweep-3dBFS_32f-to22k-dp-steep.wav -r 22050 --noPeakChunk --doubleprecision --steepLPF
$resampler_path -i $input_path/96khz_sweep-3dBFS_32f.wav -o $output_path/96khz_sweep-3dBFS_32f-to32k-dp-steep.wav -r 32000 --noPeakChunk --doubleprecision --steepLPF --mt
$resampler_path -i $input_path/96khz_sweep-3dBFS_32f.wav -o $output_path/96khz_sweep-3dBFS_32f-to44k-dp-steep.wav -r 44100 --noPeakChunk --doubleprecision --steepLPF
$resampler_path -i $input_path/96khz_sweep-3dBFS_32f.wav -o $output_path/96khz_sweep-3dBFS_32f-to48k-dp-steep.wav -r 48000 --noPeakChunk --doubleprecision --steepLPF --mt
$resampler_path -i $input_path/96khz_sweep-3dBFS_32f.wav -o $output_path/96khz_sweep-3dBFS_32f-to88k-dp-steep.wav -r 88200 --noPeakChunk --doubleprecision --steepLPF
$resampler_path -i $input_path/96khz_sweep-3dBFS_32f.wav -o $output_path/96khz_sweep-3dBFS_32f-to96k-dp-steep.wav -r 96000 --noPeakChunk --doubleprecision --steepLPF --mt
$resampler_path -i $input_path/96khz_sweep-3dBFS_32f.wav -o $output_path/96khz_sweep-3dBFS_32f-to176k-dp-steep.wav -r 176400 --noPeakChunk --doubleprecision --steepLPF
$resampler_path -i $input_path/96khz_sweep-3dBFS_32f.wav -o $output_path/96khz_sweep-3dBFS_32f-to192k-dp-steep.wav -r 192000 --noPeakChunk --doubleprecision --steepLPF --mt

# 16-bit sweeps, relaxed LPF
$resampler_path -i $input_path/96khz_sweep-3dBFS_32f.wav -o $output_path/96khz_sweep-3dBFS_16-to22k-relaxed.wav -r 22050 -b 16 --relaxedLPF
$resampler_path -i $input_path/96khz_sweep-3dBFS_32f.wav -o $output_path/96khz_sweep-3dBFS_16-to32k-relaxed.wav -r 32000 -b 16 --relaxedLPF --mt
$resampler_path -i $input_path/96khz_sweep-3dBFS_32f.wav -o $output_path/96khz_sweep-3dBFS_16-to44k-relaxed.wav -r 44100 -b 16 --relaxedLPF
$resampler_path -i $input_path/96khz_sweep-3dBFS_32f.wav -o $output_path/96khz_sweep-3dBFS_16-to48k-relaxed.wav -r 48000 -b 16 --relaxedLPF --mt
$resampler_path -i $input_path/96khz_sweep-3dBFS_32f.wav -o $output_path/96khz_sweep-3dBFS_16-to88k-relaxed.wav -r 88200 -b 16 --relaxedLPF
$resampler_path -i $input_path/96khz_sweep-3dBFS_32f.wav -o $output_path/96khz_sweep-3dBFS_16-to96k-relaxed.wav -r 96000 -b 16 --relaxedLPF --mt
$resampler_path -i $input_path/96khz_sweep-3dBFS_32f.wav -o $output_path/96khz_sweep-3dBFS_16-to176k-relaxed.wav -r 176400 -b 16 --relaxedLPF
$resampler_path -i $input_path/96khz_sweep-3dBFS_32f.wav -o $output_path/96khz_sweep-3dBFS_16-to192k-relaxed.wav -r 192000 -b 16 --relaxedLPF --mt

# check files against last results:
md5sum -c checklist.chk > results.txt

# generate new checksums:
md5sum $output_path/* > checklist.chk