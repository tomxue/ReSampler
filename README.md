## Synopsis
Resampler is a basic command-line audio sample rate conversion tool for Windows, which can convert audio file formats with a variety of different bit-depths and audio channel configurations. 

Resampler is intended to produce outstanding quality sound files, keeping aliasing and other unwanted artifacts to a minimum, as the following actual measurement graphs show:

![FIR Filter response - downsample 96k->44k](https://github.com/jniemann66/ReSampler/blob/master/147_160_FIR_Frequency-response-2016-03-30-KAISER_WINDOW.JPG)   
*Typical frequency response when downsampling from 96kHz to 44.1kHz*

![spectrogram: 0-48khz sweep ](https://github.com/jniemann66/ReSampler/blob/master/96khz_sweep-3dBFS_64f(to44k).png)  
*Spectrogram of 0-48kHz Sine Sweep @96kHz sample rate, after having been downsampled to 44kHz sample rate*

## Description of code
 
Sample Rate Conversion is accomplished using the usual interpolation/decimation techniques. However, when using complex conversion ratios (such as 44.1k <--> 48k), a rather large FIR lowpass filter is used to ensure a clean conversion.

Resampler uses the C++ wrapper of the outstanding [libsndfile](http://www.mega-nerd.com/libsndfile/) library for sound file I/O operations. (I originally embarked upon writing my own sound-file I/O library, but quickly realised the enormity of such an undertaking, and subsequently adopted libsndfile)

The FIR filter class written for this project uses SSE SIMD instructions when using single-precision floating-point to perform 4 multiply/accumulate operations simultaneously. This has been found to yield a speed improvement of approximately 3x. Some experimentation was also done with 2x double-precision SSE2 SIMD, but was found to be no faster than the basic scalar implementation, and has thus been commented-out (but not deleted) from the source.

Resampler was developed on Visual C++ 2015, as it uses some C++11 features. (Porting to other environments is intended in the future)

## Motivation
This project arose out of: 

* my own experimentation with digital filters and some basic DSP concepts
* a requirement to have a simple command-line tool to be used in a script to convert a large collection of audio files
* a need for a *quality* SRC tool, as the quality of other offerings (both commercial and free) varies wildly from terrific to appalling

Future versions of this project are anticipated to receive a dramatic speed improvement through the use of multithreading, and more efficient FIR filter designs. 

## Usage

from the command line, the main options are as follows:

**resampler.exe -i inputfile [-o outputfile] -r samplerate [-b bitformat] [-n [normalization factor]]**
 
**samplerate** is the target sample rate in Hz.

**bitformat** is the bit representation (sub format) of the data in the output file. If this option is omitted, resampler will try to deduce the intended bit format automatically. Not all bit formats are valid for a given output file type. For more details, refer to the [libsndfile](http://www.mega-nerd.com/libsndfile/) documentation. Here is a list of all subformats: 

    8			8-bit (signed or unsigned automatic, based on file type)
    s8			Signed 8 bit data
    16			Signed 16 bit data
    24			Signed 24 bit data
    32			Signed 32 bit data
    u8			Unsigned 8 bit data
    32f			32 bit float data
    64f			64 bit float data
    ulaw		U-Law encoded
    alaw		A-Law encoded
    ima-adpcm	IMA ADPCM
    ms-adpcm	Microsoft ADPCM
    gsm610		GSM 6.10 encoding
    vox-adpcm	OKI Dialogix ADPCM
    g721-32		32kbs G721 ADPCM encoding
    g723-24		24kbs G723 ADPCM encoding
    g723-40		40kbs G723 ADPCM encoding
    dwvw12		12 bit Delta Width Variable Word encoding
    dwvw16		16 bit Delta Width Variable Word encoding
    dwvw24		24 bit Delta Width Variable Word encoding
    dwvwn		N bit Delta Width Variable Word encoding
    dpcm8		8 bit differential PCM (XI only)
    dpcm16		16 bit differential PCM (XI only)
    vorbis		Xiph Vorbis encoding
    alac16		Apple Lossless Audio Codec (16 bit)
    alac20		Apple Lossless Audio Codec (20 bit)
    alac24		Apple Lossless Audio Codec (24 bit)
    alac32		Apple Lossless Audio Codec (32 bit)

*Note: the --listsubformats option will cause the program to display the valid formats for a given file-type*

**Normalization factor** is **> 0.0** and **<= 1.0**, with 1.0 producing the largest possible output level without clipping. Note: resampler will accept normalization values over 1.0, but this will certainly result in clipping, and is therefore only for experimental and testing purposes. Just using **-n** with no parameter is equivalent to **-n 1.0**

### Additional options: ###

**--doubleprecision** will force resampler to use double-precision arithmetic for its *internal calculations* and doesn't have anything to do with the file formats, although if you are working with 64-bit double-precision files, it would make sense to use double precision for calculations used in processing.

**--dither [amount]** adds **+/-amount** *bits* of dither the output file. Dithering deliberately adds a small amount of a particular type of noise (triangular pdf with noise-shaping) prior to quantization to the output file. The goal of dithering is to reduce distortion, and allow extremely quiet passages to be preserved when they would otherwise be below the threshold of the target bit depth. Usually, it only makes sense to add dither when you are converting to a lower bit depth, for example:
 
- floating-point -> 24bit
- 24bit -> 16bit
- 16bit -> 8bit

The *amount* parameter represents the number of *bits* of dither to add. The actual *level* of dithering added is equal to **+/- 2^(amount-1)** *steps* (in other words, 2\* *amount* bits peak-to-peak). The default is for *amount* is 1.0, and it doesn't need to be an integer. Values in the range 1-6 are sensible for most situations. The noise-shaping curve becomes more pronounced at higher dithering amounts.  

The effect of dithering is most noticable during extremely quiet passages (typically, in fade-outs) of the audio. If you can hear modulation effects, or "tearing" in the quietest passages of your output file, then a greater amount of dither may need to be applied. (note: in many cases, these passages are so quiet, you will need to normalize them just to hear them).

**--autoblank** when specified in conjuction with **--dither** , mute the dithering after 30,000 consecutive input samples of *silence* (< -193dB is considered silence). Dithering is re-enabled immediately upon a non-zero input sample being detected.

**--listsubformats <filetype\>** will list all valid subformats for a given *filetype*

**--version** will display the version number of the program

**--minphase** use a minimum-phase FIR filter, instead of Linear-Phase

**--flacCompression <compressionlevel\>** sets the compression level for flac output files (between 0 and 8)

**--vorbisQuality <quality\>** sets the quality level for ogg vorbis output files (between -1 and 10)

**--noClippingProtection** diables clipping protection (clipping protection is normally active by default)

 

## Supported Formats

Resampler can handle any of the file formats libsndfile can handle.  
Thus, the following file extensions are supported:  
![Supported Formats](https://github.com/jniemann66/ReSampler/blob/master/supported_formats.png)

For more information, please refer to the [libsndfile documentation](http://www.mega-nerd.com/libsndfile/)

## Additional Information

####Clipping Protection

Resampler employs a multiple-pass approach with regards to clipping detection. If clipping is detected (ie normalized signal level exceeded +/- 1.0) during processing, it will re-do the conversion with the overall gain adjusted appropriately to avoid clipping. 

####Conversion to same sampling rate

When the target sampling rate is the same as the input file (ie 1:1 ratio), sample-rate conversion is not actually performed. However, bit-depth / file format conversion and other features such as dithering and normalization are performed when requested.