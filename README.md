## Synopsis
Resampler is a basic command-line audio sample rate conversion tool for Windows, which can convert *.wav file formats with variety of different bit-depths and audio channel configurations. 

Resampler is intended to produce outstanding quality sound files, keeping aliasing and other unwanted artifacts to a minimum.

## Description of code
 
Sample Rate Conversion is accomplished using the usual interpolation/decimation techniques. However, when using complex conversion ratios (such as 44.1k <--> 48k), a rather large FIR lowpass filter is used to ensure a clean conversion.

Resampler uses the C++ wrapper of the outstanding [libsndfile](http://www.mega-nerd.com/libsndfile/) library for sound file I/O operations. (I originally embarked upon writing my own sound-file I/O library, but quickly realised the enormity of such an undertaking, and subsequently decided it was a whole world-of-hurt which would be best avoided)

Resampler needs to be compiled on Visual C++ 2015, as it uses some C++11 features. (Porting to other environments is intended in the future)

## Motivation
This project arose out of: 

* my own experimentation with digital filters and some basic DSP concepts
* a requirement to have a simple command-line tool to be used in a script to convert a large collection of audio files
* a need for a *quality* SRC tool, as the quality of other offerings (both commercial and free) varies wildly from terrific to appalling

Future versions of this project are anticipated to receive a dramatic speed improvement through the use of SIMD vectorisation, multithreading, and more efficient FIR filter designs. 

## Usage

from the command line, the main options are as follows:

**resampler.exe -i inputfile [-o outputfile] -r samplerate [-n [normalization factor]]**
 
**samplerate** is the target sample rate in Hz.

**Normalization factor** is **> 0.0** and **<= 1.0**, with 1.0 producing the largest possible output level without clipping. Note: resampler will accept normalization values over 1.0, but this will certainly result in clipping, and is therefore only for experimental and testing purposes. Just using **-n** with no parameter is equivalent to **-n 1.0**

Additional options:

**[--doubleprecision]** will force resampler to use double-precision arithmetic for its *internal calculations* and doesn't have anything to do with the file formats, although if you are working with 64-bit double-precision files, it would make sense to use double precision for calculations used in processing.