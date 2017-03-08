## Synopsis
Resampler is a high-performance command-line audio sample rate conversion tool for Windows, which can convert audio file formats with a variety of different bit-depths and audio channel configurations. 

Resampler is intended to produce outstanding quality sound files, keeping aliasing and other unwanted artifacts to a minimum, as the following actual measurement graphs show:

![FIR Filter response - downsample 96k->44k](https://github.com/jniemann66/ReSampler/blob/master/147_160_FIR_Frequency-response-2016-03-30-KAISER_WINDOW.JPG)   
*Typical frequency response when downsampling from 96kHz to 44.1kHz*

![spectrogram: 0-48khz sweep ](https://github.com/jniemann66/ReSampler/blob/master/96khz_sweep-3dBFS_64f(to44k).png)  
*Spectrogram of 0-48kHz Sine Sweep @96kHz sample rate, after having been downsampled to 44kHz sample rate*

## Motivation
This project arose out of: 

* my own experimentation with digital filters and some basic DSP concepts
* a requirement to have a simple command-line tool to be used in a script to convert a large collection of audio files
* a need for a *quality* SRC tool, as the quality of other offerings (both commercial and free) varies wildly from terrific to appalling

Future versions of this project are anticipated to receive a dramatic speed improvement through the use of multithreading, and more efficient FIR filter designs. 

## Usage

from the command line, the main options are as follows:

**resampler.exe -i inputfile [-o outputfile] -r samplerate [-b bitformat] [-n [&lt;normalization factor&gt;]]**
 
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

**Normalization factor** is a value between **0.0** and **1.0**, with 1.0 (equivalent to 100 percent) producing the largest possible output level without clipping. Note: resampler will accept normalization values over 1.0, but this will certainly result in clipping, and is therefore only for experimental and testing purposes. Just using **-n** with no parameter is equivalent to **-n 1.0**

### Additional options: ###

**--mt** Multi-Threading (since v1.2.2). This will cause ReSampler to process each channel in a separate thread. 
On a multi-core system, this makes better use of available CPU resources and results in a significant speed improvement.  

**--doubleprecision** will force resampler to use double-precision (64-bit floating point) arithmetic for its *internal calculations* and doesn't have anything to do with the file formats, although if you are working with 64-bit double-precision files, it would make sense to use double precision for calculations used in processing.

**--dither [&lt;amount&gt;]** adds **+/-amount** *bits* of dither the output file. Dithering deliberately adds a small amount of a particular type of noise (triangular pdf with noise-shaping) prior to quantization to the output file. The goal of dithering is to reduce distortion, and allow extremely quiet passages to be preserved when they would otherwise be below the threshold of the target bit depth. Usually, it only makes sense to add dither when you are converting to a lower bit depth, for example:
 
- floating-point -> 24bit
- 24bit -> 16bit
- 16bit -> 8bit

The *amount* parameter represents the number of *bits* of dither to add. The actual *level* of dithering added is equal to **+/- 2^(amount-1)** *steps*. The default is for *amount* is 1.0, and it doesn't need to be an integer. Values in the range 1-6 are sensible for most situations. The noise-shaping curve becomes more pronounced at higher dithering amounts.  

The effect of dithering is most noticable during extremely quiet passages (typically, in fade-outs) of the audio. If you can hear modulation effects, or "tearing" in the quietest passages of your output file, then a greater amount of dither may need to be applied. (note: in many cases, these passages are so quiet, you will need to normalize them just to hear them).

**--autoblank** when specified in conjuction with **--dither** causes dithering to switch-off after 30,000 consecutive input samples of *silence* (< -193dB is considered silence). Dithering is re-enabled immediately upon a non-zero input sample being detected.

**--seed &lt;n&gt;** (since v1.1.5) when specified in conjuction with **--dither** , causes the pseudo-random number generator used to generate dither noise to generate a specific sequence of noise associated with the number n.
Using the same value of n on subsequent conversions should reproduce precisely the same result. n is a signed integer in the range -2,147,483,648 through 2,147,483,647.  

**--flat-tpdf** (since v1.1.6) when specified in conjuction with **--dither** , causes the dithering to use flat tpdf noise with no noise-shaping.

**--listsubformats &lt;filetype&gt;** will list all valid subformats for a given *filetype*

**--version** will display the version number of the program

**--minphase** use a minimum-phase FIR filter, instead of Linear-Phase

**--relaxedLPF** (since v1.1.4) causes the lowpass filter to use a "late" cutoff frequency (95.45%), 
which will (theoretically) allow a small amount of aliasing, but at the same time, keep ringing to a minimum and maintain a good frequency response.

**--steepLPF** (since v1.2.0) causes the lowpass filter to use a steeper cutoff (half the standard transition width). It avoids aliasing, but may result in more ringing.

**--flacCompression  &lt;compressionlevel&gt;** sets the compression level for flac output files (between 0 and 8)

**--vorbisQuality &lt;quality&gt;** sets the quality level for ogg vorbis output files (between -1 and 10)

**--noClippingProtection** disables clipping protection (clipping protection is normally active by default)

**--rf64** (since v1.2.6) forces output .wav file to be in rf64 format. Has no effect if output file is not a .wav file.

## Supported Formats

Resampler can handle any of the file formats libsndfile can handle, plus a few extras (notably, the 1-bit DSD formats, dff and dsf - as of v1.2.0).  
Thus, the following file extensions are supported:  
![Supported Formats](https://github.com/jniemann66/ReSampler/blob/master/supported_formats.png)

For more information, please refer to the [libsndfile documentation](http://www.mega-nerd.com/libsndfile/)

## Additional Information

#### Clipping Protection

Resampler employs a multiple-pass approach with regards to clipping detection. If clipping is detected (ie normalized signal level exceeded +/- 1.0) during processing, it will re-do the conversion with the overall gain adjusted appropriately to avoid clipping. (This can be disabled with the **--noClippingProtection** option)

#### Conversion to same sampling rate

When the target sampling rate is the same as the input file (ie 1:1 ratio), sample-rate conversion is not actually performed. However, bit-depth / file format conversion and other features such as dithering and normalization are performed when requested.

#### Automatic promotion of **wav** files to **rf64** format

(since v1.2.3) For **wav** file output, if the data contained in the output file will exceed 4 Gigabytes after conversion, ReSampler will automatically "promote" the file format to **rf64** 

*Traditional wav files only have 32-bits to quantify the size of the data they contain, which means that the data size cannot exceed 4 Gigabytes (2^32 bytes). The [rf64 specification](https://tech.ebu.ch/docs/tech/tech3306-2009.pdf) was designed to extend the wav format to remove this limitation.*

Alternatively, other output formats suitable for *large* output files (exceeding 4GB) are: **w64, caf, au** or compressed formats such as **flac** or **oga**

Note: the **--rf64** option will force output .wav files to be in rf64 format.   

## Description of code
 
Sample Rate Conversion is accomplished using the usual interpolation/decimation techniques. However, when using complex conversion ratios (such as 44.1k <--> 48k), a rather large FIR lowpass filter is used to ensure a clean conversion.

Resampler uses the C++ wrapper of the outstanding [libsndfile](http://www.mega-nerd.com/libsndfile/) library for sound file I/O operations. (I originally embarked upon writing my own sound-file I/O library, but quickly realised the enormity of such an undertaking, and subsequently adopted libsndfile)

The FIR filter class written for this project uses SSE SIMD instructions when using single-precision floating-point to perform 4 multiply/accumulate operations simultaneously. This has been found to yield a speed improvement of approximately 3x. Some experimentation was also done with 2x double-precision SSE2 SIMD, but was found to be no faster than the basic scalar implementation, and has thus been commented-out (but not deleted) from the source.

(Additionally, some progress has been made recently with a build of the project using AVX instructions, on supported CPUs / OSes, to perform 8 single-precision or 4 double-precision multiply/accumulate operations at a time.) 

Resampler was developed on Visual C++ 2015, as it uses some C++11 features. (Porting to other environments is intended in the future).

#### explanation of source code files:

----------

**resampler.cpp** :	core of program

**resampler.h** : important data structures and function declarations

**FIRFilter.h** : FIR Filter DSP code

**FIRFilterAVX.h** : AVX-specific DSP code (conditional #include in AVX build)
 
**Biquad.h** : IIR Filter (used in dithering)

**Ditherer.h** : defines ditherer class, for adding dither

**dff.h** : module for reading dff files

**dsf.h** : module for reading dsf files

**alignedmalloc.h** : simple function for dynamically allocating aligned memory (AVX requires 32-byte alignment)

*(the class implementations are header-only)*

----------

## Description of Binaries included in distribution

**ReSampler/Release/ReSampler.exe** : 32-bit Windows with **SSE2** instruction set

**ReSampler/x64/Release/ReSampler.exe** : 64-bit Windows (*uses SSE2, but all 64-bit CPUs should have this*)

**ReSampler/x64/AVX_Release/ReSampler.exe** : 64-bit Windows with **AVX** instruction set (Requires Intel *Sandy Bridge* CPU or higher, AMD *Bulldozer* or higher. Requires Windows 7 SP1 or higher, Windows Server 2008 R2 SP1 or higher OS) 

----------

## Acknowledgements

The following libraries are used in ReSampler:

<table>
    <thead>
        <tr>
            <th>Name</th>
            <th>Description</th>
            <th>Reason For Inclusion</th>
            <th>Link</th>
            <th>License</th>
        </tr>
    </thead>
    <tbody>
        <tr>
            <td>libsndfile</td>
            <td>C library for reading and writing files containing sampled sound</td>
            <td>Core read/write capabilities for common file formats</td>
            <td><a href="http://www.mega-nerd.com/libsndfile/">libsndfile</a></td>
            <td>LGPL 2.1</td>
        </tr>
        <tr>
            <td>FFTW</td>
            <td>C subroutine library for computing the discrete Fourier transform (DFT)</td>
            <td>used for conversion of linear-phase impulses to minimum-phase</td>
            <td><a href="http://www.fftw.org/">fftw</a></td>
            <td>GPL (2?)</td>
        </tr>
        <tr>
            <td>CTPL</td>
            <td>Modern and efficient C++ Thread Pool Library</td>
            <td>Used for multi-threaded sample-rate conversion</td>
            <td><a href="https://github.com/vit-vit/CTPL">CTPL</a></td>
            <td>Apache 2.0</td>
        </tr>
    </tbody>
</table>