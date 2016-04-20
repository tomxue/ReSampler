// ReSampler.cpp : Audio Sample Rate Converter by Judd Niemann

#include <iostream>
#include <ostream>
#include <cassert>
#include <vector>
#include <algorithm>
#include <memory>
#include <iomanip>

#define _USE_MATH_DEFINES

#include <math.h>
#include "ReSampler.h"
#include "FIRFilter.h"
#include "ditherer.h"

////////////////////////////////////////////////////////////////////////////////////////
// This program uses the libsndfile Library, 
// available at http://www.mega-nerd.com/libsndfile/
//
// (copy of entire package included in $(ProjectDir)\libsbdfile)
//

#include "sndfile.hh"

//                                                                                    //
////////////////////////////////////////////////////////////////////////////////////////
 
int main(int argc, char * argv[])
{	
	//dumpKaiserWindow(29, 14);
	std::string sourceFilename("");
	std::string destFilename("");
	std::string outBitFormat("");
	int outFileFormat = 0;
	unsigned int OutputSampleRate;
	double NormalizeAmount = 1.0;
	double DitherAmount = 1.0;

	// parse core parameters:
	getCmdlineParam(argv, argv + argc, "-i", sourceFilename);
	getCmdlineParam(argv, argv + argc, "-o", destFilename);
	getCmdlineParam(argv, argv + argc, "-r", OutputSampleRate);
	getCmdlineParam(argv, argv + argc, "-b", outBitFormat);

	// parse help switch:
	if (findCmdlineOption(argv, argv + argc, "--help") || findCmdlineOption(argv, argv + argc, "-h")) {
		std::cout << strUsage << std::endl;
		std::cout << "Additional options:\n\n" << strExtraOptions << std::endl;
		exit(EXIT_SUCCESS);
	}
	
	// parse double precision switch:
	bool bUseDoublePrecision = findCmdlineOption(argv, argv + argc, "--doubleprecision");

	// parse bit format (sub format) parameter
	bool bListSubFormats = findCmdlineOption(argv, argv + argc, "--listsubformats");
	if (bListSubFormats) {
		std::string filetype;
		getCmdlineParam(argv, argv + argc, "--listsubformats", filetype);
		listSubFormats(filetype);
		exit(EXIT_SUCCESS);
	}
	
	// parse normalize option and parameter:
	bool bNormalize = findCmdlineOption(argv, argv + argc, "-n");
	if (bNormalize) {
		getCmdlineParam(argv, argv + argc, "-n", NormalizeAmount);
		if (NormalizeAmount <= 0.0)
			NormalizeAmount = 1.0;
		if (NormalizeAmount > 1.0)
			std::cout << "\nWarning: Normalization factor greater than 1.0 - THIS WILL CAUSE CLIPPING !!\n" << std::endl;
	}

	// parse dither option and parameter:
	bool bDither = findCmdlineOption(argv, argv + argc, "--dither");
	if (bDither) {
		getCmdlineParam(argv, argv + argc, "--dither", DitherAmount);
		if (DitherAmount <= 0.0)
			DitherAmount = 1.0;
	}
	
	bool bBadParams = false;
	if (destFilename.empty()) {
		if (sourceFilename.empty()) {
			std::cout << "Error: Input filename not specified" << std::endl;
			bBadParams = true;
		}
		else {
			std::cout << "Output filename not specified" << std::endl;
			destFilename = sourceFilename;
			if (destFilename.find(".") != std::string::npos) {
				auto dot = destFilename.find_last_of(".");
				destFilename.insert(dot, "(converted)");
			}
			else {
				destFilename.append("(converted)");
			}
			std::cout << "defaulting to: " << destFilename << "\n" << std::endl;
		}
	}
	else if (destFilename == sourceFilename) {
		std::cout << "\nError: Input and Output filenames cannot be the same" << std::endl;
		bBadParams = true;
	}

	if (OutputSampleRate == 0) {
		std::cout << "Error: Target sample rate not specified" << std::endl;
		bBadParams = true;
	}

	if (bBadParams) {
		std::cout << strUsage << std::endl;
		exit(EXIT_FAILURE);
	}

#ifdef _M_X64
	std::cout << "64-bit version\n" << std::endl;
#else
	std::cout << "32-bit version";
#if defined(USE_SSE2)
	std::cout << ", SSE2 build ... ";

	// Verify processor capabilities:

#if defined (_MSC_VER) || defined (__INTEL_COMPILER)
	bool bSSE2ok = false;
	int CPUInfo[4] = { 0,0,0,0 };
	__cpuid(CPUInfo, 0);
	if (CPUInfo[0] != 0) {
		__cpuid(CPUInfo, 1);
		if (CPUInfo[3] & (1 << 26))
			bSSE2ok = true;
	}
	if (bSSE2ok)
		std::cout << "CPU supports SSE2 (ok)";
	else {
		std::cout << "Your CPU doesn't support SSE2 - please try a non-SSE2 build on this machine" << std::endl;
		exit(EXIT_FAILURE);
	}
#endif // defined (_MSC_VER) || defined (__INTEL_COMPILER)
#endif // defined(USE_SSE2)
	std::cout << "\n" << std::endl;
#endif 

	std::cout << "Input file: " << sourceFilename << std::endl;
	std::cout << "Output file: " << destFilename << std::endl;
	
	double Limit = bNormalize ? NormalizeAmount : 1.0;

	// Isolate the file extensions
	std::string inFileExt("");
	std::string outFileExt("");
	
	if (sourceFilename.find_last_of(".") != std::string::npos)
		inFileExt = sourceFilename.substr(sourceFilename.find_last_of(".") + 1);
	
	if (destFilename.find_last_of(".") != std::string::npos)
		outFileExt = destFilename.substr(destFilename.find_last_of(".") + 1);

	if (!outBitFormat.empty()) { // new output bit format
		if (outFileFormat = determineOutputFormat(outFileExt, outBitFormat))
			std::cout << "Changing output bit format to " << outBitFormat << std::endl;
		else { // user-supplied bit format not valid; try choosing appropriate format
			determineBestBitFormat(outBitFormat, sourceFilename, destFilename); 
			if (outFileFormat = determineOutputFormat(outFileExt, outBitFormat))
				std::cout << "Changing output bit format to " << outBitFormat << std::endl;
			else {
				std::cout << "Warning: NOT Changing output file bit format !" << std::endl;
				outFileFormat = 0; // back where it started
			}
		}
	}	

	if (outFileExt != inFileExt)
	{ // file extensions differ, determine new output format: 

		if (outBitFormat.empty()) { // user changed file extension only. Attempt to choose appropriate output sub format
			std::cout << "Output Bit Format not specified" << std::endl;
			determineBestBitFormat(outBitFormat, sourceFilename, destFilename);
		}
		
		if(outFileFormat = determineOutputFormat(outFileExt, outBitFormat))
			std::cout << "Changing output file format to " << outFileExt << std::endl;
		else { // cannot determine subformat of output file
			std::cout << "Warning: NOT Changing output file format ! (extension different, but format will remain the same)" << std::endl;
		}
	}
		
	if (bUseDoublePrecision) {
		std::cout << "\nUsing double precision for calculations.\n" << std::endl;
		return (Convert<double>(sourceFilename, destFilename, OutputSampleRate, Limit, bNormalize, outFileFormat, bDither, DitherAmount)) ? EXIT_SUCCESS : EXIT_FAILURE;
	}
	else {
		return (Convert<float>(sourceFilename, destFilename, OutputSampleRate, Limit, bNormalize, outFileFormat, bDither, DitherAmount)) ? EXIT_SUCCESS : EXIT_FAILURE;
	}
}

// determineBestBitFormat() : determines the most appropriate bit format for the output file, through the following process:
// 1. Try to use infile's format and if that isn't valid for outfile, then:
// 2. use the default subformat for outfile.
// store best bit format as a string in BitFormat

bool determineBestBitFormat(std::string& BitFormat, const std::string& inFilename, const std::string& outFilename)
{
	// Inspect input file for format:
	SndfileHandle infile(inFilename, SFM_READ);
	int inFileFormat = infile.format();

	if (int e = infile.error()) {
		std::cout << "Couldn't Open Input File (" << sf_error_number(e) << ")" << std::endl;
		return false;
	}

	// get BitFormat of inFile as a string:
	for (auto& subformat : subFormats) {
		if (subformat.second == (inFileFormat & SF_FORMAT_SUBMASK)) {
			BitFormat = subformat.first;
			break;
		}
	}
	
	// get file extensions:
	std::string inFileExt("");
	if (inFilename.find_last_of(".") != std::string::npos)
		inFileExt = inFilename.substr(inFilename.find_last_of(".") + 1);

	std::string outFileExt("");
	if (outFilename.find_last_of(".") != std::string::npos)
		outFileExt = outFilename.substr(outFilename.find_last_of(".") + 1);

	// get total number of major formats:
	SF_FORMAT_INFO	formatinfo;
	int format, major_count;
	memset(&formatinfo, 0, sizeof(formatinfo));
	sf_command(NULL, SFC_GET_FORMAT_MAJOR_COUNT, &major_count, sizeof(int));
	
	// determine if inFile's subformat is valid for outFile
	for (int m = 0; m < major_count; m++)
	{
		formatinfo.format = m;
		sf_command(NULL, SFC_GET_FORMAT_MAJOR, &formatinfo, sizeof(formatinfo));	
		
		if (stricmp(formatinfo.extension, outFileExt.c_str()) == 0) {
			format = formatinfo.format | (inFileFormat & SF_FORMAT_SUBMASK); // combine outfile's major format with infile's subformat 
		
		   // Check if format / subformat combination is valid:
			SF_INFO sfinfo;
			memset(&sfinfo, 0, sizeof(sfinfo));
			sfinfo.channels = 1;
			sfinfo.format = format;
			if (!sf_format_check(&sfinfo))  { // not valid; use default format
				std::cout << "Output file format " << outFileExt << " and subformat " << BitFormat << " combination not valid ... ";
				BitFormat.clear();
				BitFormat = defaultSubFormats.find(outFileExt)->second;
				std::cout << "defaulting to " << BitFormat << std::endl;
				break;
			}	
		}
	}
	
	return true;
}

// determineOutputFormat() : returns an integer representing the output format, which libsndfile understands:
int determineOutputFormat(const std::string& outFileExt, const std::string& bitFormat)
{
	SF_FORMAT_INFO	info;
	int format = 0;
	int major_count;
	memset(&info, 0, sizeof(info));
	sf_command(NULL, SFC_GET_FORMAT_MAJOR_COUNT, &major_count, sizeof(int));
	bool bFileExtFound = false;

	// Loop through all major formats to find match for outFileExt:
	for (int m = 0; m < major_count; ++m) { 
		info.format = m;
		sf_command(NULL, SFC_GET_FORMAT_MAJOR, &info, sizeof(info));
		if (strcmpi(info.extension, outFileExt.c_str())==0) {	
			bFileExtFound = true;
			break;
		}
	}

	if (bFileExtFound) {
		// Check if subformat is recognized:
		auto sf = subFormats.find(bitFormat);
		if (sf != subFormats.end())
			format = info.format | sf->second;
		else
			std::cout << "Warning: bit format " << bitFormat << " not recognised !" << std::endl;
	}

	// Special cases:
	if (bitFormat == "8") {
		// user specified 8-bit. Determine whether it must be unsigned or signed, based on major type:
		// These formats always use unsigned 8-bit when they use 8-bit: mat rf64 voc w64 wav
		
		if ((outFileExt == "mat") || (outFileExt == "rf64") || (outFileExt == "voc") || (outFileExt == "w64") || (outFileExt == "wav"))
			format = info.format | SF_FORMAT_PCM_U8;
		else
			format = info.format | SF_FORMAT_PCM_S8;
	}

	return format;
}

// listSubFormats() - lists all valid subformats for a given file extension (without "." or "*."):
void listSubFormats(const std::string& f)
{
	SF_FORMAT_INFO	info;
	int format = 0;
	int major_count;
	memset(&info, 0, sizeof(info));
	sf_command(NULL, SFC_GET_FORMAT_MAJOR_COUNT, &major_count, sizeof(int));
	bool bFileExtFound = false;

	// Loop through all major formats to find match for outFileExt:
	for (int m = 0; m < major_count; ++m) {
		info.format = m;
		sf_command(NULL, SFC_GET_FORMAT_MAJOR, &info, sizeof(info));
		if (strcmpi(info.extension,f.c_str()) == 0) {
			bFileExtFound = true;
			break;
		}
	}
	if (bFileExtFound) {
		SF_INFO sfinfo;
		memset(&sfinfo, 0, sizeof(sfinfo));
		sfinfo.channels = 1;

		// loop through all subformats and find which ones are valid for file type:
		for (auto& subformat : subFormats) {
			sfinfo.format = (info.format & SF_FORMAT_TYPEMASK) | subformat.second;
			if (sf_format_check(&sfinfo))
				std::cout << subformat.first << std::endl;
		}
	}
	else {
		std::cout << "File extension " << f << " unknown" << std::endl;
	}
}

template<typename FloatType> 
bool Convert(const std::string& InputFilename, const std::string& OutputFilename, unsigned int OutputSampleRate, FloatType Limit, bool Normalize, int OutputFormat, bool bDither, FloatType DitherAmount)
{
	// Open input file:
	SndfileHandle infile(InputFilename, SFM_READ);
	
	if (int e=infile.error()) {
		std::cout << "Error: Couldn't Open Input File (" << sf_error_number(e) << ")" << std::endl;
		return false;
	}

	// Calculate conversion parameters, and print a summary for user:
	unsigned int nChannels = infile.channels();
	unsigned int InputSampleRate = infile.samplerate();
	
	int InputFileFormat = infile.format();
	bool bFloat = false;
	bool bDouble = false;

	// detect if input format is a floating-point format:
	switch (InputFileFormat & SF_FORMAT_SUBMASK) {
	case SF_FORMAT_FLOAT:
		bFloat = true;
		break;
	case SF_FORMAT_DOUBLE:
		bDouble = true;
		break;
	}

	size_t BufferSize = (BUFFERSIZE / nChannels) * nChannels; // round down to integer multiple of nChannels (file may have odd number of channels!)
	FloatType inbuffer[BUFFERSIZE];
	FloatType outbuffer[BUFFERSIZE];

	std::cout << "source file channels: " << nChannels << std::endl;
	std::cout << "input sample rate: " << InputSampleRate << "\noutput sample rate: " << OutputSampleRate << std::endl;
	
	for (auto& subformat : subFormats) { // scan subformats for a match:
		if (subformat.second == (InputFileFormat & SF_FORMAT_SUBMASK)) {
			std::cout << "input bit format: " << subformat.first;
			break;
		}
	}
	
	if (bFloat)
		std::cout << " (float)" << std::endl;
	if (bDouble)
		std::cout << " (double precision)" << std::endl;
	
	Fraction F = GetSimplifiedFraction(InputSampleRate, OutputSampleRate);
	FloatType ResamplingFactor = static_cast<FloatType>(OutputSampleRate) / InputSampleRate;
	std::cout << "\nConversion ratio: " << ResamplingFactor
		<< " (" << F.numerator << ":" << F.denominator << ") \n" << std::endl;
	
	FloatType PeakInputSample = 0;
	sf_count_t SamplesRead = 0i64;
	std::cout << "Scanning input file for peaks ..."; // to-do: can we read the PEAK chunk in floating-point files ?
	sf_count_t count;
	do {
		count = infile.read(inbuffer, BufferSize);
		SamplesRead += count;
		for (unsigned int s = 0; s < count; ++s) { // read all samples, without caring which channel they belong to
			PeakInputSample = max(PeakInputSample, abs(inbuffer[s]));
		}
	} while (count > 0);
	std::cout << "Done\n";
	std::cout << "Peak input sample: " << std::fixed << PeakInputSample << " (" << 20 * log10(PeakInputSample) << " dBFS)" << std::endl;
	infile.seek(0i64, SEEK_SET);

	// Calculate filter parameters:
	int OverSampFreq = InputSampleRate * F.numerator; // eg 160 * 44100
	unsigned int minSampleRate = min(InputSampleRate, OutputSampleRate);
	int TransitionWidth = minSampleRate / 22; // reasonable estimate for allowing transition width to scale with sample rate
	//TransitionWidth = max(TransitionWidth, 1700); // put a limit on how narrow the transition can be (too narrow will cause issues with low target samplerates)
	int ft = min(InputSampleRate, OutputSampleRate) / 2.0 - TransitionWidth;

	// Make some filters. Huge Filters are used for complex ratios.
	// Medium filters used for simple ratios (ie 1 in numerator or denominator)

	FloatType HugeFilterTaps[FILTERSIZE_HUGE];
	int HugeFilterSize = FILTERSIZE_HUGE;
	makeLPF<FloatType>(HugeFilterTaps, HugeFilterSize, ft, OverSampFreq);
	//applyKaiserWindow<FloatType>(HugeFilterTaps, HugeFilterSize, calcKaiserBeta(130));
	applyKaiserWindow<FloatType>(HugeFilterTaps, HugeFilterSize, 14);

	FloatType MedFilterTaps[FILTERSIZE_MEDIUM];
	int MedFilterSize = FILTERSIZE_MEDIUM;
	makeLPF<FloatType>(MedFilterTaps, MedFilterSize, ft, OverSampFreq);
	// applyKaiserWindow<FloatType>(MedFilterTaps, MedFilterSize, calcKaiserBeta(200));
	applyKaiserWindow<FloatType>(MedFilterTaps, MedFilterSize, 20);

	// make a vector of huge filters (one filter for each channel):
	std::vector<FIRFilter<FloatType, FILTERSIZE_HUGE>> HugeFilters;

	// make a vector of medium filters (one filter for each channel):
	std::vector<FIRFilter<FloatType, FILTERSIZE_MEDIUM >> MedFilters;

	for (unsigned int n = 0; n < nChannels; n++) {
		HugeFilters.emplace_back(HugeFilterTaps);
		MedFilters.emplace_back(MedFilterTaps);
	}

	// if the OutputFormat is zero, it means "No change to file format"
	// if output file format has changed, use OutputFormat. Otherwise, use same format as infile: 
	int OutputFileFormat = OutputFormat ? OutputFormat : InputFileFormat;

	// if the minor (sub) format of OutputFileFormat is not set, attempt to use minor format of input file (as a last resort)
	if ((OutputFileFormat & SF_FORMAT_SUBMASK) == 0) {
		OutputFileFormat |= (InputFileFormat & SF_FORMAT_SUBMASK); // may not be valid subformat for new file format. 
	}

	// determine number of bits in output format, for Ditherers:
	int signalBits;
	switch (OutputFileFormat & SF_FORMAT_SUBMASK) {
	case SF_FORMAT_PCM_24:
		signalBits = 24;
		break;
	case SF_FORMAT_PCM_S8:
	case SF_FORMAT_PCM_U8:
		signalBits = 8;
		break;
	default:
		signalBits = 16;
	}
	if (bDither)
		std::cout << "Generating dither for " << signalBits << "-bit output format" << std::endl;

	// make a vector of ditherers (one ditherer for each channel):
	std::vector<Ditherer<FloatType>> Ditherers;
	for (unsigned int n = 0; n < nChannels; n++) {
		Ditherers.emplace_back(signalBits, DitherAmount);
	}
	
	// Estimate initial gain:
	FloatType Gain = Normalize ? F.numerator * (Limit/PeakInputSample) : F.numerator * Limit;

	// Conditionally guard against filter overshoot:
	const FloatType OvershootCompensationFactor = 0.992; // Scaling factor to allow for filter overshoot
	if ((PeakInputSample > OvershootCompensationFactor) && !Normalize)
		Gain *= OvershootCompensationFactor;

	if (bDither) { // allow headroom for dithering:
		FloatType DitherCompensation =
			(pow(2, signalBits - 1) - pow(2, DitherAmount - 1)) / pow(2, signalBits - 1); // eg 32767/32768 = 0.999969 (-0.00027 dB)
			Gain *= DitherCompensation;
	}
	
	FloatType PeakOutputSample;
	SndfileHandle* pOutFile;
	//
	START_TIMER();

	do { // clipping detection loop (repeat if clipping detected)

		try { // Open output file:
			pOutFile = new SndfileHandle(OutputFilename, SFM_WRITE, OutputFileFormat, nChannels, OutputSampleRate); // needs to be dynamically allocated, 
			// because only way to close file is to go out of scope ... and we may need to overwrite file on 2nd pass

			if (int e = pOutFile->error()) {
				std::cout << "Error: Couldn't Open Output File (" << sf_error_number(e) << ")" << std::endl;
				return false;
			}
		}

		catch (std::bad_alloc& b) {
			std::cout << "Error: Couldn't Open Output File (memory allocation problem)" << std::endl;
			return false;
		}

		std::cout << "Converting ...";
		unsigned int DecimationIndex = 0;
		unsigned int OutBufferIndex = 0;
		PeakOutputSample = 0;
		
		if (F.numerator == 1) { // Decimate Only
			do { // Read and process blocks of samples until the end of file is reached
				count = infile.read(inbuffer, BufferSize);

				for (unsigned int s = 0; s < count; s += nChannels) {

					for (int Channel = 0; Channel < nChannels; Channel++)
						MedFilters[Channel].put(inbuffer[s + Channel]); // inject a source sample

					if (DecimationIndex == 0) { // Decimate
						for (int Channel = 0; Channel < nChannels; Channel++) {

							FloatType OutputSample = bDither ? 
								Ditherers[Channel].Dither(Gain * MedFilters[Channel].get()) : 
								Gain * MedFilters[Channel].get();

							outbuffer[OutBufferIndex + Channel] = OutputSample;
							PeakOutputSample = max(abs(PeakOutputSample), abs(OutputSample));
						}

						OutBufferIndex += nChannels;
						if (OutBufferIndex == BufferSize) {
							OutBufferIndex = 0;
							pOutFile->write(outbuffer, BufferSize);
						}
					}

					DecimationIndex++;
					if (DecimationIndex == F.denominator)
						DecimationIndex = 0;
				} // ends loop over s

			} while (count > 0);
		} // ends Decimate Only

		else if (F.denominator == 1) { // Interpolate only
			do { // Read and process blocks of samples until the end of file is reached
				count = infile.read(inbuffer, BufferSize);
				for (unsigned int s = 0; s < count; s += nChannels) {
					for (int ii = 0; ii < F.numerator; ++ii) {
						for (int Channel = 0; Channel < nChannels; Channel++) {
							if (ii == 0)
								MedFilters[Channel].put(inbuffer[s + Channel]); // inject a source sample
							else
								MedFilters[Channel].putZero(); // inject a Zero
							
							FloatType OutputSample = bDither ?
								Ditherers[Channel].Dither(Gain * MedFilters[Channel].LazyGet(F.numerator)) : 
								Gain * MedFilters[Channel].LazyGet(F.numerator);
							
							outbuffer[OutBufferIndex + Channel] = OutputSample;
							PeakOutputSample = max(PeakOutputSample, abs(OutputSample));
						}

						OutBufferIndex += nChannels;
						if (OutBufferIndex == BufferSize) {
							OutBufferIndex = 0;
							pOutFile->write(outbuffer, BufferSize);
						}
					} // ends loop over ii
				} // ends loop over s
			} while (count > 0);
		} // ends Interpolate Only

		else { // Interpolate and Decimate
			do { // Read and process blocks of samples until the end of file is reached
				count = infile.read(inbuffer, BufferSize);
				for (unsigned int s = 0; s < count; s += nChannels) {
					for (int ii = 0; ii < F.numerator; ++ii) { // (ii stands for "interpolation index")
						// Interpolate:
						if (ii == 0) { // inject a source sample
							for (int Channel = 0; Channel < nChannels; Channel++) {
								HugeFilters[Channel].put(inbuffer[s + Channel]);
							}
						}
						else { // inject a zero
							for (int Channel = 0; Channel < nChannels; Channel++) {
								HugeFilters[Channel].putZero();
							}
						}

						if (DecimationIndex == 0) { // decimate
							for (int Channel = 0; Channel < nChannels; Channel++) {

								FloatType OutputSample = bDither ?
									Ditherers[Channel].Dither(Gain * HugeFilters[Channel].LazyGet(F.numerator)) :
									Gain * HugeFilters[Channel].LazyGet(F.numerator);

								outbuffer[OutBufferIndex + Channel] = OutputSample;
								PeakOutputSample = max(PeakOutputSample, abs(OutputSample));
							}

							OutBufferIndex += nChannels;
							if (OutBufferIndex == BufferSize) {
								OutBufferIndex = 0;
								pOutFile->write(outbuffer, BufferSize);
							}
						}

						DecimationIndex++;
						if (DecimationIndex == F.denominator)
							DecimationIndex = 0;

						// To-do: showProgress();
					} // ends loop over ii
				} // ends loop over s
			} while (count > 0);
		} // ends Interpolate and Decimate

		if (OutBufferIndex != 0) {
			pOutFile->write(outbuffer, OutBufferIndex); // finish writing whatever remains in the buffer 
		}

		std::cout << "Done" << std::endl;
		std::cout << "\nPeak output sample: " << PeakOutputSample << " (" << 20 * log10(PeakOutputSample) << " dBFS)" << std::endl;
		
		delete pOutFile; // Close output file

		// To-do: Confirm assumption upon which the following statement is built:
		if (bDouble || bFloat)
			break; // Clipping is not a concern with Floating-Point formats. 

		// Test for clipping:
		if (PeakOutputSample > Limit) {
			FloatType GainAdjustment;
			
			if (bDither) {
				FloatType DitherLevel = pow(2, DitherAmount - 1) / pow(2, signalBits - 1);
				// take Ditherlevel out of calculations, as dither is NOT subject to gain control:
				GainAdjustment = (Limit-DitherLevel) / (PeakOutputSample-DitherLevel); 
			}

			else 
				GainAdjustment = Limit / PeakOutputSample;

			Gain *= GainAdjustment;
			std::cout << "\nClipping detected !" << std::endl;
			std::cout << "Re-doing with " << 20 * log10(GainAdjustment) << " dB gain adjustment" << std::endl;
			infile.seek(0i64,SEEK_SET);
		}
		
	} while (PeakOutputSample > Limit);

	STOP_TIMER();
	return true;
}

int gcd(int a, int b) {
	if (a<0) a = -a;
	if (b<0) b = -b;
	while (b != 0) {
		a %= b;
		if (a == 0) return b;
		b %= a;
	}
	return a;
}

Fraction GetSimplifiedFraction(int InputSampleRate, int OutputSampleRate)			// eg 44100, 48000
{
	Fraction f;
	f.numerator = (OutputSampleRate / gcd(InputSampleRate, OutputSampleRate));		// L (eg 160)
	f.denominator = (InputSampleRate / gcd(InputSampleRate, OutputSampleRate));		// M (eg 147)
	return f;
}

void getCmdlineParam(char** begin, char** end, const std::string& OptionName, std::string& Parameter) 
{
	Parameter = "";
	char** it = std::find(begin, end, OptionName);	
	if (it != end)	// found option
		if (++it != end) // found parameter after option
			Parameter = *it;
}

void getCmdlineParam(char** begin, char** end, const std::string& OptionName, unsigned int& nParameter)
{
	nParameter = 0;
	char** it = std::find(begin, end, OptionName);
	if (it != end)	// found option
		if (++it != end) // found parameter after option
			nParameter = atoi(*it);
}

void getCmdlineParam(char** begin, char** end, const std::string& OptionName, double& Parameter)
{
	Parameter = 0.0;
	char** it = std::find(begin, end, OptionName);
	if (it != end)	// found option
		if (++it != end) // found parameter after option
			Parameter = atof(*it);
}

bool findCmdlineOption(char** begin, char** end, const std::string& option) {
	return (std::find(begin, end, option) != end);
}