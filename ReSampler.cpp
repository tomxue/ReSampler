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
#include "biquad.h"

////////////////////////////////////////////////////////////////////////////////////////
// This program uses the following libraries:
// 1:
// libsndfile                                          
// available at http://www.mega-nerd.com/libsndfile/
//
// (copy of entire package included in $(ProjectDir)\libsbdfile)
// 
// 2:
// fftw
// http://www.fftw.org/
// 

#include "sndfile.hh"

//                                                                                    //
////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char * argv[])
{

	std::string sourceFilename("");
	std::string destFilename("");
	std::string outBitFormat("");
	int outFileFormat = 0;
	unsigned int OutputSampleRate = 44100;
	double NormalizeAmount = 1.0;
	double DitherAmount = 1.0;
	int flacCompressionLevel = 5;
	double vorbisQuality = 3;

	// parse core parameters:
	getCmdlineParam(argv, argv + argc, "-i", sourceFilename);
	getCmdlineParam(argv, argv + argc, "-o", destFilename);
	getCmdlineParam(argv, argv + argc, "-r", OutputSampleRate);
	getCmdlineParam(argv, argv + argc, "-b", outBitFormat);

	// parse version switch:
	if (findCmdlineOption(argv, argv + argc, "--version")) {
		std::cout << strVersion << std::endl;
		exit(EXIT_SUCCESS);
	}

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

	// parse auto-blanking option (for dithering):
	bool bAutoBlankingEnabled = findCmdlineOption(argv, argv + argc, "--autoblank");

	// parse minimum-phase option:
	bool bMinPhase = findCmdlineOption(argv, argv + argc, "--minphase");

	// parse flacCompression option and parameter:
	bool bSetFlacCompression = findCmdlineOption(argv, argv + argc, "--flacCompression");
	if (bSetFlacCompression) {
		getCmdlineParam(argv, argv + argc, "--flacCompression",flacCompressionLevel);
		if (flacCompressionLevel < 0)
			flacCompressionLevel = 0;
		if (flacCompressionLevel > 8)
			flacCompressionLevel = 8;
	}

	// parse vorbisQuality option and parameter:
	bool bSetVobisQuality = findCmdlineOption(argv, argv + argc, "--vorbisQuality");
	if (bSetVobisQuality) {
		getCmdlineParam(argv, argv + argc, "--vorbisQuality", vorbisQuality);
		if (vorbisQuality < -1)
			vorbisQuality = -1;
		if (vorbisQuality > 10)
			vorbisQuality = 10;
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

	std::cout << strVersion << " ";
#ifdef _M_X64
	std::cout << "64-bit version" << std::endl;
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

		if (outBitFormat.empty()) { // user changed file extension only. Attempt to choose appropriate output sub format:
			std::cout << "Output Bit Format not specified" << std::endl;
			determineBestBitFormat(outBitFormat, sourceFilename, destFilename);
		}

		if (outFileFormat = determineOutputFormat(outFileExt, outBitFormat))
			std::cout << "Changing output file format to " << outFileExt << std::endl;
		else { // cannot determine subformat of output file
			std::cout << "Warning: NOT Changing output file format ! (extension different, but format will remain the same)" << std::endl;
		}
	}
	if (bUseDoublePrecision) {
		std::cout << "Using double precision for calculations." << std::endl;
		using floattype = double;
		conversionInfo<double> ci;
		ci.InputFilename = sourceFilename;
		ci.OutputFilename = destFilename;
		ci.OutputSampleRate = OutputSampleRate;
		ci.Limit = Limit;
		ci.bNormalize = bNormalize;
		ci.OutputFormat = outFileFormat;
		ci.bDither = bDither;
		ci.DitherAmount = DitherAmount;
		ci.bAutoBlankingEnabled = bAutoBlankingEnabled;
		ci.bMinPhase = bMinPhase;
		ci.bSetFlacCompression = bSetFlacCompression;
		ci.flacCompressionLevel = flacCompressionLevel;
		ci.bSetVorbisQuality = bSetVobisQuality;
		ci.vorbisQuality = vorbisQuality;
		return Convert<double>(ci) ? EXIT_SUCCESS : EXIT_FAILURE;
	}

	else {
		conversionInfo<float> ci;
		ci.InputFilename = sourceFilename;
		ci.OutputFilename = destFilename;
		ci.OutputSampleRate = OutputSampleRate;
		ci.Limit = Limit;
		ci.bNormalize = bNormalize;
		ci.OutputFormat = outFileFormat;
		ci.bDither = bDither;
		ci.DitherAmount = DitherAmount;
		ci.bAutoBlankingEnabled = bAutoBlankingEnabled;
		ci.bMinPhase = bMinPhase;
		ci.bSetFlacCompression = bSetFlacCompression;
		ci.flacCompressionLevel = flacCompressionLevel;
		ci.bSetVorbisQuality = bSetVobisQuality;
		ci.vorbisQuality = vorbisQuality;
		return Convert<float>(ci) ? EXIT_SUCCESS : EXIT_FAILURE;
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
			if (!sf_format_check(&sfinfo)) { // not valid; use default format
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
		if (strcmpi(info.extension, outFileExt.c_str()) == 0) {
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
		if (strcmpi(info.extension, f.c_str()) == 0) {
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
bool Convert(const conversionInfo<FloatType>& ci)
{
	// Open input file:
	SndfileHandle infile(ci.InputFilename, SFM_READ);

	if (int e = infile.error()) {
		std::cout << "Error: Couldn't Open Input File (" << sf_error_number(e) << ")" << std::endl;
		return false;
	}

	// read file properties:
	unsigned int nChannels = infile.channels();
	unsigned int InputSampleRate = infile.samplerate();
	int InputFileFormat = infile.format();

	// detect if input format is a floating-point format:
	bool bFloat = false;
	bool bDouble = false;
	switch (InputFileFormat & SF_FORMAT_SUBMASK) {
	case SF_FORMAT_FLOAT:
		bFloat = true;
		break;
	case SF_FORMAT_DOUBLE:
		bDouble = true;
		break;
	}

	std::cout << "source file channels: " << nChannels << std::endl;
	std::cout << "input sample rate: " << InputSampleRate << "\noutput sample rate: " << ci.OutputSampleRate << std::endl;

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

	Fraction F = GetSimplifiedFraction(InputSampleRate, ci.OutputSampleRate);

	if (ci.bMinPhase) {
		if (F.numerator <= 4 && F.denominator <= 4) {
			if (F.numerator != F.denominator) { // oversample to improve filter performance
				F.numerator *= 8;	// 8x oversampling
				F.denominator *= 8;
			}
		}
	}

	FloatType ResamplingFactor = static_cast<FloatType>(ci.OutputSampleRate) / InputSampleRate;
	std::cout << "\nConversion ratio: " << ResamplingFactor
		<< " (" << F.numerator << ":" << F.denominator << ")" << std::endl;
	
	size_t BufferSize = (BUFFERSIZE / nChannels) * nChannels; // round down to integer multiple of nChannels (file may have odd number of channels!)
	assert(BUFFERSIZE >= BufferSize);

	FloatType inbuffer[BUFFERSIZE];
	FloatType outbuffer[BUFFERSIZE];

	sf_count_t count;
	sf_count_t SamplesRead = 0i64;
	FloatType PeakInputSample = 0.0;

	std::cout << "Scanning input file for peaks ..."; // to-do: can we read the PEAK chunk in floating-point files ?

	do {
		count = infile.read(inbuffer, BufferSize);
		SamplesRead += count;
		for (unsigned int s = 0; s < count; ++s) { // read all samples, without caring which channel they belong to
			PeakInputSample = max(PeakInputSample, abs(inbuffer[s]));
		}
	} while (count > 0);

	infile.seek(0i64, SEEK_SET); // rewind back to start of file

	std::cout << "Done\n";
	std::cout << "Peak input sample: " << std::fixed << PeakInputSample << " (" << 20 * log10(PeakInputSample) << " dBFS)" << std::endl;

	if (ci.bNormalize) { // echo Normalization settings to user
		auto prec = std::cout.precision();
		std::cout << "Normalizing to " << std::setprecision(2) << ci.Limit << std::endl;
		std::cout.precision(prec);
	}

	// Calculate filter parameters:
	int OverSampFreq = InputSampleRate * F.numerator; // eg 160 * 44100
	unsigned int minSampleRate = min(InputSampleRate, ci.OutputSampleRate);
	int TransitionWidth = minSampleRate / 22; // reasonable estimate for allowing transition width to scale with sample rate
//	TransitionWidth = max(TransitionWidth, 1700); // put a limit on how narrow the transition can be (too narrow will cause issues with low target samplerates)
	int ft = min(InputSampleRate, ci.OutputSampleRate) / 2.0 - TransitionWidth;

	// Make some filters: 
	// Huge Filters are used for complex ratios.
	// Medium filters used for simple ratios (ie 1 in numerator or denominator)

	FloatType HugeFilterTaps[FILTERSIZE_HUGE];
	int HugeFilterSize = FILTERSIZE_HUGE;
	makeLPF<FloatType>(HugeFilterTaps, HugeFilterSize, ft, OverSampFreq);
	applyKaiserWindow<FloatType>(HugeFilterTaps, HugeFilterSize, calcKaiserBeta(140));
	
	FloatType MedFilterTaps[FILTERSIZE_MEDIUM];
	int MedFilterSize = FILTERSIZE_MEDIUM;
	makeLPF<FloatType>(MedFilterTaps, MedFilterSize, ft, OverSampFreq);
	applyKaiserWindow<FloatType>(MedFilterTaps, MedFilterSize, calcKaiserBeta(195));
	
	if (ci.bMinPhase) {
		std::cout << "Using Minimum-Phase LPF" << std::endl;
		makeMinPhase<FloatType>(HugeFilterTaps, HugeFilterSize);
		makeMinPhase<FloatType>(MedFilterTaps, MedFilterSize);
	}
		
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
	int OutputFileFormat = ci.OutputFormat ? ci.OutputFormat : InputFileFormat;

	// if the minor (sub) format of OutputFileFormat is not set, attempt to use minor format of input file (as a last resort)
	if ((OutputFileFormat & SF_FORMAT_SUBMASK) == 0) {
		OutputFileFormat |= (InputFileFormat & SF_FORMAT_SUBMASK); // may not be valid subformat for new file format. 
	}

	// determine number of bits in output format, for Dithering:
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
		signalBits = 16; // to-do: what should it be for floating-point types ?
	}

	// confirm dithering options for user:
	if (ci.bDither) {
		auto prec = std::cout.precision();
		std::cout << "Generating " << std::setprecision(2) << ci.DitherAmount << " bits of dither for " << signalBits << "-bit output format";
		std::cout.precision(prec);
		if (ci.bAutoBlankingEnabled)
			std::cout << ", with auto-blanking";
		std::cout << std::endl;
	}

	// make a vector of ditherers (one ditherer for each channel):
	std::vector<Ditherer<FloatType>> Ditherers;
	for (unsigned int n = 0; n < nChannels; n++) {
		Ditherers.emplace_back(signalBits, ci.DitherAmount, ci.bAutoBlankingEnabled);
	}

	// Estimate initial gain:
	FloatType Gain = ci.bNormalize ? F.numerator * (ci.Limit / PeakInputSample) : F.numerator * ci.Limit;

	// Conditionally guard against filter overshoot:
	const FloatType OvershootCompensationFactor = 0.992; // Scaling factor to allow for filter overshoot
	if ((PeakInputSample > OvershootCompensationFactor) && !ci.bNormalize)
		Gain *= OvershootCompensationFactor;

	if (ci.bDither) { // allow headroom for dithering:
		FloatType DitherCompensation =
			(pow(2, signalBits - 1) - pow(2, ci.DitherAmount - 1)) / pow(2, signalBits - 1); // eg 32767/32768 = 0.999969 (-0.00027 dB)
		Gain *= DitherCompensation;
	}

	FloatType PeakOutputSample;
	SndfileHandle* pOutFile;
	//
	START_TIMER();

	do { // clipping detection loop (repeat if clipping detected)

		try { // Open output file:

			  // pOutFile needs to be dynamically allocated, because the only way to close file is to go out of scope 
			  // ... and we may need to overwrite file on subsequent pass:

			pOutFile = new SndfileHandle(ci.OutputFilename, SFM_WRITE, OutputFileFormat, nChannels, ci.OutputSampleRate);

			if (int e = pOutFile->error()) {
				std::cout << "Error: Couldn't Open Output File (" << sf_error_number(e) << ")" << std::endl;
				return false;
			}

			// if the minor (sub) format of OutputFileFormat is flac, and user has requested a specific compression level, set compression level:
			if (((OutputFileFormat & SF_FORMAT_FLAC) == SF_FORMAT_FLAC) && ci.bSetFlacCompression) {
				std::cout << "setting flac compression level to " << ci.flacCompressionLevel << std::endl;
				double cl = static_cast<double>(ci.flacCompressionLevel / 8.0); // there are 9 flac compression levels from 0-8. Normalize to 0-1.0
				pOutFile->command(SFC_SET_COMPRESSION_LEVEL, &cl, sizeof(cl));
			}

			// if the minor (sub) format of OutputFileFormat is vorbis, and user has requested a specific quality level, set quality level:
			if (((OutputFileFormat & SF_FORMAT_VORBIS) == SF_FORMAT_VORBIS) && ci.bSetVorbisQuality) {
				
				auto prec = std::cout.precision();
				std::cout.precision(1);
				std::cout << "setting vorbis quality level to " << ci.vorbisQuality << std::endl;
				std::cout.precision(prec);
				
				double cl = static_cast<double>((1.0-ci.vorbisQuality) / 11.0); // Normalize from (-1 to 10), to (1.0 to 0) ... why is it backwards ?
				pOutFile->command(SFC_SET_COMPRESSION_LEVEL, &cl, sizeof(cl));
			}

		}

		catch (std::bad_alloc& b) {
			std::cout << "Error: Couldn't Open Output File (memory allocation problem)" << std::endl;
			return false;
		}

		std::cout << "Converting ...";
		unsigned int DecimationIndex = 0;
		unsigned int OutBufferIndex = 0;
		PeakOutputSample = 0.0;

		if (F.numerator == 1 && F.denominator == 1) { // no change to sample rate; format conversion only

			std::cout << " No change to sample rate" << std::endl;
			do { // Read and process blocks of samples until the end of file is reached

				count = infile.read(inbuffer, BufferSize);

				for (unsigned int s = 0; s < count; s += nChannels) {

					for (int Channel = 0; Channel < nChannels; Channel++)
					{
						FloatType OutputSample = ci.bDither ?
							Ditherers[Channel].Dither(Gain * inbuffer[s + Channel]) :
							Gain * inbuffer[s + Channel];

						outbuffer[OutBufferIndex + Channel] = OutputSample;
						PeakOutputSample = max(abs(PeakOutputSample), abs(OutputSample));
					}

					OutBufferIndex += nChannels;
					if (OutBufferIndex == BufferSize) {
						OutBufferIndex = 0;
						pOutFile->write(outbuffer, BufferSize);
					}

				} // ends loop over s
			} while (count > 0);
		}


		else if (F.numerator == 1 && F.denominator != 1) { // Decimate Only
			do { // Read and process blocks of samples until the end of file is reached
				count = infile.read(inbuffer, BufferSize);

				for (unsigned int s = 0; s < count; s += nChannels) {

					for (int Channel = 0; Channel < nChannels; Channel++)
						MedFilters[Channel].put(inbuffer[s + Channel]); // inject a source sample

					if (DecimationIndex == 0) { // Decimate
						for (int Channel = 0; Channel < nChannels; Channel++) {

							FloatType OutputSample = ci.bDither ?
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

							FloatType OutputSample = ci.bDither ?
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

								FloatType OutputSample = ci.bDither ?
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

		  // Tail:
		if (OutBufferIndex != 0) {
			pOutFile->write(outbuffer, OutBufferIndex); // finish writing whatever remains in the buffer 
		}

		std::cout << "Done" << std::endl;
		
		auto prec = std::cout.precision();
		std::cout << "Peak output sample: " << std::setprecision(6) << PeakOutputSample << " (" << 20 * log10(PeakOutputSample) << " dBFS)" << std::endl;
		std::cout.precision(prec);

		delete pOutFile; // Close output file

		if (bDouble || bFloat) // To-do: Confirm assumption upon which the following statement is built:
			break; // Clipping is not a concern with Floating-Point formats. 

		// Test for clipping:
		if (PeakOutputSample > ci.Limit) { // Clipping !

			FloatType GainAdjustment;

			if (ci.bDither) {
				FloatType DitherLevel = pow(2, ci.DitherAmount - 1) / pow(2, signalBits - 1);
				// take Ditherlevel out of calculations, as dither is NOT subject to gain control:
				GainAdjustment = (ci.Limit - DitherLevel) / (PeakOutputSample - DitherLevel);
			}

			else
				GainAdjustment = ci.Limit / PeakOutputSample;

			Gain *= GainAdjustment;
			std::cout << "\nClipping detected !" << std::endl;
			std::cout << "Re-doing with " << 20 * log10(GainAdjustment) << " dB gain adjustment" << std::endl;
			infile.seek(0i64, SEEK_SET);
		}

	} while (PeakOutputSample > ci.Limit);

	STOP_TIMER();
	return true;
} // ends Convert()

// gcd() - greatest common divisor:
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

//  GetSimplifiedFraction() - turns a sample-rate ratio into a fraction:
Fraction GetSimplifiedFraction(int InputSampleRate, int OutputSampleRate)			// eg 44100, 48000
{
	Fraction f;
	f.numerator = (OutputSampleRate / gcd(InputSampleRate, OutputSampleRate));		// L (eg 160)
	f.denominator = (InputSampleRate / gcd(InputSampleRate, OutputSampleRate));		// M (eg 147)
	return f;
}


// The following functions are used for parsing commandline parameters:

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

void getCmdlineParam(char** begin, char** end, const std::string& OptionName, int& nParameter)
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

//
