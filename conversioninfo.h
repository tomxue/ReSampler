/*
* Copyright (C) 2016 - 2017 Judd Niemann - All Rights Reserved
* You may use, distribute and modify this code under the
* terms of the GNU Lesser General Public License, version 2.1
*
* You should have received a copy of GNU Lesser General Public License v2.1
* with this file. If not, please refer to: https://github.com/jniemann66/ReSampler
*/

#ifndef CONVERSIONINFO_H
#define CONVERSIONINFO_H

// conversioninfo.h:
// defines the ConversionInfo struct,
// for holding conversion parameters.

#include "ditherer.h"

#include <vector>
#include <string>

typedef enum {
	relaxed,
	normal,
	steep,
	custom
} LPFMode;

// The following functions are used for parsing commandline parameters:

void getCmdlineParam(char** begin, char** end, const std::string& optionName, std::string& parameter)
{
	parameter.clear();
	char** it = std::find(begin, end, optionName);
	if (it != end)	// found option
		if (++it != end) // found parameter after option
			parameter = *it;
}

void getCmdlineParam(char** begin, char** end, const std::string& optionName, unsigned int& nParameter)
{
	nParameter = 0;
	char** it = std::find(begin, end, optionName);
	if (it != end)	// found option
		if (++it != end) // found parameter after option
			nParameter = atoi(*it);
}

void getCmdlineParam(char** begin, char** end, const std::string& optionName, int& nParameter)
{
	nParameter = 0;
	char** it = std::find(begin, end, optionName);
	if (it != end)	// found option
		if (++it != end) // found parameter after option
			nParameter = atoi(*it);
}


void getCmdlineParam(char** begin, char** end, const std::string& optionName, double& parameter)
{
	parameter = 0.0;
	char** it = std::find(begin, end, optionName);
	if (it != end)	// found option
		if (++it != end) // found parameter after option
			parameter = atof(*it);
}

bool findCmdlineOption(char** begin, char** end, const std::string& option) {
	return (std::find(begin, end, option) != end);
}

// struct ConversionInfo : structure for holding all the parameters required for a conversion job

struct ConversionInfo
{
	std::string inputFilename;
	std::string outputFilename;
	unsigned int inputSampleRate;
	unsigned int outputSampleRate;
	double gain;
	double limit;
	bool bUseDoublePrecision;
	bool bNormalize;
	double normalizeAmount;
	int outputFormat;
	std::string outBitFormat;
	bool bDither;
	double ditherAmount;
	int ditherProfileID;
	bool bAutoBlankingEnabled;
	bool bDelayTrim;
	bool bMinPhase;
	bool bSetFlacCompression;
	int flacCompressionLevel;
	bool bSetVorbisQuality;
	double vorbisQuality;
	bool disableClippingProtection;
	LPFMode lpfMode;
	double lpfCutoff;
	double lpfTransitionWidth;
	bool bUseSeed;
	int seed;
	bool dsfInput;
	bool dffInput;
	bool bEnablePeakDetection;
	bool bMultiThreaded;
	bool bRf64;
	bool bNoPeakChunk;
	bool bWriteMetaData;
	int maxStages;
	bool bSingleStage;
	bool bShowStages;
	int overSamplingFactor;
	bool bBadParams;
	std::string appName;

	bool fromCmdLineArgs(int argc, char* argv[]);
	std::string toCmdLineArgs();
};


inline std::string ConversionInfo::toCmdLineArgs() {
	std::vector<std::string> args;
	std::string result;

	args.push_back("-i");
	args.push_back(inputFilename);
	args.push_back("-o");
	args.push_back(outputFilename);
	args.push_back("-r");
	args.push_back(std::to_string(outputSampleRate));

	if(bUseDoublePrecision)
		args.push_back("--doubleprecision");

	if(bNormalize) {
		args.push_back("-n");
		args.push_back(std::to_string(normalizeAmount));
	}

	if(bMinPhase)
		args.push_back("--minphase");

	if (lpfMode == custom) {
		args.push_back("--lpf-cutoff");
		args.push_back(std::to_string(lpfCutoff));
		args.push_back("--lpf-transition");
		args.push_back(std::to_string(lpfTransitionWidth));
	}

	if (maxStages == 1) {
		args.push_back("--maxStages");
		args.push_back(std::to_string(maxStages));
	}

	for(auto it = args.begin(); it != args.end(); it++) {
		result.append(*it);
		if(it != std::prev(args.end()))
			result.append(" ");
	}

	return result;
}

// fromCmdLineArgs()
// Return value indicates whether caller should continue execution (ie true: continue, false: terminate)
// Some commandline options (eg --version) should result in termination, but not error.
// unacceptable parameters are indicated by setting bBadParams to true

inline bool ConversionInfo::fromCmdLineArgs(int argc, char* argv[]) {
	// initialize defaults:
	inputFilename.clear();
	outputFilename.clear();
	outBitFormat.clear();
	outputFormat = 0;
	outputSampleRate = 44100;
	normalizeAmount = 1.0;
	ditherAmount = 1.0;
	flacCompressionLevel = 5;
	vorbisQuality = 3;
	ditherProfileID = DitherProfileID::standard;

	////////////////////////////////////////////////////////////////////
	// core parameters:
	getCmdlineParam(argv, argv + argc, "-i", inputFilename);
	getCmdlineParam(argv, argv + argc, "-o", outputFilename);
	getCmdlineParam(argv, argv + argc, "-r", outputSampleRate);
	getCmdlineParam(argv, argv + argc, "-b", outBitFormat);

	// gain
	if (findCmdlineOption(argv, argv + argc, "--gain")) {
		getCmdlineParam(argv, argv + argc, "--gain", gain);
	}
	else {
		gain = 1.0; // default
	}

	// double precision switch:
	bUseDoublePrecision = findCmdlineOption(argv, argv + argc, "--doubleprecision");

	// normalize option and parameter:
	bNormalize = findCmdlineOption(argv, argv + argc, "-n");
	if (bNormalize) {
		getCmdlineParam(argv, argv + argc, "-n", normalizeAmount);
		if (normalizeAmount <= 0.0)
			normalizeAmount = 1.0;
		if (normalizeAmount > 1.0)
			std::cout << "\nWarning: Normalization factor greater than 1.0 - THIS WILL CAUSE CLIPPING !!\n" << std::endl;
		limit = normalizeAmount;
	}
	else {
		limit = 1.0; // default
	}

	// dither option and parameter:
	bDither = findCmdlineOption(argv, argv + argc, "--dither");
	if (bDither) {
		getCmdlineParam(argv, argv + argc, "--dither", ditherAmount);
		if (ditherAmount <= 0.0)
			ditherAmount = 1.0;
	}

	// auto-blanking option (for dithering):
	bAutoBlankingEnabled = findCmdlineOption(argv, argv + argc, "--autoblank");

	// ns option to determine dither Profile:
	if (findCmdlineOption(argv, argv + argc, "--ns")) {
		getCmdlineParam(argv, argv + argc, "--ns", ditherProfileID);
		if (ditherProfileID < 0)
			ditherProfileID = 0;
		if (ditherProfileID >= DitherProfileID::end)
			ditherProfileID = getDefaultNoiseShape(outputSampleRate);
	}
	else {
		ditherProfileID = getDefaultNoiseShape(outputSampleRate);
	}

	// --flat-tpdf option (takes precedence over --ns)
	if (findCmdlineOption(argv, argv + argc, "--flat-tpdf")) {
		ditherProfileID = DitherProfileID::flat;
	}

	// seed option and parameter:
	bUseSeed = findCmdlineOption(argv, argv + argc, "--seed");
	seed = 0;
	if (bUseSeed) {
		getCmdlineParam(argv, argv + argc, "--seed", seed);
	}

	// delay trim (group delay compensation)
	bDelayTrim = !findCmdlineOption(argv, argv + argc, "--noDelayTrim");

	// minimum-phase option:
	bMinPhase = findCmdlineOption(argv, argv + argc, "--minphase");

	// flacCompression option and parameter:
	bSetFlacCompression = findCmdlineOption(argv, argv + argc, "--flacCompression");
	if (bSetFlacCompression) {
		getCmdlineParam(argv, argv + argc, "--flacCompression", flacCompressionLevel);
		if (flacCompressionLevel < 0)
			flacCompressionLevel = 0;
		if (flacCompressionLevel > 8)
			flacCompressionLevel = 8;
	}

	// vorbisQuality option and parameter:
	bSetVorbisQuality = findCmdlineOption(argv, argv + argc, "--vorbisQuality");
	if (bSetVorbisQuality) {
		getCmdlineParam(argv, argv + argc, "--vorbisQuality", vorbisQuality);
		if (vorbisQuality < -1)
			vorbisQuality = -1;
		if (vorbisQuality > 10)
			vorbisQuality = 10;
	}

	// noClippingProtection option:
	disableClippingProtection = findCmdlineOption(argv, argv + argc, "--noClippingProtection");

	// default cutoff and transition width:
	lpfMode = normal;
	lpfCutoff = 100.0 * (10.0 / 11.0);
	lpfTransitionWidth = 100.0 - lpfCutoff;

	// relaxedLPF option:
	if (findCmdlineOption(argv, argv + argc, "--relaxedLPF")) {
		lpfMode = relaxed;
		lpfCutoff = 100.0 * (21.0 / 22.0);				// late cutoff
		lpfTransitionWidth = 2 * (100.0 - lpfCutoff); // wide transition (double-width)  
	}

	// steepLPF option:
	if (findCmdlineOption(argv, argv + argc, "--steepLPF")) {
		lpfMode = steep;
		lpfCutoff = 100.0 * (21.0 / 22.0);				// late cutoff
		lpfTransitionWidth = 100.0 - lpfCutoff;       // steep transition  
	}

	// custom LPF cutoff frequency:
	if (findCmdlineOption(argv, argv + argc, "--lpf-cutoff")) {
		getCmdlineParam(argv, argv + argc, "--lpf-cutoff", lpfCutoff);
		lpfMode = custom;
		lpfCutoff = std::max(1.0, std::min(lpfCutoff, 99.9));

		// custom LPF transition width:
		if (findCmdlineOption(argv, argv + argc, "--lpf-transition")) {
			getCmdlineParam(argv, argv + argc, "--lpf-transition", lpfTransitionWidth);
		}
		else {
			lpfTransitionWidth = 100 - lpfCutoff; // auto mode
		}
		lpfTransitionWidth = std::max(0.1, std::min(lpfTransitionWidth, 400.0));
	}

	// multithreaded option:
	bMultiThreaded = findCmdlineOption(argv, argv + argc, "--mt");

	// rf64 option:
	bRf64 = findCmdlineOption(argv, argv + argc, "--rf64");

	// noPeakChunk option:
	bNoPeakChunk = findCmdlineOption(argv, argv + argc, "--noPeakChunk");

	// noMetadata option:
	bWriteMetaData = !findCmdlineOption(argv, argv + argc, "--noMetadata");

	// maxStages:
	if (findCmdlineOption(argv, argv + argc, "--maxStages")) {
		getCmdlineParam(argv, argv + argc, "--maxStages", maxStages);
		if (maxStages < 1)
			maxStages = 1;
		if (maxStages > 10)
			maxStages = 10;
	}
	else {
		maxStages = 3; // default;
	}

	// single stage:
	bSingleStage = findCmdlineOption(argv, argv + argc, "--singleStage");

	// showStages option:
	bShowStages = findCmdlineOption(argv, argv + argc, "--showStages");

	// test for bad parameters:
	bBadParams = false;
	if (outputFilename.empty()) {
		if (inputFilename.empty()) {
			std::cout << "Error: Input filename not specified" << std::endl;
			bBadParams = true;
		}
		else {
			std::cout << "Output filename not specified" << std::endl;
			outputFilename = inputFilename;
			if (outputFilename.find(".") != std::string::npos) {
				auto dot = outputFilename.find_last_of(".");
				outputFilename.insert(dot, "(converted)");
			}
			else {
				outputFilename.append("(converted)");
			}
			std::cout << "defaulting to: " << outputFilename << "\n" << std::endl;
		}
	}

	else if (outputFilename == inputFilename) {
		std::cout << "\nError: Input and Output filenames cannot be the same" << std::endl;
		bBadParams = true;
	}

	if (outputSampleRate == 0) {
		std::cout << "Error: Target sample rate not specified" << std::endl;
		bBadParams = true;
	}

	if (bBadParams) {
		std::cout << strUsage << std::endl;
		return false;
	}
	return true;
}

static_assert(std::is_copy_constructible<ConversionInfo>::value,
	"ConversionInfo must be copy constructible");

static_assert(std::is_copy_assignable<ConversionInfo>::value,
	"ConversionInfo must be copy assignable");

#endif // CONVERSIONINFO_H

