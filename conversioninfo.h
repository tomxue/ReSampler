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

#include <iostream>
#include <vector>
#include <string>

typedef enum {
	relaxed,
	normal,
	steep,
	custom
} LPFMode;

// The following functions are used for fetching commandline parameters:

// proposal to make the parser more permissive:
// lhs: --flatTpdf
// rhs: --flat-tpdf

// on both sides,
// remove all hyphens after first alphanum character
// convert to lowercase

// lhs: --flattpdf
// rhs: --flattpdf


// for numbers:
template<typename T>
bool getCmdlineParam(char** begin, char** end, const std::string& optionName, T& parameter) {
	bool found = false;
	char** it = std::find(begin, end, optionName);
	if (it != end) {
		found = true;
		if (++it != end)
			parameter = atof(*it);
	}
	return found;
}

// for strings:
bool getCmdlineParam(char** begin, char** end, const std::string& optionName, std::string& parameter)
{
	bool found = false;
	char** it = std::find(begin, end, optionName);
	if (it != end) {
		found = true;
		if (++it != end)
			parameter = *it;
	}
	return found;
}

// switch only (no parameter)
bool getCmdlineParam(char** begin, char** end, const std::string& optionName)
{
	return (std::find(begin, end, optionName) != end);
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
	outputSampleRate = 0;
	normalizeAmount = 1.0;
	limit = 1.0;
	ditherAmount = 1.0;
	ditherProfileID = DitherProfileID::standard;
	flacCompressionLevel = 5;
	vorbisQuality = 3;
	gain = 1.0;
	seed = 0;
	maxStages = 3;
	lpfMode = normal;
	lpfCutoff = 100.0 * (10.0 / 11.0);
	lpfTransitionWidth = 100.0 - lpfCutoff;
	dsfInput = false;
	dffInput = false;

	// get core parameters:
	getCmdlineParam(argv, argv + argc, "-i", inputFilename);
	getCmdlineParam(argv, argv + argc, "-o", outputFilename);
	getCmdlineParam(argv, argv + argc, "-r", outputSampleRate);
	getCmdlineParam(argv, argv + argc, "-b", outBitFormat);

	// get extended parameters
	getCmdlineParam(argv, argv + argc, "--gain", gain);
	bUseDoublePrecision = getCmdlineParam(argv, argv + argc, "--doubleprecision");
	bNormalize = getCmdlineParam(argv, argv + argc, "-n", normalizeAmount);
	bDither = getCmdlineParam(argv, argv + argc, "--dither", ditherAmount);
	ditherProfileID = getDefaultNoiseShape(outputSampleRate);
	getCmdlineParam(argv, argv + argc, "--ns", ditherProfileID);
	ditherProfileID = getCmdlineParam(argv, argv + argc, "--flat-tpdf") ? DitherProfileID::flat : ditherProfileID;
	bAutoBlankingEnabled = getCmdlineParam(argv, argv + argc, "--autoblank");
	bUseSeed = getCmdlineParam(argv, argv + argc, "--seed", seed);
	bDelayTrim = !getCmdlineParam(argv, argv + argc, "--noDelayTrim");
	bMinPhase = getCmdlineParam(argv, argv + argc, "--minphase");
	bSetFlacCompression = getCmdlineParam(argv, argv + argc, "--flacCompression", flacCompressionLevel);
	bSetVorbisQuality = getCmdlineParam(argv, argv + argc, "--vorbisQuality", vorbisQuality);
	bMultiThreaded = getCmdlineParam(argv, argv + argc, "--mt");
	bRf64 = getCmdlineParam(argv, argv + argc, "--rf64");
	bNoPeakChunk = getCmdlineParam(argv, argv + argc, "--noPeakChunk");
	bWriteMetaData = !getCmdlineParam(argv, argv + argc, "--noMetadata");
	getCmdlineParam(argv, argv + argc, "--maxStages", maxStages);
	bSingleStage = getCmdlineParam(argv, argv + argc, "--singleStage");
	bShowStages = getCmdlineParam(argv, argv + argc, "--showStages");

	// LPFilter settings:
	if (getCmdlineParam(argv, argv + argc, "--relaxedLPF")) {
		lpfMode = relaxed;
		lpfCutoff = 100.0 * (21.0 / 22.0);				// late cutoff
		lpfTransitionWidth = 2 * (100.0 - lpfCutoff); // wide transition (double-width)  
	}

	if (getCmdlineParam(argv, argv + argc, "--steepLPF")) {
		lpfMode = steep;
		lpfCutoff = 100.0 * (21.0 / 22.0);				// late cutoff
		lpfTransitionWidth = 100.0 - lpfCutoff;       // steep transition  
	}

	if (getCmdlineParam(argv, argv + argc, "--lpf-cutoff", lpfCutoff)) { // custom LPF cutoff frequency
		lpfMode = custom;
		if (!getCmdlineParam(argv, argv + argc, "--lpf-transition", lpfTransitionWidth)) {
			lpfTransitionWidth = 100 - lpfCutoff; // auto mode
		}
	}

	// constraining functions:
	auto constrainDouble = [](double& val, double minVal, double maxVal) {
		val = std::max(minVal, std::min(val, maxVal));
	};

	auto constrainInt = [](int& val, int minVal, int maxVal) {
		val = std::max(minVal, std::min(val, maxVal));
	};

	// set constraints:
	constrainInt(flacCompressionLevel, 0, 8);
	constrainDouble(vorbisQuality, -1, 10);
	constrainInt(maxStages, 1, 10);
	constrainDouble(lpfCutoff, 1.0, 99.9);
	constrainDouble(lpfTransitionWidth, 0.1, 400.0);

	if (bNormalize) {
		if (normalizeAmount <= 0.0)
			normalizeAmount = 1.0;
		if (normalizeAmount > 1.0)
			std::cout << "\nWarning: Normalization factor greater than 1.0 - THIS WILL CAUSE CLIPPING !!\n" << std::endl;
		limit = normalizeAmount;
	}
	
	if (bDither) {
		if (ditherAmount <= 0.0)
			ditherAmount = 1.0;
	}
	
	if (ditherProfileID < 0)
		ditherProfileID = 0;
	if (ditherProfileID >= DitherProfileID::end)
		ditherProfileID = getDefaultNoiseShape(outputSampleRate);

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

