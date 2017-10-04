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

#include <vector>
#include <string>

typedef enum {
	relaxed,
	normal,
	steep,
	custom
} LPFMode;

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
	bool bShowStages;
	int overSamplingFactor;
	std::string appName;
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

static_assert(std::is_copy_constructible<ConversionInfo>::value,
	"ConversionInfo must be copy constructible");

static_assert(std::is_copy_assignable<ConversionInfo>::value,
	"ConversionInfo must be copy assignable");

#endif // CONVERSIONINFO_H

