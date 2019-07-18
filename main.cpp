/*
* Copyright (C) 2016 - 2019 Judd Niemann - All Rights Reserved.
* You may use, distribute and modify this code under the
* terms of the GNU Lesser General Public License, version 2.1
*
* You should have received a copy of GNU Lesser General Public License v2.1
* with this file. If not, please refer to: https://github.com/jniemann66/ReSampler
*/

// main.cpp : main entry point

#include "ReSampler.h"
#include "conversioninfo.h"
#include "dsf.h"
#include "dff.h"

int main(int argc, char * argv[])
{

#ifdef COMPILING_ON_ANDROID
	std::cout.rdbuf(new androidbuf(ANDROID_LOG_INFO, "ReSampler"));
	std::cerr.rdbuf(new androidbuf(ANDROID_LOG_ERROR, "ReSampler"));
	// register proposed android cleanup function:
	// atexit(androidCleanup);
#endif

	// test for global options
	if (parseGlobalOptions(argc, argv)) {

// todo: remove these if atexit(androidClenup) works properly:
#ifdef COMPILING_ON_ANDROID
		delete std::cout.rdbuf(0);
		delete std::cerr.rdbuf(0);
#endif
// ---
		return EXIT_SUCCESS;
	}


	// ConversionInfo instance to hold parameters
	ConversionInfo ci;

	// get path/name of this app
#ifdef __APPLE__
	char pathBuf[PROC_PIDPATHINFO_MAXSIZE];
	pid_t pid = getpid();
	if (proc_pidpath(pid, pathBuf, sizeof(pathBuf)) == 0) {
		ci.appName.assign(pathBuf);
	}
#else
	ci.appName = argv[0];
#endif

	ci.overSamplingFactor = 1;

	// get conversion parameters
	ci.fromCmdLineArgs(argc, argv);
	if (ci.bBadParams) {

#ifdef COMPILING_ON_ANDROID
		delete std::cout.rdbuf(0);
		delete std::cerr.rdbuf(0);
#endif
		std::cout << strUsage << std::endl;
		return EXIT_FAILURE;
	}

	// query build version AND cpu
	if (!showBuildVersion()) {

#ifdef COMPILING_ON_ANDROID
		delete std::cout.rdbuf(0);
		delete std::cerr.rdbuf(0);
#endif
		return EXIT_FAILURE; // can't continue (CPU / build mismatch)
	}

	// echo filenames to user
	std::cout << "Input file: " << ci.inputFilename << std::endl;
	std::cout << "Output file: " << ci.outputFilename << std::endl;

	if (ci.disableClippingProtection) {
		std::cout << "clipping protection disabled " << std::endl;
	}

	// Isolate the file extensions
	std::string inFileExt;
	std::string outFileExt;
	if (ci.inputFilename.find_last_of('.') != std::string::npos)
		inFileExt = ci.inputFilename.substr(ci.inputFilename.find_last_of('.') + 1);
	if (ci.outputFilename.find_last_of('.') != std::string::npos)
		outFileExt = ci.outputFilename.substr(ci.outputFilename.find_last_of('.') + 1);

	// detect dsf or dff format
	ci.dsfInput = (inFileExt == "dsf");
	ci.dffInput = (inFileExt == "dff");

	// detect csv output
	ci.csvOutput = (outFileExt == "csv");

	if (ci.csvOutput) {
		std::cout << "Outputting to csv format" << std::endl;
	}
	else {
		if (!ci.outBitFormat.empty()) {  // new output bit format requested
			ci.outputFormat = determineOutputFormat(outFileExt, ci.outBitFormat);
			if (ci.outputFormat) {
				std::cout << "Changing output bit format to " << ci.outBitFormat << std::endl;
			}
			else { // user-supplied bit format not valid; try choosing appropriate format
				std::string outBitFormat;
				determineBestBitFormat(outBitFormat, ci);
				ci.outputFormat = determineOutputFormat(outFileExt, outBitFormat);
				if (ci.outputFormat) {
					ci.outBitFormat = outBitFormat;
					std::cout << "Changing output bit format to " << ci.outBitFormat << std::endl;
				}
				else {
					std::cout << "Warning: NOT Changing output file bit format !" << std::endl;
					ci.outputFormat = 0; // back where it started
				}
			}
		}

		if (outFileExt != inFileExt) { // file extensions differ, determine new output format:

			std::string outBitFormat{ci.outBitFormat};
			if (ci.outBitFormat.empty()) { // user changed file extension only. Attempt to choose appropriate output sub format:
				std::cout << "Output Bit Format not specified" << std::endl;
				determineBestBitFormat(outBitFormat, ci);
			}
			ci.outputFormat = determineOutputFormat(outFileExt, outBitFormat);
			if (ci.outputFormat) {
				ci.outBitFormat = outBitFormat;
				std::cout << "Changing output file format to " << outFileExt << std::endl;
			} else { // cannot determine subformat of output file
				std::cout << "Warning: NOT Changing output file format ! (extension different, but format will remain the same)" << std::endl;
			}
		}
	}
	try {

		if (ci.bUseDoublePrecision) {

#ifdef USE_QUADMATH
			std::cout << "Using quadruple-precision for calculations.\n";
#else
			std::cout << "Using double precision for calculations." << std::endl;
#endif
			if (ci.dsfInput) {
				ci.bEnablePeakDetection = false;
				return convert<DsfFile, double>(ci) ? EXIT_SUCCESS : EXIT_FAILURE;
			}
			else if (ci.dffInput) {
				ci.bEnablePeakDetection = false;
				return convert<DffFile, double>(ci) ? EXIT_SUCCESS : EXIT_FAILURE;
			}
			else {
				ci.bEnablePeakDetection = true;
				return convert<SndfileHandle, double>(ci) ? EXIT_SUCCESS : EXIT_FAILURE;
			}
		}

		else {

#ifdef USE_QUADMATH
			std::cout << "Using quadruple-precision for calculations.\n";
#endif
			if (ci.dsfInput) {
				ci.bEnablePeakDetection = false;
				return convert<DsfFile, float>(ci) ? EXIT_SUCCESS : EXIT_FAILURE;
			}
			else if (ci.dffInput) {
				ci.bEnablePeakDetection = false;
				return convert<DffFile, float>(ci) ? EXIT_SUCCESS : EXIT_FAILURE;
			}
			else {
				ci.bEnablePeakDetection = true;
				return convert<SndfileHandle, float>(ci) ? EXIT_SUCCESS : EXIT_FAILURE;
			}
		}

	} //ends try block

	catch (const std::exception& e) {
		std::cerr << "fatal error: " << e.what();

#ifdef COMPILING_ON_ANDROID
		delete std::cout.rdbuf(0);
		delete std::cerr.rdbuf(0);
#endif
		return EXIT_FAILURE;
	}
}
