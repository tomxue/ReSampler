/*
* Copyright (C) 2016 - 2017 Judd Niemann - All Rights Reserved
* You may use, distribute and modify this code under the
* terms of the GNU Lesser General Public License, version 2.1
*
* You should have received a copy of GNU Lesser General Public License v2.1
* with this file. If not, please refer to: https://github.com/jniemann66/ReSampler
*/

#ifndef DITHERER_H
#define DITHERER_H 1

// Ditherer.h
// defines Ditherer class, for adding dither and noise-shaping to input samples

// configuration:
//#define TEST_FILTER // if defined, this will result in the input signal being ignored; output the noise only. (Used for evaluating filters.)
#define MAX_FIR_FILTER_SIZE 24

#include <cmath>
#include <random>
#include "biquad.h"
#include "noiseshape.h"

typedef enum {
	bypass,
	cascadedBiquad,
	fir
} FilterType;

typedef enum {
	flatTPDF,
	slopedTPDF,
	RPDF,
	GPDF,
	impulse
} NoiseGeneratorType;

typedef enum {
	flat,
	standard,
	Wannamaker3tap,
	Wannamaker9tap,
	Wannamaker24tap,
	HighShibata44k,
	ModEWeighted44k,
	Lipshitz44k,
	ImpEWeighted44k,
	Experimental1,
	Experimental2,
	rpdf,
	end
} DitherProfileID;

typedef struct {
	DitherProfileID id;
	const char* name;
	NoiseGeneratorType noiseGeneratorType;
	FilterType filterType;
	int intendedSampleRate;
	int N;
	const double* coeffs;
	bool bUseFeedback;
} DitherProfile;

DitherProfile ditherProfileList[] = {

	// id, name, noiseGeneratorType, filterType, intendedSampleRate, N, coeffs, bUseFeedback

	{flat, "flat tpdf", flatTPDF, bypass, 44100, 1, noiseShaperPassThrough, false},
	{standard, "standard", slopedTPDF, fir, 44100, 10, std_44, true },
	{Wannamaker3tap, "Wannamaker 3-tap",flatTPDF, fir, 44100, 3, wan3, true},
	{Wannamaker9tap, "Wannamaker 9-tap",flatTPDF, fir, 44100, 9, wan9, true},
	{Wannamaker24tap, "Wannamaker 24-tap",flatTPDF, fir, 44100, 24, wan24, true},
	{HighShibata44k, "High Shibata 44k",flatTPDF, fir, 44100, 20, highShib44, true},
	{ModEWeighted44k, "Modified E-Weighted",flatTPDF, fir, 44100, 9, modew44, true},
	{Lipshitz44k, "Lipshitz",flatTPDF, fir, 44100, 5, lips44, true},
	{ImpEWeighted44k, "Improved E-Weighted",flatTPDF, fir, 44100, 9, impew44, true},
	{Experimental1, "Experimental 1",flatTPDF, fir, 44100, 20, experimental1, true },
	{Experimental2, "Experimental 2",slopedTPDF, fir, 44100, 11, experimental2, true },
	{rpdf,"flat rectangular pdf", RPDF, bypass, 44100, 1, noiseShaperPassThrough, false}
};

template<typename FloatType>
class Ditherer  
{
public:
	// Constructor:
	// signalBits: number of bits of the target bitformat
	// ditherBits: number of bits of dither to add, and doesn't have to be an integer
	// bAutoBlankingEnabled: if true, enable auto-blanking of dither (on Silence)
	// seed: seed for PRNG
	// filterID: noise-shaping filter to use

	Ditherer(unsigned int signalBits, FloatType ditherBits, bool bAutoBlankingEnabled, int seed, DitherProfileID ditherProfileID = standard) :
		signalBits(signalBits),
		ditherBits(ditherBits),
		bAutoBlankingEnabled(bAutoBlankingEnabled),
		selectedDitherProfile(ditherProfileList[ditherProfileID]),
		seed(seed),
		Z1(0),
		randGenerator(seed),		// initialize (seed) RNG
		dist(0, randMax),		// set the range of the random number distribution
		gain(1.0),
		bUseErrorFeedback(ditherProfileList[ditherProfileID].bUseFeedback),
		bPulseEmitted(false)
	{
		// general parameters:
		maxSignalMagnitude = static_cast<FloatType>((1 << (signalBits - 1)) - 1); // note the -1 : match 32767 scaling factor for 16 bit !
		reciprocalSignalMagnitude = 1.0 / maxSignalMagnitude; // value of LSB in target format
		maxDitherScaleFactor = (FloatType)pow(2, ditherBits - 1) / maxSignalMagnitude / (FloatType)randMax;
		oldRandom = 0;

		// set-up noise generator:
		switch (selectedDitherProfile.noiseGeneratorType) {
		case flatTPDF:
			noiseGenerator = &Ditherer::noiseGeneratorFlatTPDF;
			break;
		case RPDF:
			noiseGenerator = &Ditherer::noiseGeneratorRPDF;
			break;
		case GPDF:
			noiseGenerator = &Ditherer::noiseGeneratorGPDF;
			break;
		case impulse:
			noiseGenerator = &Ditherer::noiseGeneratorImpulse;
			break;
		case slopedTPDF:
		default:
			noiseGenerator = &Ditherer::noiseGeneratorSlopedTPDF;
		}

		// set-up filter type:
		switch (selectedDitherProfile.filterType) {
		case bypass:
			noiseShapingFilter = &Ditherer::noiseShaperPassThrough;
			break;
		case fir:
			noiseShapingFilter = &Ditherer::noiseShaperFIR;
			break;
		case cascadedBiquad:
		default:
			noiseShapingFilter = &Ditherer::noiseShaperCascadedBiquad;
		}

		//// IIR-specific stuff:

		f1.setCoeffs(
			1,
			-1.584512777183064e+00,
			7.706380573200580e-01,
			0,
			0
		);

		f2.setCoeffs(
			1,
			6.656411088446591e-02,  3.092585934772854e-01,  7.551346335428704e-01,  1.491544856915262e-01
		);

		f3.setCoeffs(
			1,
			6.903896840936654e-01,  6.221635814920810e-01,  1.353784987738887e+00,  6.659957897439557e-01
		);

		// FIR-specific stuff:
		
		const FloatType scale = 1.0;
		FIRLength = selectedDitherProfile.N;
		for (int n = 0; n < FIRLength; ++n) {
			FIRCoeffs[n] = scale * selectedDitherProfile.coeffs[n];
		}

		memset(FIRHistory, 0, MAX_FIR_FILTER_SIZE * sizeof(FloatType));

		// set-up Auto-blanking:
		if (bAutoBlankingEnabled) {	// initial state: silence
			ditherScaleFactor = 0.0;
		}
		else {	// initial state: dithering
			ditherScaleFactor = maxDitherScaleFactor; 
		}

		autoBlankLevelThreshold = 1.0 / pow(2, 32); // 1 LSB of 32-bit digital
		autoBlankTimeThreshold = 30000; // number of zero samples before activating autoblank
		autoBlankDecayCutoff = 0.7 * reciprocalSignalMagnitude / randMax;
		zeroCount = 0;
		
	} // Ends Constructor 

	void adjustGain(FloatType factor) {
		gain *= factor;
		maxDitherScaleFactor = gain * pow(2, ditherBits - 1) / maxSignalMagnitude / randMax;
	}

	void reset() {
		// reset filters
		f1.reset();
		f2.reset();
		f3.reset();

		memset(FIRHistory, 0, MAX_FIR_FILTER_SIZE * sizeof(FloatType));
		
		// re-seed PRNG
		randGenerator.seed(seed);
		oldRandom = 0;
		Z1 = 0;
		zeroCount = 0;
		if (bAutoBlankingEnabled) {	// initial state: silence
			ditherScaleFactor = 0.0;
		}
		else {	// initial state: dithering
			ditherScaleFactor = maxDitherScaleFactor;
		}
	}

// The Dither function ///////////////////////////////////////////////////////
//
// Ditherer Topology:
//                              Noise
//							     |
//                               v
//                    preDither [G]
//                         ^     |   +----------> preQuantize
//                         |     v   |               
//   inSample ----->+( )---+--->(+)--+->[Q]-->--+------> postQuantize
//                    -    |                    |
//                    ^    +---------->-( )+<---+
//                    |                  |               
//                 [filter]              | 
//                    |                  v
//                    +---<---[z^-1]-----+
//

FloatType Dither(FloatType inSample) {

	// Auto-Blanking
	if (bAutoBlankingEnabled) {
		if (std::abs(inSample) < autoBlankLevelThreshold) {
			++zeroCount;
			if (zeroCount > autoBlankTimeThreshold) {
				ditherScaleFactor *= autoBlankDecayFactor; // decay
				if (ditherScaleFactor < autoBlankDecayCutoff)
					ditherScaleFactor = 0.0; // decay cutoff
			}
		}
		else {
			zeroCount = 0; // reset
			ditherScaleFactor = maxDitherScaleFactor; // restore
		}
	} // ends auto-blanking

	FloatType tpdfNoise = (this->*noiseGenerator)() * ditherScaleFactor;
	
#ifdef TEST_FILTER
	//return (this->*noiseShapingFilter)(tpdfNoise); // (Output Only Filtered Noise - discard signal)
	FloatType preDither = - (this->*noiseShapingFilter)(Z1);
#else
	FloatType preDither = inSample - (this->*noiseShapingFilter)(Z1);
#endif
	FloatType preQuantize, postQuantize;
	preQuantize = preDither + tpdfNoise;
	postQuantize = reciprocalSignalMagnitude * round(maxSignalMagnitude * preQuantize); // quantize
	Z1 = (postQuantize - preDither);		
	return postQuantize;
} // ends function: Dither()

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

private:
	int oldRandom;
	int seed;
	FloatType Z1;				// last Quantization error
	FloatType maxSignalMagnitude;	// maximum integral value for signal target bit depth (for quantizing) 
	FloatType reciprocalSignalMagnitude; // for normalizing quantized signal back to +/- 1.0 
	FloatType maxDitherScaleFactor, ditherScaleFactor;	// maximum integral value for dither target bit depth
	int64_t zeroCount; // number of consecutive zeroes in input;
	FloatType autoBlankDecayCutoff;	// threshold at which ditherScaleFactor is set to zero during active blanking
	std::mt19937 randGenerator; // Mersenne Twister - one of the best random number algorithms available
	std::uniform_int_distribution<int> dist; // random number distribution
	static const int randMax = 16777215; // 2^24 - 1 */
	unsigned int signalBits;
	FloatType ditherBits;
	DitherProfile selectedDitherProfile;
	FloatType gain;
	bool bUseErrorFeedback;
	FloatType outputLimit;
	FloatType(Ditherer::*noiseShapingFilter)(FloatType); // function pointer to noise-shaping filter
	FloatType(Ditherer::*noiseGenerator)(); // function pointer to noise-generator
	bool bPulseEmitted;

	// Auto-Blanking parameters:
	bool bAutoBlankingEnabled;
	FloatType autoBlankLevelThreshold;				// input signals below this threshold are considered zero
	FloatType autoBlankTimeThreshold;				// number of zero samples before activating blanking
	const FloatType autoBlankDecayFactor = (FloatType)0.9995;	// dither level will decrease by this factor for each sample when blanking is active
	
	// IIR Filter-related stuff:
	Biquad<double> f1;
	Biquad<double> f2;
	Biquad<double> f3;
	Biquad<double> f4;

	// FIR Filter-related stuff:
	int FIRLength;
	FloatType FIRCoeffs[MAX_FIR_FILTER_SIZE];
	FloatType FIRHistory[MAX_FIR_FILTER_SIZE]; // (circular) buffer for noise history

	// --- Noise-generating functions ---

	// pure flat tpdf generator
	FloatType noiseGeneratorFlatTPDF() {
		int a = dist(randGenerator);
		int b = dist(randGenerator);
		return static_cast<FloatType>(a - b);
	}

	// the sloped TPDF generator re-uses the previous random number, giving it a "memory", 
	// which is equivalent to applying a single-pole HPF with 6dB / octave response
	// the slope is useful for providing more high-frequency emphasis (in addition to noise-shaping filter)
	FloatType noiseGeneratorSlopedTPDF() {
		int newRandom = dist(randGenerator);
		FloatType tpdfNoise = static_cast<FloatType>(newRandom - oldRandom);
		oldRandom = newRandom;
		return tpdfNoise;
	}

	FloatType noiseGeneratorRPDF() { // rectangular PDF (single PRNG)
		static constexpr int halfRand = (randMax + 1) >> 1;
		return static_cast<FloatType>(halfRand - dist(randGenerator));
	}

	FloatType noiseGeneratorGPDF() { // Gaussian PDF (n PRNGs)
		static constexpr int halfRand = (randMax + 1) >> 1;
		const int n = 3;
		FloatType r = 0;
		for (int i = 0; i < n; ++i) {
			r += dist(randGenerator);
		}
		
		return static_cast<FloatType>(halfRand - r/n);
	}

	FloatType noiseGeneratorImpulse() { // impulse - emits a single pulse at the begininng, followed by zeroes (for testing only)
		
		if (!bPulseEmitted) {
			bPulseEmitted = true;
			return static_cast<FloatType>(randMax);
		}
		
		return 0.0;
	}

	// --- Noise-shaping functions ---

	FloatType noiseShaperPassThrough(FloatType x) {
		return x;
	}

	FloatType noiseShaperCascadedBiquad(FloatType x) {
		return f3.filter(f2.filter(f1.filter(x)));
	}

	FloatType noiseShaperFIR(FloatType x) { // very simple FIR ...

		// put sample at end of buffer:
		FloatType* historyPtr = &FIRHistory[FIRLength - 1];
		*historyPtr = x;
		
		FloatType filterOutput = 0.0;
		
		// macc with coefficients:
		for (size_t k = 0; k < FIRLength; k++) {
			filterOutput += *historyPtr-- * FIRCoeffs[k];
		}

		// shift buffer backwards for next time:
		memmove(FIRHistory, &FIRHistory[1],
			(FIRLength - 1) * sizeof(FloatType));

		return filterOutput;
	}
};

// *Minimally Audible Noise Shaping
// STANLEY P. LIPSHITZ,JOHN VANDERKOOY, ROBERT A. WANNAMAKER
// J.AudioEng.Soc.,Vol.39,No.11,1991November

// **Psychoacoustically Optimal Noise Shaping
// Robert. A. Wannamaker
// Journal of the Audio Engineering Society 40(7 / 8) : 611 - 620 · July 1992

#endif // !DITHERER_H