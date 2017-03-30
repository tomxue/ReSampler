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
	impulse,
	legacyTPDF
} NoiseGeneratorType;

typedef enum {
	flat,
	legacy,
	flat_f,
	ModEWeighted44k,
	Lipshitz44k,
	standard,
	standard88,
	standard96,
	standard176,
	standard192,
	Wannamaker24tap,
	Wannamaker9tap,
	smooth,
	slick,
	ImpEWeighted44k,
	HighShibata44k,
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
	
	{ flat, "flat tpdf", flatTPDF, bypass, 44100, 1, noiseShaperPassThrough, false },
	{ legacy, "classic", legacyTPDF, bypass, 44100, 10, noiseShaperPassThrough, true },
	{ flat_f, "flat tpdf (with error-correction feedback)", flatTPDF, fir, 44100, 1, noiseShaperPassThrough, true },
	{ ModEWeighted44k, "Modified E-Weighted",flatTPDF, fir, 44100, 9, modew44, true },
	{ Lipshitz44k, "Lipshitz",flatTPDF, fir, 44100, 5, lips44, true },
	{ standard, "standard", slopedTPDF, fir, 44100, 10, std_44, true },
	{ standard88, "standard (88k)", slopedTPDF, fir, 44100, 12, standard_88, true },
	{ standard96, "standard (96k)", slopedTPDF, fir, 44100, 12, standard_96, true },
	{ standard176, "standard (176k)", slopedTPDF, fir, 44100, 10, standard_176, true },
	{ standard192, "standard (192k)", slopedTPDF, fir, 44100, 10, standard_192, true },
	{ Wannamaker24tap, "Wannamaker 24-tap",flatTPDF, fir, 44100, 24, wan24, true },
	{ Wannamaker9tap, "Wannamaker 9-tap",flatTPDF, fir, 44100, 9, wan9, true },
	{ smooth, "smooth", slopedTPDF, fir, 44100, 10, smooth_44, true },
	{ slick, "slick", slopedTPDF, fir, 44100, 10, notch12250_2_44, true },
	{ ImpEWeighted44k, "Improved E-Weighted",flatTPDF, fir, 44100, 9, impew44, true },
	{ HighShibata44k, "High Shibata 44k",flatTPDF, fir, 44100, 20, highShib44, true },
	{ Experimental1, "Experimental 1",slopedTPDF, fir, 44100, 12, standard_88, true },
	{ Experimental2, "Experimental 2",slopedTPDF, fir, 44100, 6, standard_96, true },
	{ rpdf,"flat rectangular pdf", RPDF, bypass, 44100, 1, noiseShaperPassThrough, false }
	

	/*
	{ flat, "flat tpdf", impulse, bypass, 44100, 1, noiseShaperPassThrough, false },
	{ legacy, "classic", impulse, bypass, 44100, 10, noiseShaperPassThrough, true },
	{ Wannamaker3tap, "Wannamaker 3-tap",impulse, fir, 44100, 3, wan3, true },
	{ ModEWeighted44k, "Modified E-Weighted",impulse, fir, 44100, 9, modew44, true },
	{ Lipshitz44k, "Lipshitz",impulse, fir, 44100, 5, lips44, true },
	{ standard, "standard", impulse, fir, 44100, 10, std_44, true },
	{ standard88, "standard (88k)", impulse, fir, 44100, 12, standard_88, true },
	{ standard96, "standard (96k)", impulse, fir, 44100, 12, standard_96, true },
	{ standard176, "standard (176k)", impulse, fir, 44100, 10, standard_176, true },
	{ standard192, "standard (192k)", impulse, fir, 44100, 10, standard_192, true },
	{ Wannamaker24tap, "Wannamaker 24-tap",impulse, fir, 44100, 24, wan24, true },
	{ Wannamaker9tap, "Wannamaker 9-tap",impulse, fir, 44100, 9, wan9, true },
	{ smooth, "smooth", impulse, fir, 44100, 10, smooth_44, true },
	{ slick, "slick", impulse, fir, 44100, 10, notch12250_2_44, true },
	{ ImpEWeighted44k, "Improved E-Weighted",impulse, fir, 44100, 9, impew44, true },
	{ HighShibata44k, "High Shibata 44k",impulse, fir, 44100, 20, highShib44, true },
	{ Experimental1, "Experimental 1",impulse, fir, 44100, 12, standard_88, true },
	{ Experimental2, "Experimental 2",impulse, fir, 44100, 6, standard_96, true },
	{ rpdf,"flat rectangular pdf", impulse, bypass, 44100, 1, noiseShaperPassThrough, false }
	*/


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
		bPulseEmitted(false),
		masterVolume(1.0)
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
		case legacyTPDF:
			noiseGenerator = &Ditherer::noiseGeneratorLegacy;
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
		if (ditherBits < 1.5)
		{
			// IIR noise-shaping filter (2 biquads) - flatter response; more energy in spectrum
			f1.setCoeffs(0.798141839881378,
				-0.7040563852194521,
				0.15341541599754416,
				0.3060312586301247,
				0.02511886431509577);

			f2.setCoeffs(0.5,
				-0.7215722413008345,
				0.23235922079486643,
				-1.5531272249269004,
				0.7943282347242815);
		}
		else
		{
			// IIR noise-shaping filter (2 biquads)
			f1.setCoeffs(0.1872346691747817,
				-0.1651633303505913,
				0.03598944852318585,
				1.2861600144545022,
				0.49000000000000016);

			f2.setCoeffs(0.5,
				-0.7215722413008345,
				0.23235922079486643,
				-1.2511963408503206,
				0.5328999999999999);
		}

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
		autoBlankDecayCutoff = 0.25 * reciprocalSignalMagnitude / randMax;
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
		masterVolume = 1.0;
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
//                    preDither [G1]
//                         ^     |   +----------> preQuantize
//                         |     v   |               
//   inSample ----->+( )---+--->(+)--+--[G2]-->[Q]-->--+------> postQuantize
//                    -    |                           |
//                    ^    +---------->-( )+<----------+
//                    |                  |               
//                 [filter]              | 
//                    |                  v
//                    +---<---[z^-1]-----+
//
//  Gain Stages:
//	G1 = ditherScaleFactor
//  G2 = masterVolume
//

FloatType Dither(FloatType inSample) {

	// Auto-Blanking
	if (bAutoBlankingEnabled) {
		if (std::abs(inSample) < autoBlankLevelThreshold) {
			++zeroCount;
			if (zeroCount > autoBlankTimeThreshold) {
				ditherScaleFactor *= autoBlankDecayFactor; // decay
				if (ditherScaleFactor < autoBlankDecayCutoff) {
					ditherScaleFactor = 0.0; // decay cutoff
					masterVolume = 0.0; // mute
				}
			}
		}
		else {
			zeroCount = 0; // reset
			ditherScaleFactor = maxDitherScaleFactor; // restore
			masterVolume = 1.0;
		}
	} // ends auto-blanking

	FloatType tpdfNoise = (this->*noiseGenerator)() * ditherScaleFactor;
	FloatType preDither = bUseErrorFeedback ? inSample - (this->*noiseShapingFilter)(Z1) : inSample;
	FloatType preQuantize, postQuantize;
	preQuantize = masterVolume * (preDither + tpdfNoise);
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
	FloatType masterVolume;
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
		const int n = 5;
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

	FloatType noiseGeneratorLegacy() { // legacy noise generator (from previous version of ReSampler) - applies filter to noise _before_ injection into dither engine
		int newRandom = dist(randGenerator);
		FloatType tpdfNoise = static_cast<FloatType>(newRandom - oldRandom); // sloped TDPF
		oldRandom = newRandom;
		return f2.filter(f1.filter(tpdfNoise));
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