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
#include "biquad.h"
#include <random>

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
	sloped,
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
	rpdf
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

//////////////////////////
//
// Filter Coefficients
//
//////////////////////////

const double noiseShaperPassThrough[1] = {
	1
};

// filters based on E-weighted curves
// from 'Minimally Audible Noise Shaping' (*)

const double modew44[] = { // Modified E-weighted (appendix: 2)
	1.6620, -1.2630, 0.4827,
	-0.2913, 0.1268,-0.1124,
	0.03252, -0.01265, -0.03524
};

const double lips44[] = { // improved E-weighted (appendix: 5)
	2.033, -2.165, 1.959,
	-1.590, 0.6149
};

const double impew44[] = { // improved E-weighted 9 coeff (appendix: 6)
	2.847, -4.685,  6.214,
	-7.184, 6.639, -5.032,
	3.263, -1.632, 0.4191
};

// filters based on F-weighted curves
// from 'Psychoacoustically Optimal Noise Shaping' (**)
// this filter is the "F-Weighted" noise filter described by Wannamaker
// It is designed to produce minimum audibility:

const double wan3[] = { // Table 3; 3 Coefficients
	1.623, -0.982, 0.109
};

const double wan9[] = { // Table 3; 9 Coefficients ('f-weighted' in SoX)
	2.4120002321781  , -3.370000324394779,  3.937000378973959,
	-4.174000401787478,  3.353000322758366, -2.205000212252369,
	1.281000123308519, -0.569000054771701,  0.084700008153185
};

const double wan24[] = { // Table 4; 24 Coefficients
	2.391510032751124, -3.284444044979632,  3.679506050389904,
	-3.635044049781009,  2.524185034568077, -1.146701015703782,
	0.115354001579743,  0.51374500703561 , -0.749277010261162,
	0.512386007016999, -0.188997002588268, -0.043705000598528,
	0.149843002052063, -0.151186002070453,  0.076302001044937,
	-0.012070000165296, -0.021127000289329,  0.025232000345547,
	-0.016121000220773,  0.004453000060982,  0.000876000011999,
	-0.001799000024635,  0.0007740000106  , -0.000128000001755
};

const double highShib44[20] = { // High-Shibata 44k (20 taps)
	
	3.0259189605712890625, -6.0268716812133789062,   9.195003509521484375,
	-11.824929237365722656, 12.767142295837402344, -11.917946815490722656,
	9.1739168167114257812,  -5.3712320327758789062, 1.1393624544143676758,
	2.4484779834747314453,  -4.9719839096069335938,   6.0392003059387207031,
	-5.9359521865844726562,  4.903278350830078125,   -3.5527443885803222656,
	2.1909697055816650391, -1.1672389507293701172,  0.4903914332389831543,
	-0.16519790887832641602,  0.023217858746647834778
	
};

const double experimental1[] = {
	1.552893042564392, -1.459592938423157,  0.557879626750946,
	0.004247912205756, -0.247446924448013, -0.010727410204709,
	0.024524994194508, -0.084559261798859, -0.021680492907763,
	-0.007851422764361, -0.010477361269295, -0.006129192188382,
	-0.005679284688085, -0.00112113065552 , -0.001433005789295,
	-0.000978846335784, -0.000214568979573,  0.000104520288005,
	0.001021763193421,  0.001436339691281
};

const double experimental2[] = {
	/*
	-0.156273992376491,
	0.288411025954531,
	-0.197304884034655,
	0.129645690377117,
	-0.070106447409447,
	-0.041085952214833,
	0.050192058358397,
	0.012898104533912,
	0.015359661780711,
	0.020719908866823,
	-0.011813660625704,
	0.007166224946099,
	-0.005038759518396,
	-0.006848775277087,
	0.004246829537086,
	0.003585617408335
	*/
	/*
	- 0.113267579111506,
	0.215432144212954,
	-0.206263340877997,
	0.246064185480709,
	-0.054824084366412,
	0.179481292975895,
	-0.002556243072638,
	0.013174521887528,
	-0.000238789304217,
	0.122047711532922,
	0.088340424932406,
	0.032665078075013
	*/

	3.528870344161987,  -7.397897720336914,  11.062552452087402,
	-12.915676116943359,  12.005429267883301,  -8.976534843444824,
	5.175451755523682,  -2.103446006774902,   0.348601877689362,
	0.211675748229027,  -0.179740116000175,   0.030822610482574,
	0.024646494537592,  -0.010657157748938,  -0.007913443259895,
	0.004375995136797,   0.004815561696887,  -0.00556468218565 ,
	0.001390703488141,   0.002223829738796,  -0.000689988140948,
	-0.000269168231171,   0.000755135435611,  -0.

};

DitherProfile ditherProfileList[] = {

	// id, name, noiseGeneratorType, filterType, intendedSampleRate, N, coeffs, bUseFeedback

	{flat, "flat tpdf", flatTPDF, bypass, 44100, 1, noiseShaperPassThrough, false},
	{sloped,"sloped tpdf", slopedTPDF, bypass, 44100, 1, noiseShaperPassThrough, true },
	{standard, "standard", slopedTPDF, cascadedBiquad, 44100, 1, noiseShaperPassThrough, true},
	{Wannamaker3tap, "Wannamaker 3-tap",flatTPDF, fir, 44100, 3, wan3, true},
	{Wannamaker9tap, "Wannamaker 9-tap",flatTPDF, fir, 44100, 9, wan9, true},
	{Wannamaker24tap, "Wannamaker 24-tap",flatTPDF, fir, 44100, 24, wan24, true},
	{HighShibata44k, "High Shibata 44k",slopedTPDF, fir, 44100, 20, highShib44, true},
	{ModEWeighted44k, "Modified E-Weighted",flatTPDF, fir, 44100, 9, modew44, true},
	{Lipshitz44k, "Lipshitz",flatTPDF, fir, 44100, 5, lips44, true},
	{ImpEWeighted44k, "Improved E-Weighted",flatTPDF, fir, 44100, 9, impew44, true},
	{Experimental1, "Experimental 1",flatTPDF, fir, 44100, 20, experimental1, true },
	{Experimental2, "Experimental 2",flatTPDF, fir, 44100, 20, experimental2, true },
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

		/*

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

		*/

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
//
//							 tpdfNoise
//                               |
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