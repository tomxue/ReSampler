// Ditherer.h
//
// (C) J.Niemann, 2016, 2017
//
// defines Ditherer class, for adding tpdf dither to input samples
//
//

#ifndef DITHERER_H
#define DITHERER_H 1

// configuration:
#define DITHER_TOPOLOGY 1
//#define TEST_FILTER // if defined, this will result in ditherer outputing the tpdf noise only. (Used for evaluating filters.)
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
	GPDF
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

// filters based on F-weighted curves
// from 'Psychoacoustically Optimal Noise Shaping' (*)
// this filter is the "F-Weighted" noise filter described by Wannamaker
// It is designed to produce minimum audibility:

const double wan3[] = {
	0.916014598760224, -0.554236805904214,  0.061519156663502
};

const double wan9[] = { // f-weighted used in SoX
	0.481319388145356, -0.672490189904582,  0.785636165476065,
	-0.832929985953034,  0.669097806157289, -0.440012127222434,
	0.255626092957795, -0.113545079541753,  0.016902053140925
};

const double wan24[] = {
	0.514957949278981, -0.707231224942256,  0.792298950922098,
	-0.782725058134344,  0.543526529766033, -0.246916297818599,
	0.024838892281917,  0.110623443620278, -0.161339968204987,
	0.110330813502457, -0.040696257820323, -0.009410889845009,
	0.032265323579584, -0.032554508456871,  0.01642992144958 ,
	-0.00259900332752 , -0.00454922479706 ,  0.005433144321457,
	-0.003471295165116,  0.000958853506002,  0.000188626919214,
	-0.000387374232494,  0.000166663510812, -0.000027561924269
};

// filters based on E-weighted curves
// from 'Minimally Audible Noise Shaping' (**)

const double modew44[] = { // Modified E-weighted (appendix: 2)
	0.866358574385245, -0.658369963567126,  0.251619304365678,
	-0.151847324138641,  0.066097633713628, -0.058591277834477,
	0.016951853693747, -0.006594125129948, -0.018369720915365
};

const double lips44[] = { // improved E-weighted (appendix: 5)
	0.690593301813357, -0.735432611129325,  0.665456113257436,
	-0.54010986221507 ,  0.208876449230218
};

const double impew44[] = { // improved E-weighted 9 coeff (appendix: 6)
	0.357422249248531, -0.588171140754959,  0.780127101099534,
	-0.901904263646452,  0.8334830743804  , -0.631734723645454,
	0.409648331330508, -0.204886937398526,  0.052615266828261
};

const double highShib44[20] = { // High-Shibata 44k (20 taps)
	0.210994777646055, -0.420248680433543,  0.64115984143207 ,
	-0.824542344864141,  0.890242066951182, -0.83102838215377 ,
	0.639689484122861, -0.374531480587619,  0.07944678322847 ,
	0.170730305215106, -0.346692245607667,  0.421108345040065,
	-0.41390894073698 ,  0.341901464927132, -0.247729870518492,
	0.15277447013028 , -0.081390587810871,  0.034194581138503,
	-0.011519110890133,  0.001618961712954
};

const double experimental1[] = {
	0.028745003564836, -0.075090008890252,  0.193254173536814,
	-0.404962386520082,  0.647462529770936, -0.817964064209564,
	0.789543428217277, -0.543123896114409,  0.127792565113188,
	0.288920282305045, -0.586693858174364,  0.677365569962496,
	-0.600226811847569,  0.421995454718511, -0.243835392469087,
	0.111111065840466, -0.03997384353328 ,  0.010306705271183,
	-0.002058234039559,  0.000220560246701
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

};

DitherProfile ditherProfileList[] = {

	{flat, "flat tpdf", flatTPDF, bypass, 44100, 1, noiseShaperPassThrough, false},
	{sloped,"sloped tpdf", slopedTPDF, bypass, 44100, 1, noiseShaperPassThrough, true },
	{standard, "standard", slopedTPDF, cascadedBiquad, 44100, 1, noiseShaperPassThrough, true},
	{Wannamaker3tap, "Wannamaker 3-tap",flatTPDF, fir, 44100, 3, wan3, true},
	{Wannamaker9tap, "Wannamaker 9-tap",flatTPDF, fir, 44100, 9, wan9, true},
	{Wannamaker24tap, "Wannamaker 24-tap",flatTPDF, fir, 44100, 24, wan24, true},
	{HighShibata44k, "High Shibata 44k",flatTPDF, fir, 44100, 20, highShib44, true},
	{ModEWeighted44k, "Modified E-Weighted",flatTPDF, fir, 44100, 9, modew44, true},
	{Lipshitz44k, "Lipshitz",flatTPDF, fir, 44100, 5, lips44, true},
	{ImpEWeighted44k, "Improved E-Weighted",flatTPDF, fir, 44100, 9, impew44, true},
	{Experimental1, "Experimental 1",flatTPDF, fir, 44100, 20, experimental1, true },
	{Experimental2, "Experimental 2",flatTPDF, fir, 44100, 12, experimental2, true },
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
		bUseErrorFeedback(ditherProfileList[ditherProfileID].bUseFeedback)
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
		currentIndex = FIRLength - 1;
		memset(noise, 0, MAX_FIR_FILTER_SIZE * sizeof(FloatType));

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
		currentIndex = FIRLength - 1;
		memset(noise, 0, MAX_FIR_FILTER_SIZE * sizeof(FloatType));
		
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

#if (DITHER_TOPOLOGY == 1)

// The Dither function ///////////////////////////////////////////////////////

// 1. Ditherer Topology:
//
//          tpdfNoise --[G]-->[filter]
//                               |
//                    preDither  |    +-----> preQuantize
//                        ^      v    |
//   inSample ----->+( )--+---->(+)---+--->[Q]-->-+--> postQuantize
//                    -   |                       |
//                    ^   +-------------->-( )+<--+
//                    |                     |
//                    +-------[z^-1]---<----+
//                      1.00               
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

	FloatType tpdfNoise = (this->*noiseGenerator)();
	FloatType preDither = bUseErrorFeedback ? inSample - Z1 : inSample;
	FloatType shapedNoise = (this->*noiseShapingFilter)(tpdfNoise) * ditherScaleFactor;

#ifdef TEST_FILTER
	return shapedNoise; // (Output Only Filtered Noise - discard signal)
#endif

	FloatType preQuantize = preDither + shapedNoise;
	FloatType postQuantize = reciprocalSignalMagnitude * round(maxSignalMagnitude * preQuantize); // quantize
	
	Z1 = postQuantize - preDither; // calculate error 
	return postQuantize;
} // ends function: Dither()

#endif //(DITHER_TOPOLOGY == 1)

#if (DITHER_TOPOLOGY == 2)

// The Dither function  - Topology #2 ///////////////////////////////////////////////////////

// 2. Ditherer Topology:
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
	return (this->*noiseShapingFilter)(tpdfNoise); // (Output Only Filtered Noise - discard signal)
#endif

	FloatType preDither = inSample - (this->*noiseShapingFilter)(Z1);
	FloatType preQuantize, postQuantize;
	preQuantize = preDither + tpdfNoise;
	postQuantize = reciprocalSignalMagnitude * round(maxSignalMagnitude * preQuantize); // quantize
	Z1 = (postQuantize - preDither);		
	return postQuantize;
} // ends function: Dither()

#endif //(DITHER_TOPOLOGY == 2)

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

	// Auto-Blanking parameters:
	bool bAutoBlankingEnabled;
	FloatType autoBlankLevelThreshold;				// input signals below this threshold are considered zero
	FloatType autoBlankTimeThreshold;				// number of zero samples before activating blanking
	const FloatType autoBlankDecayFactor = (FloatType)0.9995;	// dither level will decrease by this factor for each sample when blanking is active
	
	// IIR Filter-related stuff:
	Biquad<double> f1;
	Biquad<double> f2;

	// FIR Filter-related stuff:
	int currentIndex;
	FloatType FIRCoeffs[MAX_FIR_FILTER_SIZE];
	int FIRLength;
	FloatType noise[MAX_FIR_FILTER_SIZE]; // (circular) buffer for noise history

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

	// --- Noise-shaping functions ---

	FloatType noiseShaperPassThrough(FloatType x) {
		return x;
	}

	FloatType noiseShaperCascadedBiquad(FloatType x) {
		return f2.filter(f1.filter(x));
	}

	FloatType noiseShaperFIR(FloatType x) {
		// put x into history buffer (goes in "backwards"):
		noise[currentIndex--] = x;
		if (currentIndex < 0) {
			currentIndex = FIRLength - 1;
		}

		// get result from FIR:
		FloatType filterOutput = 0.0;
		int index = currentIndex;
		for (int i = 0; i < FIRLength; ++i) {
			if (++index == FIRLength) {
				index = 0;
			}
			filterOutput += noise[index] * FIRCoeffs[i];
		}
		return filterOutput;
	}
};

// *Psychoacoustically Optimal Noise Shaping
// Robert. A. Wannamaker
// Journal of the Audio Engineering Society 40(7 / 8) : 611 - 620 · July 1992

// **Minimally Audible Noise Shaping
// STANLEY P. LIPSHITZ,JOHN VANDERKOOY, ROBERT A. WANNAMAKER
// J.AudioEng.Soc.,Vol.39,No.11,1991November

#endif // !DITHERER_H