// Ditherer.h
//
// by J.Niemann, 2016
//
// defines Ditherer class, for adding tpdf dither to input samples
//
// signalBits is the number of bits of the target bitformat
// ditherBits is the number of bits of dither to add, and doesn't have to be an integer
// input and output samples are of type FloatType
//

#ifndef DITHERER_H
#define DITHERER_H 1

#define FIR_NOISE_SHAPING_FILTER_SIZE 9

#define USE_IIR // if defined, use IIR Filter for noise shaping, otherwise use FIR. 
// Note: through experimentation, it has proven difficult to get the desired response curve using FIR
// that is why IIR was substituted ...

#define USE_HIGH_QUALITY_RANDOM // if defined, uses C++ std library, instead of rand(), for "better" random numbers.
// Note: the audible difference in quality between rand() and MT is VERY noticable, so high-quality random numbers are important.

#define DITHER_USE_SATURATION  // restrict output amplitude to +/- 0.999 (guards against excessive dither levels causing clipping)

#include <cmath>
#include "biquad.h"
#include <random>

template<typename FloatType>
class Ditherer {
public:
	unsigned int signalBits;
	FloatType ditherBits;

	// Auto-Blanking parameters:
	bool bAutoBlankingEnabled; 
	FloatType autoBlankLevelThreshold;				// input signals below this threshold are considered zero
	FloatType autoBlankTimeThreshold;				// number of zero samples before activating blanking
	const FloatType autoBlankDecayFactor = (FloatType)0.9995;	// dither level will decrease by this factor for each sample when blanking is active

// Constructor:

	Ditherer(unsigned int signalBits, FloatType ditherBits, bool bAutoBlankingEnabled, int seed)
		: signalBits(signalBits), ditherBits(ditherBits), bAutoBlankingEnabled(bAutoBlankingEnabled), seed(seed), E1(0), E2(0)

#ifdef USE_HIGH_QUALITY_RANDOM
		,randGenerator(seed)		// initialize (seed) RNG
		,dist(0,RAND_MAX)		// set the range of the random number distribution
#endif

	{

#ifdef USE_IIR
		//// IIR-related stuff:
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

		// super-slick HPF (1 biquad):
		// f1.setCoeffs(0.008978326844454853, -0.017956653688909707, 0.008978326844454853, 1.6865949922346122, 0.7224999999999999);
		
#else 
		// FIR-related stuff:
		currentIndex = FIR_NOISE_SHAPING_FILTER_SIZE - 1;
		memset(noise, 0, FIR_NOISE_SHAPING_FILTER_SIZE * sizeof(FloatType));

#endif // USE_IIR
		
		signalMagnitude = static_cast<FloatType>((1 << (signalBits - 1)) - 1); // note the -1 : match 32767 scaling factor for 16 bit !
		
		reciprocalSignalMagnitude = 1.0 / signalMagnitude;
		maxDitherMagnitude = pow(2,ditherBits-1) / signalMagnitude / RAND_MAX;
		oldRandom = newRandom = 0;

		if (bAutoBlankingEnabled) {	// initial state: silence
			ditherMagnitude = 0.0;
		}
		else {	// initial state: dithering
			ditherMagnitude = maxDitherMagnitude; 
		}

		autoBlankLevelThreshold = 1.0 / pow(2, 32); // 1 LSB of 32-bit digital
		autoBlankTimeThreshold = 30000; // number of zero samples before activating autoblank
		autoBlankDecayCutoff = 0.7 * reciprocalSignalMagnitude / RAND_MAX;
		zeroCount = 0;
		
	} // Ends Constructor 

// The Dither function ///////////////////////////////////////////////////////

	FloatType Dither(FloatType inSample) {

		// Auto-Blanking
		if (bAutoBlankingEnabled) {
			if (abs(inSample) < autoBlankLevelThreshold) {
				++zeroCount;
				if (zeroCount > autoBlankTimeThreshold) {
					ditherMagnitude *= autoBlankDecayFactor; // decay
					if (ditherMagnitude < autoBlankDecayCutoff)
						ditherMagnitude = 0.0; // decay cutoff
				}
			}
			else {
				zeroCount = 0; // reset
				ditherMagnitude = maxDitherMagnitude; // restore
			}
		} // ends auto-blanking

#ifdef USE_HIGH_QUALITY_RANDOM
		newRandom = dist(randGenerator);
#else
		newRandom = rand();
#endif	

FloatType tpdfNoise = ditherMagnitude * static_cast<FloatType>(newRandom - oldRandom);
FloatType shapedNoise = 0.0;


// 1. Ditherer Topology:
//
//          tpdf Dither------>[filter]
//                               |
//                               |    +-----> outSample
//                               v    |
//   inSample ----->+( )------->(+)---+--->[Q]-->-+--> quantizedOutSample
//                    +               |           |
//                    ^               +-->-( )+<--+
//                    |                     |
//                    +-------[z^-1]--------+
//                    | -1.00               |
//                    +-------[z^-2]--------+
//                      +0.043
//
//

		inSample -= E1;				// apply 1st order Error Feedback
		inSample += E2*0.043;		// apply 2nd order Error Feedback
	
#ifndef USE_IIR
		// FIR Noise Shaping:
		// put tpdf noise into history buffer (noise goes in "backwards"):
		noise[currentIndex] = tpdfNoise;
		currentIndex = currentIndex ? currentIndex - 1 : FIR_NOISE_SHAPING_FILTER_SIZE - 1;

		// Filter the noise:

		int index = currentIndex;
		for (int i = 0; i < FIR_NOISE_SHAPING_FILTER_SIZE; ++i) {
			shapedNoise += noise[index] * NoiseShapeCoeffs[i];
			index = (index == FIR_NOISE_SHAPING_FILTER_SIZE - 1) ? 0 : index + 1;
		}
#else
		// IIR Noise Shaping:
		shapedNoise = f2.filter(f1.filter(tpdfNoise)); // filter the triangular noise with two cascaded biQuads 
		//	shapedNoise = f1.filter(tpdfNoise); // filter the triangular noise with one biQuad 
#endif
		outSample = inSample + shapedNoise;
		
		// Calculate the quantization error. This needs to exactly model the behavior of the I/O library writing samples to outfile, 
		// otherwise nasty quantization distortion will result.
		
		quantizedOutSample = reciprocalSignalMagnitude * round(signalMagnitude * outSample); // quantize
			
		E2 = E1;
		E1 = quantizedOutSample - inSample; // calculate error 
		
		oldRandom = newRandom;
		
#ifdef DITHER_USE_SATURATION
		// branchless clipping - restrict to +/- 1.0 (0.0000 dB):
		quantizedOutSample = 0.5*(fabs(quantizedOutSample + 1.0) - fabs(quantizedOutSample - 1.0));
#endif
		return quantizedOutSample;
	} // ends function: Dither()

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

private:
	int oldRandom, newRandom;
	int seed;
	FloatType E1,E2;				// last Quantization error
	FloatType quantizedOutSample;
	FloatType outSample;
	FloatType signalMagnitude;	// maximum integral value for signal target bit depth (for quantizing) 
	FloatType reciprocalSignalMagnitude; // for normalizing quantized signal back to +/- 1.0 
	FloatType maxDitherMagnitude, ditherMagnitude;	// maximum integral value for dither target bit depth
	__int64 zeroCount; // number of consecutive zeroes in input;
	FloatType autoBlankDecayCutoff;	// threshold at which ditherMagnitude is set to zero during active blanking

#ifdef USE_HIGH_QUALITY_RANDOM
	std::mt19937 randGenerator; // Mersenne Twister - one of the best random number algorithms available
	std::uniform_int_distribution<int> dist; // random number distribution
#endif
	
#ifdef USE_IIR
	// IIR Filters:
	Biquad<double> f1;
	Biquad<double> f2;
#else	
	// FIR Filter-related stuff:
	int currentIndex;
	// small FIR filter for noise-shaping:
	double NoiseShapeCoeffs[FIR_NOISE_SHAPING_FILTER_SIZE] = {
	
		// 11-tap:
		/*-0.014702063883960252,
		-0.0010319367876352055,
		0.06696663418581869,
		0.0010013618187379699,
		-0.30273448956854543,
		0.49822670684302905,
		-0.30273448956854543,
		0.0010013618187379699,
		0.06696663418581869,
		-0.0010319367876351854,
		-0.014702063883960252*/

		// Psychoacoustically Optimal Noise Shaping:
		// this filter is the "F-Weighted" noise filter described by Wannamaker.
		// It is designed to produce minimum audibility, but I personally don't like the sound of it. 
		// However, YMMV ...
		// From experimentation, it seems clear that 
		// "audibility" and "pleastantness/unpleasantness" of noise are two very different things !

		// 3-tap:
	/*	1.623,
		-0.982,
		0.109*/

		// 9-tap:
			2.412,
		-3.370,
		3.937,
		-4.174,
		3.353,
		-2.205,
		1.281,
		-0.569,
		0.0847 
	
		//0.09648,
		//- 0.13480,
		//0.15748,
		//- 0.16696,
		//0.13412,
		//- 0.08820,
		//0.05124,
		//- 0.02276,
		//0.00339 // normalized to 0dB peak (x0.04)
			
		
		//0.00339,  // reversed, normalized to 0dB peak (x0.04)
		//- 0.02276,
		//0.05124,
		//- 0.08820,
		//0.13412,
		//- 0.16696,
		//0.15748,
		//- 0.13480,
		//0.09648

		//24-tap (needs normalization):
	/*	2.391510,
		-3.284444,
		3.679506,
		-3.635044,
		2.524185,
		-1.146701,
		0.115354,
		0.513745,
		-0.749277,
		0.512386,
		-0.188997,
		-0.043705,
		0.149843,
		-0.151186,
		0.076302,
		-0.012070,
		-0.021127,
		0.025232,
		-0.016121,
		0.004453,
		0.000876,
		-0.001799,
		0.000774,
		-0.000128*/
		
		// some previous filter attempts: 

		// 9-tap:
	/*	-0.044563530870540866,
		-0.0512547777405835,
		0.01658357145118836,
		-0.205550523954609,
		0.5489106896304642,
		-0.20555052395460902,
		0.016583571451188384,
		-0.0512547777405835,
		-0.04456353087054085 */

		// 7-tap:
	/*	-0.021149859750950312,
		0.033161869936279724,
		-0.15368643534442833,
		0.5731512044421722,
		-0.15368643534442833,
		0.033161869936279724,
		-0.021149859750950312 */
	
	};

	// (circular) buffer for noise history:
	FloatType noise[FIR_NOISE_SHAPING_FILTER_SIZE];
#endif
};

// *Psychoacoustically Optimal Noise Shaping
// Robert. A. Wannamaker
// Journal of the Audio Engineering Society 40(7 / 8) : 611 - 620 · July 1992

#endif // !DITHERER_H