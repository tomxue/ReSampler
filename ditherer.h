// Ditherer.h
// defines Ditherer class, for adding tpdf dither to input samples
//
// signalBits is the number of bits of the target bitformat
// ditherBits is the number of bits of dither to add, and doesn't have to be an integer
// input and output samples are of type FloatType
//

#ifndef DITHERER_H
#define DITHERER_H 1

 #define USE_IIR // if defined, use IIR Filter for noise shaping, otherwise use FIR. 
// Note: through experimentation, it has proven difficult to get the desired response curve using FIR
// that is why IIR was substituted ...

#include <cmath>
#include "biquad.h"

#ifndef USE_IIR
	#define NOISE_SHAPING_FILTER_SIZE 11
#endif // !USE_IIR

template<typename FloatType>
class Ditherer {
public:
	unsigned int signalBits;
	FloatType ditherBits;
	
	Ditherer(unsigned int signalBits, FloatType ditherBits)
		: signalBits(signalBits), ditherBits(ditherBits), E(0)
		
	{

#ifdef USE_IIR
		// IIR-related stuff:
		if (ditherBits < 1.5)
		{
			// IIR noise-shaping filter (2 biquads) - flatter response; more energy in spectrum
			f1.setCoeffs(0.798141839881378, -0.7040563852194521, 0.15341541599754416, 0.3060312586301247, 0.02511886431509577);
			f2.setCoeffs(0.5, -0.7215722413008345, 0.23235922079486643, -1.5531272249269004, 0.7943282347242815);
		}
		else
		{	
			// IIR noise-shaping filter (2 biquads)
			f1.setCoeffs(0.1872346691747817, -0.1651633303505913, 0.03598944852318585, 1.2861600144545022, 0.49000000000000016);
			f2.setCoeffs(0.5, -0.7215722413008345, 0.23235922079486643, -1.2511963408503206, 0.5328999999999999);
		}
#else 
		// FIR-related stuff:
		currentIndex = NOISE_SHAPING_FILTER_SIZE - 1;
		memset(noise, 0, NOISE_SHAPING_FILTER_SIZE * sizeof(FloatType));
#endif // USE_IIR
		
		signalMagnitude = static_cast<FloatType>(1 << (signalBits - 1));
		reciprocalSignalMagnitude = 1.0 / signalMagnitude;
		ditherMagnitude = pow(2,ditherBits-1) / signalMagnitude / RAND_MAX;
		oldRandom = newRandom = 0;
		
	}
	
	FloatType Dither(FloatType inSample) {
		FloatType shapedNoise = 0.0;
		newRandom = rand(); // to-do: better quality randoms ?
		inSample -= E; // apply Error Feedback
		FloatType tpdfNoise = ditherMagnitude * static_cast<FloatType>(newRandom - oldRandom);

#ifndef USE_IIR
		// FIR Noise Shaping:
		// put tpdf noise into history buffer (noise goes in "backwards"):
		noise[currentIndex] = tpdfNoise;
		currentIndex = currentIndex ? currentIndex - 1 : NOISE_SHAPING_FILTER_SIZE - 1;

		// Filter the noise:

		int index = currentIndex;
		for (int i = 0; i < NOISE_SHAPING_FILTER_SIZE; ++i) {
			shapedNoise += noise[index] * NoiseShapeCoeffs[i];
			index = (index == NOISE_SHAPING_FILTER_SIZE - 1) ? 0 : index + 1;
		}
#else
		// IIR Noise Shaping:
		shapedNoise = 1.05 * // tweak factor
			f2.filter(f1.filter(tpdfNoise)); // filter the triangular noise with two cascaded biQuads 
#endif

		outSample = inSample + shapedNoise;
		
		// Calculate the quantization error. This needs to exactly model the behavior of the I/O library writing samples to outfile, 
		// which in our case, appears to use round() ( ... as opposed to floor(), or cast-to-integer ...)
		
		quantizedOutSample = round(signalMagnitude * outSample); // quantize
		E = reciprocalSignalMagnitude * quantizedOutSample - inSample; // calculate error 
		
		oldRandom = newRandom;
		return outSample;
	}

private:
	int oldRandom, newRandom;
	FloatType E;				// last Quantization error
	int quantizedOutSample;
	FloatType outSample;
	FloatType signalMagnitude;	// maximum integral value for signal target bit depth (for quantizing) 
	FloatType reciprocalSignalMagnitude; // for normalizing quantized signal back to +/- 1.0 
	FloatType ditherMagnitude;	// maximum integral value for dither target bit depth
	
#ifdef USE_IIR
	// IIR Filters:
	Biquad<double> f1;
	Biquad<double> f2;
#else	
	// FIR Filter-related stuff:
	int currentIndex;
	// small FIR filter for noise-shaping:
	double NoiseShapeCoeffs[NOISE_SHAPING_FILTER_SIZE] = {
	
		// 11-tap:
		-0.014702063883960252,
		-0.0010319367876352055,
		0.06696663418581869,
		0.0010013618187379699,
		-0.30273448956854543,
		0.49822670684302905,
		-0.30273448956854543,
		0.0010013618187379699,
		0.06696663418581869,
		-0.0010319367876351854,
		-0.014702063883960252

		// super-duper 9-tap (Psychoacoustically Optimal Noise Shaping):
		/*2.412,
		-3.370,
		3.937,
		-4.174,
		3.353,
		-2.205,
		1.281,
		-0.569,
		0.0847*/



		// previous filter attempts: 

		// 9-tap:
	/*	-0.044563530870540866,
		-0.0512547777405835,
		0.01658357145118836,
		-0.205550523954609,
		0.5489106896304642,
		-0.20555052395460902,
		0.016583571451188384,
		-0.0512547777405835,
		-0.04456353087054085*/

		// 7-tap:
	/*	-0.021149859750950312,
		0.033161869936279724,
		-0.15368643534442833,
		0.5731512044421722,
		-0.15368643534442833,
		0.033161869936279724,
		-0.021149859750950312*/
	
	};

	// (circular) buffer for noise history:
	FloatType noise[NOISE_SHAPING_FILTER_SIZE];
#endif
};

#endif // !DITHERER_H