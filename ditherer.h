// Ditherer.h
// defines Ditherer class, for adding tpdf dither to input samples
//
// signalBits is the number of bits of the target bitformat
// ditherBits is the number of bits of dither to add, and doesn't have to be an integer
// input and output samples are of type FloatType
//

#ifndef DITHERER_H
#define DITHERER_H 1

#include <cmath>

#define NOISE_SHAPING_FILTER_SIZE 11

template<typename FloatType>
class Ditherer {
public:
	unsigned int signalBits;
	FloatType ditherBits;
	
	Ditherer(unsigned int signalBits, FloatType ditherBits) 
		: signalBits(signalBits), ditherBits(ditherBits)
	{
		//
		signalMagnitude = static_cast<FloatType>(1 << (signalBits - 1));
		ditherMagnitude = pow(2,ditherBits-1) / signalMagnitude / RAND_MAX;
		oldRandom = newRandom = 0;
		currentIndex = NOISE_SHAPING_FILTER_SIZE - 1;
		memset(noise, 0, NOISE_SHAPING_FILTER_SIZE * sizeof(FloatType));
	}
	
	FloatType Dither(FloatType inSample) {
		newRandom = rand();
		
		// put tpdf noise into history buffer (noise goes in "backwards"):
		noise[currentIndex] = ditherMagnitude * static_cast<FloatType>(newRandom - oldRandom);
		currentIndex = currentIndex ? currentIndex-1 : NOISE_SHAPING_FILTER_SIZE - 1; 
	
		// Filter the noise:
		FloatType shapedNoise = 0.0;
		int index = currentIndex;
		for (int i = 0; i <  NOISE_SHAPING_FILTER_SIZE; ++i) {
			shapedNoise += noise[index] * NoiseShapeCoeffs[i];
			index = (index == NOISE_SHAPING_FILTER_SIZE - 1) ? 0 : index + 1;
		}

		oldRandom = newRandom;
		outSample = inSample + shapedNoise;
		
		return outSample;
	}

private:
	int oldRandom, newRandom;
		
	FloatType outSample;
	FloatType signalMagnitude;
	FloatType ditherMagnitude;
	int currentIndex;
	
	// to-do: feed quantize error signal into filter
	// to-do: FIR doesn't really cut it; need IIR to get the right curve !!
	
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

	// (circular) buffer for tpdf noise history:
	FloatType noise[NOISE_SHAPING_FILTER_SIZE];
};

#endif // !DITHERER_H