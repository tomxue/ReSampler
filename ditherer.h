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

#define NOISE_SHAPING_FILTER_SIZE 19

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
		if (currentIndex == 0)
			currentIndex = NOISE_SHAPING_FILTER_SIZE - 1; // wrap
		else
			--currentIndex;

		// Filter the noise:
		FloatType shapedNoise = 0.0;
		int index = currentIndex;
		for (int i = 0; i <  NOISE_SHAPING_FILTER_SIZE; ++i) {
			shapedNoise += noise[index] * NoiseShapeCoeffs[i];
			if (index == NOISE_SHAPING_FILTER_SIZE -1)
				index = 0; // wrap
			else
				index++;
		}

		oldRandom = newRandom;
		return inSample + shapedNoise;
	}

private:
	int oldRandom, newRandom;
	FloatType signalMagnitude;
	FloatType ditherMagnitude;
	int currentIndex;
	
	// small FIR filter for noise-shaping:
	double NoiseShapeCoeffs[NOISE_SHAPING_FILTER_SIZE] = {
		
		0.00015310801577591,
		- 0.00058049470516376,
		0.00263188880890734,
		- 0.00210302799542606,
		- 0.00008889366607690,
		0.02589640750256570,
		- 0.05980358305724670,
		0.12852161421719700,
		- 0.16645794147960500,
		0.20276672371195100,
		- 0.16645794147960500,
		0.12852161421719700,
		- 0.05980358305724670,
		0.02589640750256570,
		- 0.00008889366607690,
		- 0.00210302799542606,
		0.00263188880890734,
		- 0.00058049470516376,
		0.00015310801577591,

	};

	// (circular) buffer for tpdf noise history:
	FloatType noise[NOISE_SHAPING_FILTER_SIZE];
};

#endif // !DITHERER_H