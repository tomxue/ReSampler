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

#define NOISE_SHAPING_FILTER_SIZE 21

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
		-0.0046834126581082255,
		0.0007635507792379351,
		0.01017176600030326,
		-0.016897478734828913,
		0.0007449345856288483,
		0.0357406681188413,
		-0.052843819844849316,
		0.00007713549280345215,
		0.12912288457582494,
		-0.2711134166929588,
		0.3331163915773193,
		-0.2711134166929588,
		0.12912288457582494,
		0.00007713549280344157,
		-0.0528438198448493,
		0.0357406681188413,
		0.0007449345856288483,
		-0.016897478734828913,
		0.01017176600030326,
		0.0007635507792379351,
		-0.0046834126581082255
	};

	// (circular) buffer for tpdf noise history:
	FloatType noise[NOISE_SHAPING_FILTER_SIZE];
};

#endif // !DITHERER_H