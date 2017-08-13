#include "factorize.h"

// gcd() - greatest common divisor:
int gcd(int a, int b) {
	if (a<0) a = -a;
	if (b<0) b = -b;
	while (b != 0) {
		a %= b;
		if (a == 0) return b;
		b %= a;
	}
	return a;
}

std::vector<int> factorize(int n) {
	std::vector<int> factors;
	int maxFactor = std::sqrt(n);

	for (int factor = 2; factor <= maxFactor; factor++) {
		while (n % factor == 0) {
			factors.push_back(factor);
			n /= factor;
		}
		if (n == 1)
			return factors;
	}

	factors.push_back(n);
	return factors;
}

//  getSimplifiedFraction() - turns a sample-rate ratio into a fraction:
Fraction getSimplifiedFraction(int inputSampleRate, int outputSampleRate)			// eg 44100, 48000
{
	Fraction f;
	f.numerator = (outputSampleRate / gcd(inputSampleRate, outputSampleRate));		// L (eg 160)
	f.denominator = (inputSampleRate / gcd(inputSampleRate, outputSampleRate));		// M (eg 147)
	return f;
}