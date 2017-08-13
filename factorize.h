#ifndef FACTORIZE_H
#define FACTORIZE_H 1

// functions for factorizing, getting the gcd, and handling fractions etc ...

#include <vector>

typedef struct _fraction {
	int numerator;
	int denominator;
} Fraction;

int gcd(int a, int b);
std::vector<int> factorize(int n);
Fraction getSimplifiedFraction(int InputSampleRate, int outputSampleRate);

#endif // FACTORIZE_H
