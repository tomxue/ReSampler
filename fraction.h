/*
* Copyright (C) 2016 - 2017 Judd Niemann - All Rights Reserved
* You may use, distribute and modify this code under the
* terms of the GNU Lesser General Public License, version 2.1
*
* You should have received a copy of GNU Lesser General Public License v2.1
* with this file. If not, please refer to: https://github.com/jniemann66/ReSampler
*/

#ifndef FRACTION_H
#define FRACTION_H

#include <vector>

// fraction.h
// defines Fraction type, and functions for obtaining gcd, simplified fractions, and prime factors of integers

typedef struct fraction {
	int numerator;
	int denominator;
} Fraction;

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

//  getSimplifiedFraction() - turns a ratio into a fraction:
// eg: 96000, 44100 => 147 / 320
Fraction getSimplifiedFraction(int numerator, int denominator) 
{
	Fraction f;
	f.numerator = (denominator / gcd(numerator, denominator));		// L (eg 160)
	f.denominator = (numerator / gcd(numerator, denominator));		// M (eg 147)
	return f;
}

// factorize(n) - returns a vector of prime factors of n
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

#endif // FRACTION_H
