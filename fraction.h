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

std::vector<Fraction> decomposeFraction(Fraction f, int maxStages) {
	std::vector<int> numerators;
	std::vector<int> denominators;

	// to-do: improve algorithm and remove these special cases:
	if (maxStages == 3 && f.numerator == 80) {
		numerators = { 2, 4, 10 };
	}
	else if (maxStages == 3 && f.numerator == 160) {
		numerators = { 2, 8, 10 };
	}
	else if (maxStages == 3 && f.numerator == 320) {
		numerators = { 4, 8, 10 };
	}
	else if (maxStages == 3 && f.numerator == 640) {
		numerators = { 4, 8, 20 };
	}
	else {
		numerators = factorize(f.numerator);
	}
	//
	if (maxStages == 3 && f.denominator == 80) {
		denominators = { 4, 4, 5 };
	}
	else if (maxStages == 3 && f.denominator == 160) {
		denominators = { 2, 8, 10 }; // still slower than single-stage
	}
	else if (maxStages == 3 && f.denominator == 320) {
		denominators = { 4, 5, 16 };
	}
	else if (maxStages == 3 && f.denominator == 640) {
		denominators = { 5, 8, 16 };
	}
	else if (maxStages == 3 && f.denominator == 64) {
		denominators = { 4, 4, 4 };
	}
	else {
		denominators = factorize(f.denominator);
	}

	// if too many items, consolidate into maxStages items - to-do: algorithm is very crude and produces suboptimal results - fix !
	while (numerators.size() > maxStages) {
		numerators[maxStages - 1] *= numerators.back();
		numerators.pop_back();
	}
	while (denominators.size() > maxStages) {
		denominators[maxStages - 1] *= denominators.back();
		denominators.pop_back();
	}

	int n = std::max(numerators.size(), denominators.size());
	assert(n <= maxStages);

	// if not enough items, fill with 1's at the front
	if (numerators.size() < n) {
		numerators.insert(numerators.begin(), n - numerators.size(), 1);
	}
	if (denominators.size() < n) {
		denominators.insert(denominators.begin(), n - denominators.size(), 1);
	}

	assert(numerators.size() == denominators.size());

	// pack numerators and denominators into vector of fractions
	std::vector<Fraction> fractions;
	for (int i = 0; i < n; i++) {
		Fraction f;
		f.numerator = numerators[i];
		f.denominator = denominators[i];
		fractions.push_back(f);
	}

	return fractions;
}


#endif // FRACTION_H
