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

struct Fraction {
	int numerator;
	int denominator;
};

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

//  getSimplifiedFraction() - reduces a fraction to its simplest form:
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
	
	std::vector<Fraction> fractions; // return value
	if (maxStages <= 1) {
		fractions.push_back(f);
		return fractions;
	}
		
	std::vector<int> numerators = factorize(f.numerator);
	std::vector<int> denominators = factorize(f.denominator);

	int n = std::max(numerators.size(), denominators.size());
	assert(n <= maxStages);

	// if not enough items, fill numerators with 1's at the end
	if (numerators.size() < n) {
	//	numerators.insert(numerators.begin(), n - numerators.size(), 1); // push front
		numerators.push_back(1);
	}

	// if not enough items, file denominators with 1's at the end
	if (denominators.size() < n) {
		//denominators.insert(denominators.begin(), n - denominators.size(), 1); // push front
		denominators.push_back(1);
	}

	assert(numerators.size() == denominators.size());

	// if too many items, consolidate into maxStages items - to-do: algorithm is very crude and produces suboptimal results - fix !
	while (numerators.size() > maxStages) {
		numerators[maxStages - 1] *= numerators.back();
		numerators.pop_back();
	}
	while (denominators.size() > maxStages) {
		denominators[maxStages - 1] *= denominators.back();
		denominators.pop_back();
	}

	// pack numerators and denominators into vector of fractions
	
	for (int i = 0; i < n; i++) {
		Fraction f;
		f.numerator = numerators[i];
		f.denominator = denominators[i];
		fractions.push_back(f);
	}

	return fractions;
}

//if (maxStages == 3 && f.numerator == 80) {
//	numerators = { 2, 4, 10 };
//}
//else if (maxStages == 3 && f.numerator == 160) {
//	numerators = { 2, 8, 10 };
//}
//else if (maxStages == 3 && f.numerator == 320) {
//	numerators = { 4, 8, 10 };
//}
//else if (maxStages == 3 && f.numerator == 640) {
//	numerators = { 4, 8, 20 };
//}
//else {
//	numerators = factorize(f.numerator);
//}
////
//if (maxStages == 3 && f.denominator == 80) {
//	denominators = { 4, 4, 5 };
//}
//else if (maxStages == 3 && f.denominator == 160) {
//	denominators = { 2, 8, 10 }; // still slower than single-stage
//}
//else if (maxStages == 3 && f.denominator == 320) {
//	denominators = { 4, 5, 16 };
//}
//else if (maxStages == 3 && f.denominator == 640) {
//	denominators = { 5, 8, 16 };
//}
//else if (maxStages == 3 && f.denominator == 64) {
//	denominators = { 4, 4, 4 };
//}
//else {
//	denominators = factorize(f.denominator);
//}

std::vector<Fraction> getPresetFractions(Fraction f, int maxStages) {

	if (maxStages <= 1) {
		std::vector<Fraction> fractions;
		fractions.push_back(f);
		return fractions;
	}

	struct PresetFractionSet {
		Fraction master;
		std::vector<Fraction> components;
	}; 
	
	// hardcoded table of known presets
	const std::vector<PresetFractionSet> presetList {
		{{147,320}, {{3,4}, {7,5}, {7,16}}},
		{{147,640 }, {{3,5},{7,8},{7,16}}},
	};
	
	// search for f in table
	for (auto& preset : presetList) {
		if (preset.master.numerator == f.numerator && preset.master.denominator == f.denominator) {
			return preset.components;
		}
	}

	// unknown fraction
	return decomposeFraction(f, maxStages); // decompose algorithmically
}


#endif // FRACTION_H
