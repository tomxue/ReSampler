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
#include <set>
#include <numeric>

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

// getnFactors() - factorize a number x into numFactors factors.
// return a set of vectors representing possible solutions
std::set<std::vector<int>> getnFactors(int x, int numFactors) {

	std::vector<int> primes = factorize(x);
    std::set<std::vector<int>> results;
    std::vector<int> currentFactors(numFactors,1);

    std::function<void(std::vector<int>, int)> recursiveFunc =
    [&results, &currentFactors, &recursiveFunc](std::vector<int> primeFactors, int numFactors) {
        if(numFactors == 1) { // leaf node
            currentFactors[0] = std::accumulate(primeFactors.begin(), primeFactors.end(), 1, std::multiplies<int>());
            std::vector<int> newFactors = currentFactors;
            std::sort(newFactors.begin(), newFactors.end() , std::less<int>());
            results.insert(newFactors);
            return;
        }

        int maxFirstItems = primeFactors.size() - (numFactors - 1);
        for(int j = 1; j <= maxFirstItems; j++) {
            currentFactors[numFactors-1] = std::accumulate(primeFactors.begin(), primeFactors.begin() + j, 1, std::multiplies<int>());
            std::vector<int> remainingItems(primeFactors.begin() + j, primeFactors.end());
            recursiveFunc(remainingItems, numFactors - 1);
        }
        return;
    }; // ends recursiveFunc

    recursiveFunc(primes, numFactors);

    return results;
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

// getDecompositionCandidates() : returns a vector of groups of fractions, with
// each group representing a possible decomposition of the input fraction into maxStages stages

std::vector<std::vector<Fraction>> getDecompositionCandidates(Fraction f, int maxStages) {

	auto numeratorGroups = getnFactors(f.numerator, maxStages);
	auto denominatorGroups = getnFactors(f.denominator, maxStages);
	std::vector<Fraction> tempFractionGroup(maxStages, Fraction{ 1,1 });
	std::vector<std::vector<Fraction>> decompositionCandidates; // the return value
	double minRatio = std::min(1.0, static_cast<double>(f.numerator) / f.denominator);

	for (auto& numeratorGroup : numeratorGroups) {
		for (auto& denominatorGroup : denominatorGroups) {
			double ratio = 1.0;
			for (int stage = 0; stage < maxStages; stage++) {
				tempFractionGroup.at(stage) = Fraction{ numeratorGroup.at(stage), denominatorGroup.at(stage) };
				ratio *= static_cast<double>(tempFractionGroup.at(stage).numerator) / tempFractionGroup.at(stage).denominator;
				if (ratio >= minRatio) {
					decompositionCandidates.push_back(tempFractionGroup);
				}
			}
		}
	}

	return decompositionCandidates;
}


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
	const std::vector<PresetFractionSet> presetList{
		{ { 147,40 },{ { 3,2 },{ 7,4 },{ 7,5 } } },
		{{147,80},{{3,2},{7,5},{7,8}}},
		{{147,160},{{3,2},{7,8},{7,10}}},
		{{147,320}, {{3,5}, {7,8}, {7,8}}},
	//	{ { 147,320 },{ { 147,320 } }},
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
