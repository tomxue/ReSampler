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

#include <utility>
#include <vector>
#include <set>
#include <numeric>
#include <iostream>

// fraction.h
// defines Fraction type, and functions for obtaining gcd, simplified fractions, and prime factors of integers

// single-stage policies
static const bool singleStageOnDecimateOnly = false;
static const bool singleStageOnInterpolateOnly = false;
// ---

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

//  getFractionFromSamplerates() - reduces an input and output sample rate to a simplified faction
// eg: 96000, 44100 => 147 / 320
Fraction getFractionFromSamplerates(int inputRate, int outputRate) 
{
	Fraction f;
	f.numerator = (outputRate / gcd(inputRate, outputRate));		// L (eg 147)
	f.denominator = (inputRate / gcd(inputRate, outputRate));		// M (eg 320)
	return f;
}

// factorize(n) - returns a vector of prime factors of n
std::vector<int> factorize(int n) {
	std::vector<int> factors;
	int maxFactor = static_cast<int>(std::sqrt(n));

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

// getnFactors() - take a vector of prime factors, and consolidate into
// a set of vectors of factors (each with size <= maxFactors) representing possible solutions
std::set<std::vector<int>> getnFactors(const std::vector<int> &primes, int maxFactors) {

    std::set<std::vector<int>> solutions; // the retval
    std::vector<int> currentFactors(maxFactors,1);

    std::function<void(std::vector<int>, int)> recursiveFunc =
    [&solutions, &currentFactors, &recursiveFunc](std::vector<int> primeFactors, int numFactors) {
        if(numFactors == 1) { // leaf node
            currentFactors[0] = std::accumulate(primeFactors.begin(), primeFactors.end(), 1, std::multiplies<int>());
            std::vector<int> newFactors = currentFactors;
            std::sort(newFactors.begin(), newFactors.end() , std::less<int>());
            solutions.insert(newFactors);
            return;
        }

        int maxFirstItems = static_cast<int>(primeFactors.size() - (numFactors - 1));
        for(int j = 1; j <= maxFirstItems; j++) {
            currentFactors[numFactors-1] = std::accumulate(primeFactors.begin(), primeFactors.begin() + j, 1, std::multiplies<int>());
            std::vector<int> remainingItems(primeFactors.begin() + j, primeFactors.end());
            recursiveFunc(remainingItems, numFactors - 1);
        }
        return;
    }; // ends recursiveFunc

    recursiveFunc(std::move(primes), maxFactors);

    return solutions;
}

// getnFactors() - factorize a number x into maxFactors (or less) factors.
// return a set of vectors representing possible solutions
std::set<std::vector<int>> getnFactors(int x, int maxFactors) {
	std::vector<int> primes = factorize(x);
	return getnFactors(primes, maxFactors);
}

// getDecompositionCandidates() : returns a vector of groups of fractions, with
// each group representing a possible decomposition of the input fraction into <= maxStages stages
// may also return an empty set if suitable solution is not possible with given value of maxStages

std::vector<std::vector<Fraction>> getDecompositionCandidates(Fraction f, int maxStages) {
	auto numeratorPrimes = factorize(f.numerator);
	auto denominatorPrimes = factorize(f.denominator);
	int maxPossibleStages = (int) std::max(numeratorPrimes.size(), denominatorPrimes.size()); // determines just how many stages can be formed
	int numStages = std::max(1, std::min(maxStages, maxPossibleStages)); // determines exactly how many stages we will have 
	
	while (numeratorPrimes.size() < numStages) { // pad with 1s at front
		numeratorPrimes.insert(numeratorPrimes.begin(), numStages - numeratorPrimes.size(), 1);
	}

	while (denominatorPrimes.size() < numStages) { // pad with 1s at front
		denominatorPrimes.insert(denominatorPrimes.begin(), numStages - denominatorPrimes.size(), 1);
	}

	auto numeratorGroups = getnFactors(numeratorPrimes, numStages);
	auto denominatorGroups = getnFactors(denominatorPrimes, numStages);
	std::vector<Fraction> tempFractionGroup(numStages, Fraction{ 1,1 });
	std::vector<std::vector<Fraction>> decompositionCandidates; // the return value
	double minRatio = std::min(1.0, static_cast<double>(f.numerator) / f.denominator); // to be a viable candidate, conversion ratio must be >= this value at all stages

	for (auto& numeratorGroup : numeratorGroups) {
		for (auto& denominatorGroup : denominatorGroups) {
			double ratio = 1.0;
			bool badRatio = false;

			for (int stage = 0; stage < numStages; stage++) {
				tempFractionGroup.at(stage) = Fraction{ numeratorGroup.at(stage), denominatorGroup.at(stage) };
				ratio *= static_cast<double>(tempFractionGroup.at(stage).numerator) / tempFractionGroup.at(stage).denominator;
				if (ratio < minRatio && stage != numStages - 1) {
					badRatio = true;
					break;
				}
			}

			if (!badRatio) {
				decompositionCandidates.push_back(tempFractionGroup);	
			}

		} // ends loop over denominatorGroups
	} // ends loop over numeratorGroups
	
	return decompositionCandidates;
}

std::vector<Fraction> decomposeFraction(Fraction f, int maxStages) {

	std::vector<Fraction> fractions; // return value
	if (maxStages <= 1) {
		fractions.push_back(f);
		return fractions;
	}

	std::vector<std::vector<Fraction>> solutions;

	// large values of maxStages may not produce a solution. 
	// Therefore, keep decreasing maxStages until solution is obtained
	// (solution is guaranteed for maxStages <= 1)

	do { 
		solutions = getDecompositionCandidates(f, maxStages);
		maxStages--;
	} while (solutions.empty() && maxStages > 0);
		
	return *solutions.rbegin(); // last is best
}

std::vector<Fraction> getPresetFractions(Fraction f, int maxStages) {

	// apply single-stage policies:
	if (maxStages <= 1) {
		return std::vector<Fraction> {f}; // single-stage conversion
	}
	
	if (f.numerator == 1 && singleStageOnDecimateOnly) {
		return std::vector<Fraction> {f}; // single-stage conversion
	}
	
	if (f.denominator == 1 && singleStageOnInterpolateOnly) {
		return std::vector<Fraction> {f}; // single-stage conversion
	}
	// ---

	struct PresetFractionSet {
		Fraction master;
		std::vector<Fraction> components;
	}; 
	
	// hardcoded table of known presets
	const std::vector<PresetFractionSet> presetList{

		{ { 5,147 },{ { 1,3 },{ 1,7 },{ 5,7 } } },
		{ { 147,40 },{ { 3,2 },{ 7,2 },{ 7,10 } } },
		{ { 147,80 },{ { 3,2 },{ 7,4 },{ 7,10 } } }, // cleanest 3-stage
		{{147,160},{{3,2},{7,8},{7,10}}},
		{{147,320}, {{3,5}, {7,8}, {7,8}}},
		{{147,640 }, {{3,5},{7,8},{7,16}}},
		{{4,1},{{4,1}}}
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

// utility functions:
void dumpFractionList(std::vector<Fraction> fractions) {
	for(auto fractionIt = fractions.begin(); fractionIt != fractions.end(); fractionIt++) {
		std::cout << fractionIt->numerator << "/" << fractionIt->denominator;
		if (fractionIt != std::prev(fractions.end())) {
			std::cout << " , ";
		}
	}
}

void dumpDecompositionCandidates(std::vector<std::vector<Fraction>> candidates) {
	for (std::vector<Fraction>& candidate : candidates) {
		dumpFractionList(candidate);
		std::cout << "\n";
	}
}

// test functions:
void testDecomposition(int numStages, bool unique = true) {
	std::vector<int> rates{8000, 11025, 16000, 22050, 32000, 37800, 44056, 44100, 47250, 48000, 50000, 50400, 88200, 96000, 176400, 192000, 352800, 384000, 2822400, 5644800};
	struct Result {
		Fraction fraction;
		std::vector<Fraction> fractionList;
	};

	std::vector<Result> decompositionList; 
	for (int i : rates) {
		for (int o : rates) {
			Result d;
			d.fraction = getFractionFromSamplerates(i, o);
			d.fractionList = decomposeFraction(d.fraction, numStages);
			decompositionList.push_back(d);
		}
	}

	if (unique) {
		struct Cmp { // comparison function object
			bool operator()(const Result& lhs, const Result& rhs) const {
				double r1 = static_cast<double>(lhs.fraction.numerator) / lhs.fraction.denominator;
				double r2 = static_cast<double>(rhs.fraction.numerator) / rhs.fraction.denominator; 
				return r1 < r2; 
			}
		};

		std::set<Result, Cmp> u(decompositionList.begin(), decompositionList.end()); // sort & de-dupe
		decompositionList.assign(u.begin(), u.end()); // convert back to vector again
	}

	for(auto& d : decompositionList) {
		double r = static_cast<double>(d.fraction.numerator) / d.fraction.denominator;
		std::cout << r << " , " << d.fraction.numerator << " , " << d.fraction.denominator << " , " << "\"";
		dumpFractionList(d.fractionList);
		std::cout << "\"\n";
	}

	std::cout << std::endl;
}

#endif // FRACTION_H
