/*
* Copyright (C) 2016 - 2017 Judd Niemann - All Rights Reserved
* You may use, distribute and modify this code under the
* terms of the GNU Lesser General Public License, version 2.1
*
* You should have received a copy of GNU Lesser General Public License v2.1
* with this file. If not, please refer to: https://github.com/jniemann66/ReSampler
*/

#ifndef SRCONVERT_H
#define SRCONVERT_H 1

// srconvert.h : core sample rate conversion code.

#include "FIRFilter.h"
#include "conversioninfo.h"
#include "fraction.h"

static_assert(std::is_copy_constructible<ConversionInfo>::value, "ConversionInfo needs to be copy Constructible");
static_assert(std::is_copy_assignable<ConversionInfo>::value, "ConversionInfo needs to be copy Assignable");

template<typename FloatType>
std::vector<FloatType> makeFilterCoefficients(const ConversionInfo& ci, Fraction fraction) {
	// determine base filter size
	int baseFilterSize;
	int overSamplingFactor = 1;
	Fraction f = fraction;
	if ((fraction.numerator != fraction.denominator) && (fraction.numerator <= 4 || fraction.denominator <= 4)) { // simple ratios
		baseFilterSize = FILTERSIZE_MEDIUM * std::max(fraction.denominator, fraction.numerator) / 2;
		if (ci.bMinPhase) { // oversample to improve filter performance
			overSamplingFactor = 8;
			f.numerator *= overSamplingFactor;
			f.denominator *= overSamplingFactor;
		}
	}
	else { // complex ratios
		baseFilterSize = FILTERSIZE_HUGE * std::max(fraction.denominator, fraction.numerator) / 320;
	}

	// determine cutoff frequency and steepness
	double targetNyquist = std::min(ci.inputSampleRate, ci.outputSampleRate) / 2.0;
	double ft = (ci.lpfCutoff / 100.0) * targetNyquist;
	double steepness = steepness = 0.090909091 / (ci.lpfTransitionWidth / 100.0);

	// scale the filter size, according to selected options:
	int filterSize = std::min(static_cast<int>(overSamplingFactor * baseFilterSize * steepness), FILTERSIZE_LIMIT)
		| static_cast<int>(1);	// ensure that filter length is always odd

	// determine sidelobe attenuation
	int sidelobeAtten = ((fraction.numerator == 1) || (fraction.denominator == 1)) ?
		195 :
		160;

	// Make some filter coefficients:
	int overSampFreq = ci.inputSampleRate * f.numerator;
	std::vector<FloatType> filterTaps(filterSize, 0);
	FloatType* pFilterTaps = &filterTaps[0];
	makeLPF<FloatType>(pFilterTaps, filterSize, ft, overSampFreq);
	applyKaiserWindow<FloatType>(pFilterTaps, filterSize, calcKaiserBeta(sidelobeAtten));

	// conditionally convert filter coefficients to minimum-phase:
	if (ci.bMinPhase) {
		makeMinPhase<FloatType>(pFilterTaps, filterSize);
	}

	return filterTaps;
}

template<typename FloatType>
class ResamplingStage
{
public:
    ResamplingStage(int L, int M, FIRFilter<FloatType>& filter, bool bypassMode = false)
        : L(L), M(M), filter(filter), bypassMode(bypassMode), m(0)
    {
		SetConvertFunction();
    }

    void convert(FloatType* outBuffer, size_t& outBufferSize, const FloatType* inBuffer, const size_t& inBufferSize) {
		(this->*convertFn)(outBuffer, outBufferSize, inBuffer, inBufferSize);
    }

	void setBypassMode(bool bypassMode) {
		ResamplingStage::bypassMode = bypassMode;
		SetConvertFunction();
	}

private:
    int L;	// interpoLation factor
    int M;	// deciMation factor
	int m;	// decimation index
    FIRFilter<FloatType> filter;
	bool bypassMode;
    
	// The following typedef defines the type 'ConvertFunction' which is a pointer to any of the member functions which 
	// take the arguments (FloatType* outBuffer, size_t& outBufferSize, const FloatType* inBuffer, const size_t& inBufferSize) ...
	typedef void (ResamplingStage::*ConvertFunction) (FloatType* outBuffer, size_t& outBufferSize, const FloatType* inBuffer, const size_t& inBufferSize); // see https://isocpp.org/wiki/faq/pointers-to-members
	ConvertFunction convertFn;

	// passThrough() - just copies input straight to output (used in bypassMode mode)
    void passThrough(FloatType* outBuffer, size_t& outBufferSize, const FloatType* inBuffer, const size_t& inBufferSize) {
        memcpy(outBuffer, inBuffer, inBufferSize * sizeof(FloatType));
        outBufferSize = inBufferSize;
    }

	// filerOnly() - keeps 1:1 conversion ratio, but applies filter
    void filterOnly(FloatType* outBuffer, size_t& outBufferSize, const FloatType* inBuffer, const size_t& inBufferSize) {
        size_t o = 0;
        for(size_t i = 0; i < inBufferSize; ++i) {
            filter.put(inBuffer[i]);
            outBuffer[o++] = filter.get();
        }
        outBufferSize = inBufferSize;
    }

	// interpolate() - interpolate (zero-stuffing) and apply filter:
    void interpolate(FloatType* outBuffer, size_t& outBufferSize, const FloatType* inBuffer, const size_t& inBufferSize) {
        size_t o = 0;
        for (size_t i = 0; i < inBufferSize; ++i) {
            for(int l = 0; l < L; ++l) {
				//filter.put((l == 0) ? inBuffer[i] : 0);
				((l == 0) ? filter.put(inBuffer[i]) : filter.putZero());
				//outBuffer[o++] = filter.get();
				outBuffer[o++] = filter.lazyGet(L);
			}
        }
        outBufferSize = o;   
    }

	// decimate() - decimate and apply filter
    void decimate(FloatType* outBuffer, size_t& outBufferSize, const FloatType* inBuffer, const size_t& inBufferSize) {
        size_t o = 0;
		int localm = m;
        for (size_t i = 0; i < inBufferSize; ++i) {
			filter.put(inBuffer[i]);
            if (localm == 0) {
                outBuffer[o++] = filter.get();    
            }
            if(++localm == M) {
                localm = 0;
            }
        }
        outBufferSize = o;
		m = localm;
    }
    
	// interpolateAndDecimate()
	void interpolateAndDecimate(FloatType* outBuffer, size_t& outBufferSize, const FloatType* inBuffer, const size_t& inBufferSize) {
		size_t o = 0;
		int localm = m;
		for (size_t i = 0; i < inBufferSize; ++i) {
			for(int l = 0; l < L; ++l) {
				((l == 0) ? filter.put(inBuffer[i]) : filter.putZero());
				if (localm == 0) {
					outBuffer[o++] = filter.lazyGet(L);
				}
				if (++localm == M) {
					localm = 0;
				}
			}
		}
        outBufferSize = o;
		m = localm;
	}

	void SetConvertFunction() {
		if (bypassMode) {
			convertFn = &ResamplingStage::passThrough;
		}
		else if (L == 1 && M == 1) {
			convertFn = &ResamplingStage::filterOnly;
		}
		else if (L != 1 && M == 1) {
			convertFn = &ResamplingStage::interpolate;
		}
		else if (L == 1 && M != 1) {
			convertFn = &ResamplingStage::decimate;
		}
		else {
			convertFn = &ResamplingStage::interpolateAndDecimate;
		}
	}
};

template <typename FloatType>
class AbstractResampler
{
public:
	virtual void convert(FloatType* outBuffer, size_t& outBufferSize, const FloatType* inBuffer, const size_t& inBufferSize) = 0;
	int getGroupDelay() {
		return groupDelay;
	}
protected:
	AbstractResampler(const ConversionInfo& ci) : ci(ci), groupDelay(0) {}
	ConversionInfo ci;
	int groupDelay;
	std::vector<ResamplingStage<FloatType>> convertStages;
};

template <typename FloatType>
class SingleStageResampler : public AbstractResampler<FloatType>
{
public:
	SingleStageResampler(const ConversionInfo& ci) : AbstractResampler<FloatType>(ci) {
		Fraction f = getSimplifiedFraction(ci.inputSampleRate, ci.outputSampleRate);
		std::vector<FloatType> filterTaps = makeFilterCoefficients<FloatType>(ci, f);
		bool bypassMode = (f.numerator == 1 && f.denominator == 1);

		int overSamplingFactor = ci.bMinPhase && (f.numerator != f.denominator) && (f.numerator <= 4 || f.denominator <= 4) ? 8 : 1;
		f.numerator *= overSamplingFactor;
		f.denominator *= overSamplingFactor;

		FIRFilter<FloatType> firFilter(filterTaps.data(), filterTaps.size());
		AbstractResampler<FloatType>::convertStages.emplace_back(f.numerator, f.denominator, firFilter, bypassMode);
		AbstractResampler<FloatType>::groupDelay = (ci.bMinPhase || !ci.bDelayTrim) ? 0 : (filterTaps.size() - 1) / 2 / f.denominator;
		if (f.numerator == 1 && f.denominator == 1) {
			AbstractResampler<FloatType>::groupDelay = 0;
		}
	}
	void convert(FloatType* outBuffer, size_t& outBufferSize, const FloatType* inBuffer, const size_t& inBufferSize) {
		AbstractResampler<FloatType>::convertStages[0].convert(outBuffer, outBufferSize, inBuffer, inBufferSize);
	}
};

template <typename FloatType>
class MultiStageResampler : public AbstractResampler<FloatType>
{
public:
	MultiStageResampler(const ConversionInfo& ci) : AbstractResampler<FloatType>(ci), numStages(3), indexOfLastStage(2) {
		makeConversionParams();
	}

	void convert(FloatType* outBuffer, size_t& outBufferSize, const FloatType* inBuffer, const size_t& inBufferSize) {
		for (int i = 0; i < numStages; i++) {
			// first stage reads directly from inBuffer. Subsequent stages read from output of previous stage
			const FloatType* in = (i == 0) ? inBuffer : intermediateOutputBuffers[i - 1].data(); 
			size_t inSize = (i == 0) ? inBufferSize : intermediateOutputBuffers[i - 1].size();
			FloatType* out = (i == indexOfLastStage) ? outBuffer : intermediateOutputBuffers[i].data(); // last stage writes straight to outBuffer;
			convertStages[i].convert(out, outBufferSize, in, inSize);
		}
	}

private:
	/*const std::vector<Fraction> ratios {
		{3, 4},
		{7, 8},
		{7, 10}
	};*/

	int numStages;
	int indexOfLastStage;
	std::vector<std::vector<FloatType>> intermediateOutputBuffers;	// intermediate output buffer for each ConvertStage;
	
	std::vector<Fraction> getPartialRatios(Fraction f, int maxStages) {
		std::vector<int> numerators;
		std::vector<int> denominators;
		
		// to-do: improve algorithm and remove these special cases:
		if (maxStages == 3  && f.numerator == 80) {
			numerators = { 2, 4, 10 };
		}
		else if (maxStages == 3 && f.numerator == 160) {
			numerators = { 4, 4, 10 };
		}
		else if (maxStages == 3 && f.numerator == 320) {
			numerators = { 4, 8, 10 };
		}
		else {
			numerators = factorize(f.numerator);
		}
		//
		if (maxStages == 3 && f.denominator == 80) {
			denominators = { 2, 4, 10 };
		}
		else if (maxStages == 3 && f.denominator == 160) {
			denominators = { 4, 4, 10 };
		}
		else if (maxStages == 3 && f.denominator == 320) {
			denominators = { 4, 8, 10 };
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

		int nStages = std::max(numerators.size(), denominators.size());
		assert(nStages <= maxStages);

		// if not enough items, fill with 1's at the front
		if (numerators.size() < nStages) {
			numerators.insert(numerators.begin(), nStages - numerators.size(), 1);
		}
		if (denominators.size() < nStages) {
			denominators.insert(denominators.begin(), nStages - denominators.size(), 1);
		}
	
		assert(numerators.size() == denominators.size());
		
		// pack numerators and denominators into vector of fractions
		std::vector<Fraction> fractions;
		for (int i = 0; i < numStages; i++) {
			Fraction f;
			f.numerator = numerators[i];
			f.denominator = denominators[i];

			std::cout << "numerator: " << f.numerator << std::endl;
			std::cout << "denominator: " << f.denominator << std::endl;
			fractions.push_back(f);
		}

		return fractions;
	}

	void makeConversionParams() {
		Fraction masterConversionRatio = getSimplifiedFraction(ci.inputSampleRate, ci.outputSampleRate);
		auto ratios = getPartialRatios(masterConversionRatio, 3);
		numStages = ratios.size();
		int inputRate = ci.inputSampleRate;
		double guarantee = inputRate / 2.0;
		double ft = ci.lpfCutoff/100 * ci.outputSampleRate / 2.0;
		
		for (int i = 0; i < numStages; i++) {
			ConversionInfo newCi = ci;
			newCi.inputSampleRate = inputRate;
			newCi.outputSampleRate = inputRate * ratios[i].numerator / ratios[i].denominator;
			decltype(newCi.inputSampleRate) minSampleRate = std::min(newCi.inputSampleRate, newCi.outputSampleRate);
			double stopFreq = std::max(minSampleRate / 2.0, minSampleRate - guarantee);
			double widthAdjust = ((stopFreq - ft) / (minSampleRate / 2.0 - ft));
			newCi.lpfTransitionWidth *= widthAdjust;
			guarantee = std::min(guarantee, stopFreq);
			
			// dump stage parameters:
			std::cout << "Stage: " << 1 + i << "\n";
			std::cout << "inputRate: " << newCi.inputSampleRate << "\n";
			std::cout << "outputRate: " << newCi.outputSampleRate << "\n";
			std::cout << "ft: " << ft << "\n";
			std::cout << "stopFreq: " << stopFreq << "\n";
			std::cout << "widthAdjust: " << widthAdjust << "\n";
			std::cout << "guarantee: " << guarantee << "\n";
			
			// make the filter coefficients
			std::vector<FloatType> filterTaps = makeFilterCoefficients<FloatType>(newCi, ratios[i]);
			std::cout << "Generated Filter Size: " << filterTaps.size() << "\n";
			
			// make the filter
			FIRFilter<FloatType> firFilter(filterTaps.data(), filterTaps.size());

			// make the ConvertStage:
			Fraction f = ratios[i];
			bool bypassMode = (f.numerator == 1 && f.denominator == 1);
			int overSamplingFactor = ci.bMinPhase && (f.numerator != f.denominator) && (f.numerator <= 4 || f.denominator <= 4) ? 8 : 1;
			f.numerator *= overSamplingFactor;
			f.denominator *= overSamplingFactor;
			convertStages.emplace_back(f.numerator, f.denominator, firFilter, bypassMode);
			
			// add Group Delay:
			groupDelay += (ci.bMinPhase || !ci.bDelayTrim) ? 0 : (filterTaps.size() - 1) / 2 / f.denominator;
			if (f.numerator == 1 && f.denominator == 1) {
				groupDelay = 0;
			}

			// calculate size of output buffer for this stage:
			double cumulativeNumerator = 1.0;
			double cumulativeDenominator = 1.0;
			for (int j = 0; j <= i; j++) {
				cumulativeNumerator *= ratios[j].numerator;
				cumulativeDenominator *= ratios[j].denominator;
			}
			size_t outBufferSize = std::ceil(BUFFERSIZE * cumulativeNumerator / cumulativeDenominator);
			std::cout << "Output Buffer Size: " << outBufferSize << "\n\n" << std::endl;

			// make the outputBuffer (last stage doesn't need one)
			if (i != numStages - 1) {	
				intermediateOutputBuffers.emplace_back(std::vector<FloatType>(outBufferSize, 0));
			}
			
			// set input rate of next stage
			inputRate = newCi.outputSampleRate; 
		}
	}
};

#endif // SRCONVERT_H