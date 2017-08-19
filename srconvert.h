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
	AbstractResampler(const ConversionInfo& ci) : ci(ci) {}
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

#endif // SRCONVERT_H