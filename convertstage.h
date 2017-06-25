#ifndef CONVERTSTAGE_H
#define CONVERTSTAGE_H 1

#include "FIRFilter.h"

template<typename FloatType>
class ConvertStage
{
public:
    ConvertStage(int L, int M, FIRFilter<FloatType>& filter, bool passThrough = false)
        : L(L), M(M), filter(filter), l(0), m(0)
    {
        // to-do: check that FIRFilter has copy constructor
        if(passThrough) {
            convertFn = &ConvertStage::passThrough;
        }
        else if(L == 1 && M == 1) {
            convertFn = &ConvertStage::filterOnly;
        }
        else if (L != 1 && M == 1) {
            convertFn = &ConvertStage::interpolate;
        }
        else if (L == 1 && M != 1) {
            convertFn = &ConvertStage::decimate;
        }
        else {
            convertFn = &ConvertStage::interpolateAndDecimate;
        }
    }

    void convert(FloatType* outBuffer, size_t& outBufferSize, const FloatType* inBuffer, const size_t& inBufferSize) {
        (*convertFn)(outBuffer, outBufferSize, inBuffer, inBufferSize);
    }    

private:
    int L;
    int M;
    FIRFilter<FloatType> filter;
    int l;
    int m;

    void (*convertFn)(FloatType* outBuffer, size_t& outBufferSize, const FloatType* inBuffer, const size_t& inBufferSize);
    
    void passThrough(FloatType* outBuffer, size_t& outBufferSize, const FloatType* inBuffer, const size_t& inBufferSize) {
        memcpy(outBuffer, inBuffer, inBufferSize);
        outBufferSize = inBufferSize;
    }
    void filterOnly(FloatType* outBuffer, size_t& outBufferSize, const FloatType* inBuffer, const size_t& inBufferSize) {

    }
    void interpolate(FloatType* outBuffer, size_t& outBufferSize, const FloatType* inBuffer, const size_t& inBufferSize) {

    }
    void decimate(FloatType* outBuffer, size_t& outBufferSize, const FloatType* inBuffer, const size_t& inBufferSize) {

    }
	void interpolateAndDecimate(FloatType* outBuffer, size_t& outBufferSize, const FloatType* inBuffer, const size_t& inBufferSize) {
		sitze_t o = 0;
		for (size_t s = 0; s < inBufferSize; ++s) {
			fot(l = 0; l < L; ++i) {
				//(ii == 0) ? filter.put(s) : filter.putZero();
				filter.put((ii == 0) ? s : 0);
				if (m == 0) {
					outBuffer[o++] = filter.lazyGet(M);
				}
				if (m++ == L) {
					m = 0;
				}
			}
		}
	}
};

#endif