#ifndef CONVERTSTAGE_H
#define CONVERTSTAGE_H 1

#include "FIRFilter.h"

template<typename FloatType>
class ConvertStage
{
public:
    ConvertStage(int L, int M, FIRFilter<FloatType>& filter)
        : L(L), M(M), filter(filter), l(0), m(0), bypassMode(false)
    {
		SetConvertFunction();
    }

    void convert(FloatType* outBuffer, size_t& outBufferSize, const FloatType* inBuffer, const size_t& inBufferSize) {
        (this->*convertFn)(outBuffer, outBufferSize, inBuffer, inBufferSize);
    }

	void setBypassMode(bool bypassMode) {
		ConvertStage::bypassMode = bypassMode;
		SetConvertFunction();
	}

private:
    int L;	// interpoLation factor
    int M;	// deciMation factor
    FIRFilter<FloatType> filter;
    int l;	// interpolation index
    int m;	// decimation index
	bool bypassMode;
    
	// The following typedef defines the type 'ConvertFunction' which is a pointer to any of the member functions which 
	// take the arguments (FloatType* outBuffer, size_t& outBufferSize, const FloatType* inBuffer, const size_t& inBufferSize) ...
	typedef void (ConvertStage::*ConvertFunction) (FloatType* outBuffer, size_t& outBufferSize, const FloatType* inBuffer, const size_t& inBufferSize); // see https://isocpp.org/wiki/faq/pointers-to-members
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
            for(l = 0; l < L; ++l) {
				filter.put((l == 0) ? inBuffer[i] : 0);
				outBuffer[o++] = filter.get();
			}
        }
        outBufferSize = o;   
    }

	// decimate() - decimate and apply filter
    void decimate(FloatType* outBuffer, size_t& outBufferSize, const FloatType* inBuffer, const size_t& inBufferSize) {
        size_t o = 0;
        for (size_t i = 0; i < inBufferSize; ++i) {
			filter.put(inBuffer[i]);
            if (m == 0) {
                outBuffer[o++] = filter.get();    
            }
            if(++m == M) {
                m = 0;
            }
        }
        outBufferSize = o;
    }
    
	// interpolateAndDecimate()
	void interpolateAndDecimate(FloatType* outBuffer, size_t& outBufferSize, const FloatType* inBuffer, const size_t& inBufferSize) {
		size_t o = 0;
		for (size_t i = 0; i < inBufferSize; ++i) {
			for(l = 0; l < L; ++l) {
			//	filter.put((l == 0) ? inBuffer[i] : 0);
				((l == 0) ? filter.put(inBuffer[i]) : filter.putZero());
				if (m == 0) {
					outBuffer[o++] = filter.LazyGet(L);
					//outBuffer[o++] = filter.get();
				}
				if (++m == M) {
					m = 0;
				}
			}
		}
        outBufferSize = o;
	}

	void SetConvertFunction() {
		if (bypassMode) {
			convertFn = &ConvertStage::passThrough;
			std::cout << "convert mode: bypass" << std::endl;
		}
		else if (L == 1 && M == 1) {
			convertFn = &ConvertStage::filterOnly;
			std::cout << "convert mode: filter only" << std::endl;
		}
		else if (L != 1 && M == 1) {
			convertFn = &ConvertStage::interpolate;
			std::cout << "convert mode: interpolate" << std::endl;
		}
		else if (L == 1 && M != 1) {
			convertFn = &ConvertStage::decimate;
			std::cout << "convert mode: decimate" << std::endl;
		}
		else {
			convertFn = &ConvertStage::interpolateAndDecimate;
			std::cout << "convert mode: interpolate & decimate" << std::endl;
		}
	}
};

#endif