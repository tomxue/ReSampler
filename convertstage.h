#ifndef CONVERTSTAGE_H
#define CONVERTSTAGE_H 1

#include "FIRFilter.h"

template<typename FloatType>
class ConvertStage
{
public:
    ConvertStage(int L, int M, FIRFilter<FloatType>& filter)
        : L(L), M(M), filter(filter)
    {
        // to-do: check that FIRFilter has copy constructor
        if(L == 1 && M == 1) {
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

    }
};

#endif