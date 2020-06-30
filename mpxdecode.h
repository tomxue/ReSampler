#ifndef MPXDECODE_H
#define MPXDECODE_H

#include <vector>

#include "FIRFilter.h"

class MpxDecoder
{

    // make 19khz bandpass filter for the Pilot Tone
    template <typename FloatType>
    static std::vector<FloatType> make19KhzBandpass(int sampleRate)
    {
        // determine cutoff frequency and steepness
        double nyquist = sampleRate / 2.0;
        double ft1 = 18500.0 / nyquist;
        double ft2 = 19500.0 / nyquist;
        double steepness = 0.090909091 / (1000.0 / nyquist);

        // determine filtersize
        int filterSize = static_cast<int>(
            std::min<int>(FILTERSIZE_BASE * steepness, FILTERSIZE_LIMIT)
            | 1 // ensure that filter length is always odd
        );

        // determine sidelobe attenuation
        int sidelobeAtten = 160;

        // Make some filter coefficients:

        // lower transition
        std::vector<FloatType> filterTaps1(filterSize, 0);
        FloatType* pFilterTaps1 = &filterTaps1[0];
        ReSampler::makeLPF<FloatType>(pFilterTaps1, filterSize, ft1, sampleRate);

        // upper transition
        std::vector<FloatType> filterTaps2(filterSize, 0);
        FloatType* pFilterTaps2 = &filterTaps2[0];
        ReSampler::makeLPF<FloatType>(pFilterTaps2, filterSize, ft2, sampleRate);

        // make bandpass
        for(int i = 0; i < filterSize; i++) {
            filterTaps2[i] -= filterTaps1.at(i);
        }

        ReSampler::applyKaiserWindow<FloatType>(pFilterTaps2, filterSize, ReSampler::calcKaiserBeta(sidelobeAtten));
        return filterTaps2;
    }
};

#endif // MPXDECODE_H
