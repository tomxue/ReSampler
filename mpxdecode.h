/*
* Copyright (C) 2020 Judd Niemann - All Rights Reserved.
* You may use, distribute and modify this code under the
* terms of the GNU Lesser General Public License, version 2.1
*
* You should have received a copy of GNU Lesser General Public License v2.1
* with this file. If not, please refer to: https://github.com/jniemann66/ReSampler
*/

#ifndef MPXDECODE_H
#define MPXDECODE_H

// FM Broadcast Multiplex Decoder (WIP)

// todo: 19khz bandpass for Pilot Tone
// todo: 23-53khz bandpass for audio subcarrier
// todo: audio deemphasis filter (IIR ?)
// todo: 19khz -> 38khz frequency doubler
// todo: subcarrier multiplier / demodulator
// todo: mid / side matrix decoder
// todo: 15khz lpf for audio channels

#include <vector>
#include <sndfile.h>
#include <sndfile.hh>

#include "FIRFilter.h"

class MpxDecoder
{
    // make 19khz bandpass filter for the Pilot Tone
    template <typename FloatType>
    static std::vector<FloatType> make19KhzBandpass(int sampleRate)
    {
        // determine cutoff frequency and steepness
        double nyquist = sampleRate / 2.0;
        double ft1 = 18500.0;
        double ft2 = 19500.0;
        double steepness = 0.090909091 / (1000.0 / nyquist);

        // determine filtersize
        int filterSize = static_cast<int>(
            std::min<int>(FILTERSIZE_BASE * steepness, FILTERSIZE_LIMIT)
            | 1 // ensure that filter length is always odd
        );

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

        // determine sidelobe attenuation
        int sidelobeAtten = 60;
        ReSampler::applyKaiserWindow<FloatType>(pFilterTaps2, filterSize, ReSampler::calcKaiserBeta(sidelobeAtten));
        return filterTaps2;
    }

public:
    // function for testing / tweaking the performance of the 19khz bandpass
    static void save19khzFilter(const std::string& filename)
    {
        std::vector<double> filt = make19KhzBandpass<double>(192000);
        SndfileHandle sndfile(filename, SFM_WRITE, SF_FORMAT_WAV | SF_FORMAT_FLOAT, 1, 192000);
        std::cout << "19khz filter size " << filt.size() << std::endl;
        sndfile.write(filt.data(), filt.size());
    }
};

#endif // MPXDECODE_H
