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
    MpxDecoder(int sampleRate)
    {
        // todo: emplace filters

        auto f1 = make19KhzBandpass<double>(sampleRate);
        auto f2 = make19KhzBandpass<double>(sampleRate);
        auto f3 = make57KhzBandpass<double>(sampleRate);
        std::vector<double> f0(f1.size(), 0);
        f0[(f1.size() - 1) / 2] = 1.0; // single impulse at halfway point
        filters.emplace_back(f0.data(), f0.size());
        filters.emplace_back(f1.data(), f1.size());
        filters.emplace_back(f2.data(), f2.size());
        filters.emplace_back(f3.data(), f3.size());
    }

    template <typename FloatType>
    static std::vector<FloatType> makeBandpass(int sampleRate, double ft1, double ft2)
    {
        // determine cutoff frequency and steepness
        double nyquist = sampleRate / 2.0;
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

    // 19khz bandpass filter for the Pilot Tone
    template<typename FloatType>
    static std::vector<FloatType> make19KhzBandpass(int sampleRate)
    {
        return makeBandpass<FloatType>(sampleRate, 18500, 19500);
    }

    // 38khz bandpass filter for the Audio Subcarrier
    template<typename FloatType>
    static std::vector<FloatType> make38KhzBandpass(int sampleRate)
    {
        return makeBandpass<FloatType>(sampleRate, 22000, 54000);
    }

    // 57khz bandpass filter for RDS / RBDS
    template<typename FloatType>
    static std::vector<FloatType> make57KhzBandpass(int sampleRate)
    {
        return makeBandpass<FloatType>(sampleRate, 54000, 60000);
    }

public:
    // function for testing / tweaking the performance of the bandpass filters
    static void saveFilters1(const std::string& filename)
    {
        std::vector<double> filt1 = make19KhzBandpass<double>(192000);
        std::vector<double> filt2 = make38KhzBandpass<double>(192000);
        SndfileHandle sndfile(filename, SFM_WRITE, SF_FORMAT_WAV | SF_FORMAT_FLOAT, 2, 192000);
        std::cout << "filter size " << filt1.size() << std::endl;
        std::vector<double> interleaved;
        interleaved.reserve(2 * filt1.size());
        for(int i = 0; i < filt1.size(); i++) {
            interleaved.push_back(filt1.at(i));
            interleaved.push_back(filt2.at(i));
        }
        sndfile.writef(interleaved.data(), filt1.size());
    }

private:
    std::vector<ReSampler::FIRFilter<double>> filters;

};

#endif // MPXDECODE_H
