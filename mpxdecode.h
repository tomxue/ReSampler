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

#include <vector>
#include <sndfile.h>
#include <sndfile.hh>

#include "FIRFilter.h"

//#define MPXDECODER_TUNE_PILOT_AGC

class MpxDecoder
{
public:
    MpxDecoder(int sampleRate)
    {
        // create filters
        auto f1 = make19KhzBandpass<double>(sampleRate);
        auto f2 = make38KhzBandpass<double>(sampleRate);
        auto f3 = make57KhzBandpass<double>(sampleRate);
        auto lpf = make15khzLowpass<double>(sampleRate);
        std::vector<double> f0(f1.size(), 0);
        f0[(f1.size() - 1) / 2] = 1.0; // single impulse at halfway point

        length = f1.size();
                delayLine.resize(length, 0.0);
                centerTap = (length - 1) / 2;
                currentIndex = length - 1;

        filters.emplace_back(f0.data(), f0.size());
        filters.emplace_back(f1.data(), f1.size());
        filters.emplace_back(f2.data(), f2.size());
        filters.emplace_back(f3.data(), f3.size());
        filters.emplace_back(lpf.data(), lpf.size()); // left lowpass
        filters.emplace_back(lpf.data(), lpf.size()); // right lowpass

        // these values determined experimentally:
        // (#define MPXDECODER_TUNE_PILOT_AGC to debug & tweak)
        decreaseRate = std::pow(10.0, /* dB per sec = */ -3200.0 / sampleRate / 20.0);
        peakDecreaseRate = std::pow(10.0, -450.0 / sampleRate / 20.0);
        increaseRate = std::pow(10.0, 3200.0 / sampleRate / 20.0);
    }

#ifdef MPXDECODER_TUNE_PILOT_AGC
    ~MpxDecoder()
    {
        int64_t totalCount = plusCount + stableCount + minusCount;
        std::cout << 100.0 * plusCount / totalCount << "%, "
                  << 100.0 * stableCount / totalCount << "%, "
                  << 100.0 * minusCount / totalCount << "%" << std::endl;

        std::cout << "Peak Pilot Gain: " << peakPilotGain << std::endl;
    }
#endif

    template<typename FloatType>
    std::pair<FloatType, FloatType> decode(FloatType input)
    {

        //todo: clean up

#ifdef USE_FIR_AS_DELAY
        filters.at(0).put(input);
        FloatType monoRaw = filters.at(0).get();
#else
        delayLine[currentIndex] = input; // place input into history
        if(currentIndex == 0) {
            currentIndex = length - 1;
        } else {
            currentIndex--;
        }
        int d = currentIndex + centerTap;
        FloatType monoRaw = delayLine[d >= length ? d - length : d];
#endif

        // --- //

        filters.at(1).put(input);
        filters.at(2).put(input);

        FloatType pilotRaw = filters.at(1).get();
        FloatType sideRaw = filters.at(2).get();
        FloatType pilot = pilotRaw * pilotGain;
        FloatType pilotAbs = std::fabs(pilot);

        if(pilotAbs > pilotPeak) {
            pilotPeak = pilotAbs;
        }

#ifdef MPXDECODER_TUNE_PILOT_AGC
        std::cout << pilotGain << ", " << pilotPeak << "\n";
#endif

        if(pilotPeak < pilotStableLow) { // pilot too quiet
            pilotGain *= increaseRate;
            if(pilotGain >= pilotMaxGain) {
                pilotGain = pilotMaxGain;
            }

#ifdef MPXDECODER_TUNE_PILOT_AGC
            plusCount++;
#endif

        } else if(pilotPeak > pilotStableHigh) { // pilot too loud
            pilotGain *= decreaseRate;

#ifdef MPXDECODER_TUNE_PILOT_AGC
            minusCount++;
#endif

        } else { // stable pilot tone

#ifdef MPXDECODER_TUNE_PILOT_AGC
            stableCount++;
#endif

        }

#ifdef MPXDECODER_TUNE_PILOT_AGC
        peakPilotGain = std::max(pilotGain, peakPilotGain);
#endif

        // decay the peak hold
        pilotPeak *= peakDecreaseRate;

        // double pilot frequency. Note: amplitude approx 1/2 of full-scale (canonical doubler is 2x^2 - 1)
        FloatType doubledPilot = pilot * pilot - doublerDcOffset;

        // do the spectrum shift
        constexpr double scaling = 2.5 * 2 * 2; // 10.0
        FloatType side = scaling * doubledPilot * sideRaw;

		// separate L, R stereo channels
		FloatType left = 0.5 * (monoRaw + side);
		FloatType right = 0.5 * (monoRaw - side);

		if(!lowpassEnabled) {
			return {left, right};
		}

		// filter & return outputs
		filters.at(4).put(left);
		filters.at(5).put(right);
		return {filters.at(4).get(), filters.at(5).get()};
    }

    template <typename FloatType>
    static std::vector<FloatType> makeBandpass(int sampleRate, double ft1, double ft2)
    {
        // determine cutoff frequency and steepness
        double nyquist = sampleRate / 2.0;
        double steepness = 0.090909091 / (2000.0 / nyquist);

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

    // 15khz lowpass for audio output
    template<typename FloatType>
    static std::vector<FloatType> make15khzLowpass(int sampleRate)
    {
		// determine filter steepness
        double nyquist = sampleRate / 2.0;
		double steepness = 0.090909091 / (lpfW / nyquist);

        // determine filtersize
        int filterSize = static_cast<int>(
            std::min<int>(FILTERSIZE_BASE * steepness, FILTERSIZE_LIMIT)
            | 1 // ensure that filter length is always odd
        );

        std::vector<FloatType> filterTaps1(filterSize, 0);
        FloatType* pFilterTaps1 = &filterTaps1[0];
		ReSampler::makeLPF<FloatType>(pFilterTaps1, filterSize, lpfT, sampleRate);
        int sidelobeAtten = 160;
        ReSampler::applyKaiserWindow<FloatType>(pFilterTaps1, filterSize, ReSampler::calcKaiserBeta(sidelobeAtten));
        return filterTaps1;
    }

    // 19khz bandpass filter for the Pilot Tone
    template<typename FloatType>
    static std::vector<FloatType> make19KhzBandpass(int sampleRate)
    {
        return makeBandpass<FloatType>(sampleRate, 18900, 19100);
    }

    // 38khz bandpass filter for the Audio Subcarrier
    template<typename FloatType>
    static std::vector<FloatType> make38KhzBandpass(int sampleRate)
    {
        return makeBandpass<FloatType>(sampleRate, 38000, 53000); // we only want half of it
    }

    // 57khz bandpass filter for RDS / RBDS
    template<typename FloatType>
    static std::vector<FloatType> make57KhzBandpass(int sampleRate)
    {
        return makeBandpass<FloatType>(sampleRate, 54000, 60000);
    }

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

	static double getLpfT()
	{
		return lpfT;
	}

	static double getLpfW()
	{
		return lpfW;
	}

	bool getLowpassEnabled() const
	{
		return lowpassEnabled;
	}

	void setLowpassEnabled(bool value)
	{
		lowpassEnabled = value;
	}

private:
	std::vector<ReSampler::FIRFilter<double>> filters;

	static constexpr double lpfT = 15500.0;	// LPF transition freq (Hz)
	static constexpr double lpfW = 3500.0;	// LPF transition width (Hz)
    static constexpr double pilotStableLow = 0.98;
    static constexpr double pilotStableHigh = 0.99;
    static constexpr double doublerDcOffset = 0.5 * (pilotStableLow + pilotStableHigh) / 2;

    // if more gain than this is needed, then something is wrong with the Pilot Tone:
    static constexpr double pilotMaxGain = 25.0;

    std::vector<double> delayLine;
    int length;
    int currentIndex;
    int centerTap;

	bool lowpassEnabled; //{true}; // do final stereo 15khz LPF or not ?
    double pilotPeak{0.0};
    double pilotGain{1.0};
    double increaseRate;
    double decreaseRate;
    double peakDecreaseRate;

#ifdef MPXDECODER_TUNE_PILOT_AGC
    int64_t plusCount{0};
    int64_t minusCount{0};
    int64_t stableCount{0};
    double peakPilotGain{0.0};
#endif

};

#endif // MPXDECODE_H








