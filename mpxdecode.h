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

template<typename FloatType>
struct FMBasebandFilterResult
{
	FloatType r0;
	FloatType r1;
	FloatType r2;
	FloatType r3;
};

template<typename FloatType>
class FMBasebandFilter
{
public:
	FMBasebandFilter(int sampleRate)
	{
		// FIR filters must all be linear-phase and of odd, equal length
		f1 = make19KhzBandpass(sampleRate);
		f2 = make38KhzBandpass(sampleRate);
		f3 = make57KhzBandpass(sampleRate);
		length = f1.size();
        h.resize(length, 0.0);
		centerTap = (length - 1) / 2;
		currentIndex = length - 1;
	}

	FMBasebandFilterResult<FloatType> filter(FloatType input)
	{
		FloatType s0 = h.at(centerTap);
		FloatType s1{0.0};
		FloatType s2{0.0};
		FloatType s3{0.0};

		h[currentIndex] = input; // place input into history
        int d = currentIndex + centerTap;
        s0 = h[d >= length ? d - length : d]; // delay only
		int p = currentIndex;
		for(int j = 0 ; j < length; j++) {
			FloatType v = h.at(p++);
            if(p == length) {
                p = 0;
            }
            s1 += f1.at(j) * v;
			s2 += f2.at(j) * v;

		//	s3 += f3.at(j) * v;
		}

		if(currentIndex == 0) {
			currentIndex = length - 1;
		} else {
			currentIndex--;
		}

		return {s0, s1, s2, s3};
	}

private:
	int length;
	int currentIndex;
	int centerTap;
	std::vector<FloatType> f1;
	std::vector<FloatType> f2;
	std::vector<FloatType> f3;
	std::vector<FloatType> h;

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
	static std::vector<FloatType> make19KhzBandpass(int sampleRate)
	{
		return makeBandpass(sampleRate, 18900, 19100);
	}

	// 38khz bandpass filter for the Audio Subcarrier
	static std::vector<FloatType> make38KhzBandpass(int sampleRate)
	{
		return makeBandpass(sampleRate, 38000, 53000); // we only want half of it
	}

	// 57khz bandpass filter for RDS / RBDS
	static std::vector<FloatType> make57KhzBandpass(int sampleRate)
	{
		return makeBandpass(sampleRate, 54000, 60000);
	}
public:
	// function for testing / tweaking the performance of the bandpass filters
	static void saveFilters1(const std::string& filename)
	{
		std::vector<FloatType> filt1 = make19KhzBandpass(192000);
		std::vector<FloatType> filt2 = make38KhzBandpass(192000);
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
};

template<typename FloatType>
class FMAudioFilter
{
public:
	FMAudioFilter(int sampleRate)
	{
		make15khzLowpass(sampleRate);
	}

	std::pair<FloatType, FloatType> filter(FloatType left, FloatType right)
	{
        filters[0].put(left);
        filters[1].put(right);
        return {filters[0].get(), filters[1].get()};
	}

private:
    std::vector<ReSampler::FIRFilter<FloatType>> filters;
	void make15khzLowpass(int sampleRate)
	{
		// determine cutoff frequency and steepness
		double nyquist = sampleRate / 2.0;
		double steepness = 0.090909091 / (1000.0 / nyquist);

		// determine filtersize
		int filterSize = static_cast<int>(
			std::min<int>(FILTERSIZE_BASE * steepness, FILTERSIZE_LIMIT)
			| 1 // ensure that filter length is always odd
		);

        std::vector<FloatType> lpf(filterSize, 0.0);
		ReSampler::makeLPF<FloatType>(lpf.data(), filterSize, 15500, sampleRate);
		int sidelobeAtten = 160;
		ReSampler::applyKaiserWindow<FloatType>(lpf.data(), filterSize, ReSampler::calcKaiserBeta(sidelobeAtten));
        filters.emplace_back(lpf.data(), lpf.size()); // left lowpass
        filters.emplace_back(lpf.data(), lpf.size()); // right lowpass
	}
};

class MpxDecoder
{
public:
	MpxDecoder(int sampleRate) : basebandFilter(sampleRate), audioFilter(sampleRate)
    { 
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
		FMBasebandFilterResult<FloatType> filtered = basebandFilter.filter(input);
		FloatType monoRaw = filtered.r0;
		FloatType pilotRaw = filtered.r1;
		FloatType sideRaw = filtered.r2;
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

		// separate L, R and return filtered result
		return audioFilter.filter(0.5 * (monoRaw + side), 0.5 * (monoRaw - side));
    }

private:
	FMBasebandFilter<double> basebandFilter;
    std::vector<ReSampler::FIRFilter<double>> audioFilters;
	FMAudioFilter<double> audioFilter;

    static constexpr double pilotStableLow = 0.98;
    static constexpr double pilotStableHigh = 0.99;
    static constexpr double doublerDcOffset = 0.5 * (pilotStableLow + pilotStableHigh) / 2;

    // if more gain than this is needed, then something is wrong with the Pilot Tone:
    static constexpr double pilotMaxGain = 25.0;

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
