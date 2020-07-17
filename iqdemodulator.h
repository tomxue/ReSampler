/*
* Copyright (C) 2020 Judd Niemann - All Rights Reserved.
* You may use, distribute and modify this code under the
* terms of the GNU Lesser General Public License, version 2.1
*
* You should have received a copy of GNU Lesser General Public License v2.1
* with this file. If not, please refer to: https://github.com/jniemann66/ReSampler
*
*/

#ifndef IQDEMODULATOR_H
#define IQDEMODULATOR_H

#include <string>
#include <cstdint>
#include <memory>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <map>
#include <iostream>

#include <sndfile.h>
#include <sndfile.hh>

#include "biquad.h"
#include "mpxdecode.h"

namespace  ReSampler {

enum ModulationType
{
	ModulationTypeNone = 0,
	NFM,
	AM,
	LSB,
	USB,
	WFM,
	WFM_NO_LOWPASS,
	DSB,
	CW
};

static const std::map<std::string, ModulationType> modulationTypeMap
{
	{"NONE", ModulationType::ModulationTypeNone},
	{"NFM", ModulationType::NFM},
	{"AM", ModulationType::AM},
	{"LSB", ModulationType::LSB},
	{"USB", ModulationType::USB},
	{"WFM", ModulationType::WFM},
	{"DSB", ModulationType::DSB},
	{"CW", ModulationType::CW}
};

class IQFile
{
public:

	IQFile(const std::string& fileName) : sndfileHandle(new SndfileHandle(fileName))
	{
	}

	IQFile(const std::string& fileName, int infileMode, int infileFormat, int infileChannels, int infileRate) :
		sndfileHandle(new SndfileHandle(fileName, infileMode, infileFormat & 0xFFFF00FF, infileChannels, infileRate))
	{
		// Extract modulation type from 2nd-last byte of file format code.
		// (Note: libsndfile has this for the subformat mask:
		// SF_FORMAT_SUBMASK = 0x0000FFFF
		// So far, only the least-significant byte has been used. ie: 0x000000FF.
		// If they ever add more formats in the future which use the upper byte,
		// then this strategy may need reevaluation ...)

		modulationType = static_cast<ModulationType>((infileFormat & 0x0000FF00) >>  8);
		bool enableLowpass = true;
		if(modulationType == WFM_NO_LOWPASS) {
			enableLowpass = false;
			modulationType = WFM;
		}

		if(modulationType == WFM) {
			if(samplerate() != 0) {
				int sampleRate = sndfileHandle->samplerate();
				setDeEmphasisTc(2, sampleRate, 50);
				mpxDecoder = std::unique_ptr<MpxDecoder>(new MpxDecoder(sampleRate));
				mpxDecoder->setLowpassEnabled(enableLowpass);
			}
		}
	}

	bool error() {

		if(sndfileHandle == nullptr) {
			return true;
		}

		if(sndfileHandle->samplerate() == 0) {
			return true;
		}

		if(modulationType == WFM && sndfileHandle->samplerate() < 116000) {
			std::cout << "Sample Rate not high enough for WFM" << std::endl;
			return true;
		}

		if(sndfileHandle->channels() != 2) {
			std::cout << "2 channels expected for an I/Q input file !" << std::endl;
			return true;
		}

		return sndfileHandle->error();
	}

	int channels() {
		if(modulationType == WFM) {
			return 2; // FM stereo
		} else {
			return 1;
		}
	}

	int samplerate() {
		return sndfileHandle == nullptr ? 0 : sndfileHandle->samplerate();
	}

	int64_t frames() {
		return sndfileHandle == nullptr ? 0LL : sndfileHandle->frames();
	}

	int format() {
		return sndfileHandle == nullptr ? 0 : sndfileHandle->format();
	}

	template<typename FloatType>
	int64_t read(FloatType* inbuffer, int64_t count) {

		if(sndfileHandle == nullptr) {
			std::cout << "nullptr" << std::endl;
			return 0LL;
		}

		if(wavBuffer.size() < count) {
			wavBuffer.resize(count);
		}

		int64_t samplesRead = sndfileHandle->read(wavBuffer.data(), count);
		int64_t j = 0;

		switch(modulationType) {
		case AM:
			// Amplitude Modulation
			for(int64_t i = 0; i < samplesRead; i += 2) {
				inbuffer [j++] = demodulateAM(wavBuffer.at(i), wavBuffer.at(i + 1));
			}
			break;
		case LSB:
		case USB:
			// Single Side Band
			for(int64_t i = 0; i < samplesRead; i += 2) {
				// just copy I-component
				inbuffer[j++] = wavBuffer.at(i);
			}
			break;
		case WFM:
			// Wideband FM:
			for(int64_t i = 0; i < samplesRead; i += 2) {
				// demodulate, decode, deemphasize
				std::pair<FloatType, FloatType> decoded = mpxDecoder->decode(demodulateFM(wavBuffer.at(i), wavBuffer.at(i + 1)));
				inbuffer[j++] = deEmphasisFilters[0].filter(decoded.first);
				inbuffer[j++] = deEmphasisFilters[1].filter(decoded.second);
			}
			break;
		default:
			// Narrowband FM
			for(int64_t i = 0; i < samplesRead; i += 2) {
				inbuffer[j++] = demodulateFM(wavBuffer.at(i), wavBuffer.at(i + 1));
			}
		}

		return j;
	}

	sf_count_t seek(int64_t frames, int whence) {
		if(sndfileHandle == nullptr) {
			return 0LL;
		}

		return sndfileHandle->seek(frames, whence);
	}

	// getters
	ModulationType getModulationType() const
	{
		return modulationType;
	}

	// setters
	void setModulationType(const ModulationType &value)
	{
		modulationType = value;
	}

private:
	template<typename FloatType>
	FloatType demodulateFM(FloatType i, FloatType q)
	{
		static constexpr double maxGain = 60.0;
		static const double c = std::pow(10.0, -(maxGain / 20.0));

		// this is actually quite simple, thanks to some clever calculus tricks.
		// see https://www.embedded.com/dsp-tricks-frequency-demodulation-algorithms/

		i2 = i1;
		i1 = i0;
		i0 = i;
		q2 = q1;
		q1 = q0;
		q0 = q;

		double gain = 1.0 / (c + i * i + q * q);
		return gain * (((q0 - q2) * i1) - ((i0 - i2) * q1));
	}

    const std::vector<double> differentiator_coeffs{
        0.0035,
        -0.0140,
        0.0401,
        -0.1321,
        1.2639,
        -1.2639,
        0.1321,
        -0.0401,
        0.0140,
        -0.0035
    };


	template<typename FloatType>
	FloatType demodulateAM(FloatType i, FloatType q)
	{
		static constexpr FloatType scale = 0.7071; // << 1/sqrt(2)
		return scale * std::sqrt(i * i + q * q);
	}

	//	double tau = 1/(2*pi*f); // Hz to time constant
	//	double f = 2122.1; // 75 us
	//  double f = 3183.1; // 50 us

	void setDeEmphasisTc(int channels, int sampleRate, double tc = 50.0 /* microseconds */)
	{
		deEmphasisFilters.resize(channels);
		double p1 = -exp(-1.0 / (sampleRate * tc * 0.000001));
		double z1 = (1 + p1) / 5.0;
		for(auto& biquad : deEmphasisFilters)
		{
			biquad.setCoeffs(z1, z1, 0.0, p1, 0.0);
			biquad.reset();
		}
	}

private:
	// resources
	std::unique_ptr<SndfileHandle> sndfileHandle;
	std::unique_ptr<MpxDecoder> mpxDecoder;
	std::vector<double> wavBuffer;
	std::vector<Biquad<double>> deEmphasisFilters;

	// properties
	ModulationType modulationType{ModulationType::NFM};

	// registers used for demodulating FM
	double i0{0.0};
	double i1{0.0};
	double i2{0.0};
	double q0{0.0};
	double q1{0.0};
	double q2{0.0};

};

} // namespace  ReSampler

#endif // IQDEMODULATOR_H
