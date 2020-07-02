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

namespace  ReSampler {

enum ModulationType
{
	ModulationTypeNone = 0,
	NFM,
	AM,
	LSB,
	USB,
	WFM,
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
		if(modulationType == WFM) {
			setDeEmphasisTc(2, sndfileHandle->samplerate(), 50);
		}
	}

	bool error() {

		if(sndfileHandle == nullptr) {
			return true;
		}

		if(sndfileHandle->channels() != 2) {
			std::cout << "2 channels expected for an I/Q input file !" << std::endl;
			return true;
		}

		return sndfileHandle->error();
	}

	int channels() {
		return 1;  // I & Q inputs always get demodulated into a single channel
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

		if(wavBuffer.size() < count * 2) {
			wavBuffer.resize(count * 2);
		}

		int64_t samplesRead = sndfileHandle->read(wavBuffer.data(), count * 2);

		int64_t j = 0;
		for(int64_t i = 0; i < samplesRead; i += 2) {
			switch(modulationType) {
			case AM:
				inbuffer[j++] = demodulateAM(wavBuffer.at(i), wavBuffer.at(i + 1));
				break;
            case LSB:
            case USB:
                // SSB : just copy I-component
                inbuffer[j++] = wavBuffer.at(i);
                break;
			case WFM:
				// wideband FM
				// todo: multiplex decode
				inbuffer[j++] = deEmphasisFilters[0].filter(
						demodulateFM(wavBuffer.at(i), wavBuffer.at(i + 1))
				);
				break;
			default:
                // Narrowband FM
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

	template<typename FloatType>
	FloatType demodulateAM(FloatType i, FloatType q)
	{
		static constexpr FloatType scale = 0.7071; // << 1/sqrt(2)
		return scale * std::sqrt(i * i + q * q);
	}

	//		double f = 2122.1; // 75 us
	//		double T1 = 1/(2*pi*fc); // Hz to time constant

	void setDeEmphasisTc(int channels, int sampleRate, double tc = 50.0 /* microseconds */)
	{
		deEmphasisFilters.resize(channels);
		double p1 = -exp(-1.0 / (sampleRate * tc * 0.000001));
		double z1 = (1 + p1) / 5.0;
		for(auto& biquad : deEmphasisFilters)
		{
			biquad.setCoeffs(z1, z1, 0, p1, 0);
			std::cout << sampleRate << "," << std::setprecision(9) << z1 << "," << p1 << std::endl;
			// 0.06501945611827269, 0.06501945611827269, 0, 0.8699610877634546, 0// b0, b1, b2, a1, a2
			biquad.reset();
		}
	}

private:
	// resources
	std::unique_ptr<SndfileHandle> sndfileHandle;
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
