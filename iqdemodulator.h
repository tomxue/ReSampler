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

//#define COLLECT_IQ_STATS

#define ERROR_IQFILE_WFM_SAMPLERATE_TOO_LOW (0xff01)
#define ERROR_IQFILE_TWO_CHANNELS_EXPECTED (0xff02)

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

enum DeEmphasisType
{
	NoDeEmphasis = 0,
	DeEmphasis50 = 0x10,
	DeEmphasis75 = 0x20
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

static const std::map<std::string, DeEmphasisType> deEmphasisTypeMap
{
	{"NONE", NoDeEmphasis},
	{"50", DeEmphasis50},
	{"75", DeEmphasis75}
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

		int format = (infileFormat & 0x0000FF00) >>  8;
		modulationType = static_cast<ModulationType>(format & 0x0f);
		deEmphasisType = static_cast<DeEmphasisType>(format & 0x30);

        differentiatorLength = differentiatorCoeffs.size();
        differentiatorDelay = differentiatorLength / 2;
        historyI.resize(differentiatorLength, 0.0);
        historyQ.resize(differentiatorLength, 0.0);
        differentiatorIndex = differentiatorLength - 1;
        //std::cout << differentiatorLength << std::endl;
        //std::cout << differentiatorLength / 2 << std::endl;

		bool enableLowpass = true;
		if(modulationType == WFM_NO_LOWPASS) {
			enableLowpass = false;
			modulationType = WFM;
		}

		if(modulationType == WFM) {
			if(samplerate() != 0) {
				int sampleRate = sndfileHandle->samplerate();

				switch (deEmphasisType) {
				case NoDeEmphasis:
					break;
				case DeEmphasis50:
					setDeEmphasisTc(2, sampleRate, 50);
					break;
				case DeEmphasis75:
					setDeEmphasisTc(2, sampleRate, 75);
					break;
				}

				mpxDecoder = std::unique_ptr<MpxDecoder>(new MpxDecoder(sampleRate));
				mpxDecoder->setLowpassEnabled(enableLowpass);
			}

		}
	}

#ifdef COLLECT_IQ_STATS
	~IQFile()
	{
		double iAvg = sumI / framesRead;
		double qAvg = sumQ / framesRead;
		double iRMS = std::sqrt(sumISquared / framesRead);
		double qRMS = std::sqrt(sumQSquared / framesRead);
		std::cout << "I: peak=" << peakI << ", avg=" << iAvg << ", RMS=" << iRMS << std::endl;
		std::cout << "Q: peak=" << peakQ << ", avg=" << qAvg << ", RMS=" << qRMS << std::endl;
	}
#endif

	int error() {

		if(sndfileHandle == nullptr) {
			return SF_ERR_UNRECOGNISED_FORMAT;
		}

		if(sndfileHandle->samplerate() == 0) {
			return SF_ERR_UNRECOGNISED_FORMAT;
		}

		if(modulationType == WFM && sndfileHandle->samplerate() < 116000) {
			return ERROR_IQFILE_WFM_SAMPLERATE_TOO_LOW;
		}

		if(sndfileHandle->channels() != 2) {
			return ERROR_IQFILE_TWO_CHANNELS_EXPECTED;
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
			if(deEmphasisType == NoDeEmphasis) {
				for(int64_t i = 0; i < samplesRead; i += 2) {
					// demodulate, decode
					std::pair<FloatType, FloatType> decoded = mpxDecoder->decode(demodulateFM(wavBuffer.at(i), wavBuffer.at(i + 1)));
					inbuffer[j++] = decoded.first;
					inbuffer[j++] = decoded.second;
				}
			} else {
				for(int64_t i = 0; i < samplesRead; i += 2) {
					// demodulate, decode, deemphasize
					double iVal = wavBuffer.at(i);
					double qVal = wavBuffer.at(i + 1);

#ifdef COLLECT_IQ_STATS
					peakI = std::max(peakI, std::abs(iVal));
					peakQ = std::max(peakQ, std::abs(qVal));
					sumI += iVal;
					sumQ += qVal;
					sumISquared += (iVal * iVal);
					sumQSquared += (qVal * qVal);
					framesRead++;
#endif

                    std::pair<FloatType, FloatType> decoded = mpxDecoder->decode(demodulateFM(iVal, qVal));
					inbuffer[j++] = deEmphasisFilters[0].filter(decoded.first);
					inbuffer[j++] = deEmphasisFilters[1].filter(decoded.second);
				}
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
		static constexpr double threshold = -45.0;
		static const double c = std::pow(10.0, threshold / 20.0);

		// this is actually quite simple, thanks to some clever calculus tricks.
		// see https://www.embedded.com/dsp-tricks-frequency-demodulation-algorithms/

		i2 = i1;
		i1 = i0;
		i0 = i;
		q2 = q1;
		q1 = q0;
		q0 = q;

		double gain = 2.0 / (c + i1 * i1 + q1 * q1);
		return gain * (((q0 - q2) * i1) - ((i0 - i2) * q1));
	}

	template<typename FloatType>
	FloatType demodulateFM3(FloatType i, FloatType q)
	{
		static constexpr double threshold = -40.0; // dB
		static const double c = std::pow(10.0, threshold / 20.0);
		static const double gainTrim = 2.75494098472591;

		// determine magnitude and gain
		double iSquared = i * i;
		double qSquared = q * q;
		double a = std::sqrt(iSquared + qSquared);
		double g = 2.0 / std::max(a, c);

		// place input into history
		z0.real(i * g);
		z0.imag(q * g);

		auto dz = z0 * std::conj(z1);
		z1 = z0;
		return gainTrim * std::arg(dz);
	}

    template<typename FloatType>
    FloatType demodulateFM4(FloatType i, FloatType q)
    {
        static constexpr double threshold = -90.0;
        static const double c = std::pow(10.0, threshold / 20.0);
        FloatType dI{0.0}; // differentiated I
        FloatType dQ{0.0}; // differentiated Q

        // place input into history
        historyI[differentiatorIndex] = i;
        historyQ[differentiatorIndex] = q;

        // get position of delay tap
        FloatType delayedI;
        FloatType delayedQ;
        int delayIndex = differentiatorIndex + differentiatorDelay;
        if(delayIndex >= differentiatorLength) {
            delayIndex -= differentiatorLength;
        }

        // get delayed values from history
        delayedI = historyI.at(delayIndex);
        delayedQ = historyQ.at(delayIndex);

        // perform the convolution
        int p = differentiatorIndex;
        for(int j = 0 ; j < differentiatorLength; j++) {
            FloatType vI = historyI.at(p);
            FloatType vQ = historyQ.at(p);
            if(++p == differentiatorLength) {
                p = 0; // wrap
            }
            dI += differentiatorCoeffs[j] * vI;
            dQ += differentiatorCoeffs[j] * vQ;
        }

        // update the current index
        if(differentiatorIndex == 0) {
            differentiatorIndex = differentiatorLength - 1; // wrap
        } else {
            differentiatorIndex--;
        }

        double gain = 2.0 / (c + delayedI * delayedI + delayedQ * delayedQ);
        return gain * (dQ * delayedI - dI  * delayedQ);
    }

	template<typename FloatType>
	FloatType demodulateAM(FloatType i, FloatType q)
	{
		static constexpr FloatType scale = 0.7071; // << 1/sqrt(2)
		return scale * std::sqrt(i * i + q * q);
	}

	//	double tau = 1/(2*pi*f); // Hz to time constant
	//	double f = 2122.1; // 75 us
	//  double f = 3183.1; // 50 us

	void setDeEmphasisHz(int channels, int sampleRate, double freqHz)
	{
		setDeEmphasisTc(channels, sampleRate, (1.0 / (2.0 * M_PI * freqHz)));
	}

	void setDeEmphasisTc(int channels, int sampleRate, double tc = 50.0 /* microseconds */)
	{
		deEmphasisFilters.resize(channels);
		double p1 = -exp(-1.0 / (sampleRate * tc * 0.000001));
		double z1 = (1 + p1) / 5.0;
		for(auto& biquad : deEmphasisFilters) {
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
	DeEmphasisType deEmphasisType{DeEmphasis50};

	// registers used for demodulating FM (atan2-less)
	double i0{0.0};
	double i1{0.0};
	double i2{0.0};
	double q0{0.0};
	double q1{0.0};
	double q2{0.0};

	// registers for demodulating FM (atan2 version)
	std::complex<double> z0{0.0};
	std::complex<double> z1{0.0};

    //
    const std::vector<double> differentiatorCoeffs
    {
        0.0035,
        0.0,
        -0.0140,
        0.0,
        0.0401,
        0.0,
        -0.1321,
        0.0,
        1.2639,
        0.0,
        -1.2639,
        0.0,
        0.1321,
        0.0,
        -0.0401,
        0.0,
        0.0140,
        0.0,
        -0.0035
    };

    int differentiatorLength;
    int differentiatorDelay;
    int differentiatorIndex;
    std::vector<double> historyI;
    std::vector<double> historyQ;


#ifdef COLLECT_IQ_STATS
	int64_t framesRead{0ll};
	double peakI{0.0};
	double peakQ{0.0};
	double sumI{0.0};
	double sumQ{0.0};
	double sumISquared{0.0};
	double sumQSquared{0.0};
#endif

};

} // namespace  ReSampler

#endif // IQDEMODULATOR_H
