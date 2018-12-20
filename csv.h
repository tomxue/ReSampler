/*
* Copyright (C) 2016 - 2018 Judd Niemann - All Rights Reserved.
* You may use, distribute and modify this code under the
* terms of the GNU Lesser General Public License, version 2.1
*
* You should have received a copy of GNU Lesser General Public License v2.1
* with this file. If not, please refer to: https://github.com/jniemann66/ReSampler
*/

#ifndef RESAMPLER_CSV_H
#define RESAMPLER_CSV_H

// csv.h : defines module for exporting audio data as a csv file

#include <iostream>
#include <cassert>
#include <cstdint>
#include <string>
#include <fstream>
#include <cmath>

enum CsvOpenMode {
    csv_read,
    csv_write
};

enum CsvSignedness {
	Signed,
	Unsigned
};

enum CsvNumericFormat {
	Integer,
	FloatingPoint,
	Fixed,
	Scientific
};

enum CsvNumericBase {
	Binary = 2,
	Octal = 8,
	Decimal = 10,
	Hexadecimal = 16
};

enum IntegerWriteScalingStyle {
    Pow2Clip,
    Pow2Minus1
};

class CsvFile {
public:
    CsvFile(const std::string& path, CsvOpenMode mode = csv_write) : path(path), mode(mode), signedness(Signed), numericBase(Decimal), numBits(16), numSignificantDigits(10), integerWriteScalingStyle(Pow2Minus1)
    {
        file.exceptions(std::ifstream::failbit | std::ifstream::badbit);
        currentChannel = 0;

        switch (mode) {
            case csv_read:
                try {
                    file.open(path, std::ios::in | std::ios::binary);
                    err = false;
                }
                catch (std::ios_base::failure& e) {
                    e.what();
                    err = true;
                    return;
                }
                break;

            case csv_write:
                try {
                    file.open(path, std::ios::out | std::ios::binary);
                    err = false;
                }
                catch (std::ios_base::failure& e) {
                    e.what();
                    err = true;
                    return;
                }
                break;
        }
    }

    ~CsvFile() {
        if(file.is_open()) {
            file.close();
        }
    }

    template <typename T>
    int64_t write(const T* buffer, int64_t count) {
        int64_t i;
        for(i = 0; i < count; i++) {
            switch(numericFormat) {
                case CsvNumericFormat::FloatingPoint:
                    file << buffer[i];
                    break;

                default:
                    file << scaleToInt(buffer[i]);
                    break;
            }

            if(++currentChannel < numChannels) {
                file << ",";
            } else {
                file  << "\r\n";
                currentChannel = 0;
            }
        }
        return i;
    }

private:
    std::string path;
    CsvOpenMode mode;
    std::fstream file;
	int numChannels;
	CsvNumericFormat numericFormat;
	double scaleFactor;
    CsvSignedness signedness;
    CsvNumericBase numericBase;
    int numBits;
    int numSignificantDigits;
    IntegerWriteScalingStyle integerWriteScalingStyle;
    int currentChannel;
    bool err;

	template <typename IntType, typename FloatType>
	IntType scaleToInt(FloatType x) {
	    return static_cast<IntType>(std::round(scaleFactor * x));
	}

	void setStreamFormat() {
	    if(file.is_open()) {
	        if(numericFormat == FloatingPoint) {
	            file.unsetf(std::ios::fixed);
	            file << std::setprecision(numSignificantDigits);
	        } else {
	            file.setf(std::ios::fixed);
	            file << std::setprecision(0);
	        }
	    }
	}

public:
    CsvNumericFormat getNumericFormat() const {
        return numericFormat;
    }

    void setNumericFormat(CsvNumericFormat numericFormat) {
        CsvFile::numericFormat = numericFormat;
        setStreamFormat();
    }

    CsvSignedness getSignedness() const {
        return signedness;
    }

    void setSignedness(CsvSignedness signedness) {
        CsvFile::signedness = signedness;
        setStreamFormat();
    }

    CsvNumericBase getNumericBase() const {
        return numericBase;
    }

    void setNumericBase(CsvNumericBase numericBase) {
        CsvFile::numericBase = numericBase;
        setStreamFormat();
    }

    int getNumBits() const {
        return numBits;
    }

    void setNumBits(int numBits) {

        scaleFactor = 1 << (numBits - 1);
        if(integerWriteScalingStyle == Pow2Minus1) {
            scaleFactor -= 1;
        }

        CsvFile::numBits = numBits;
        setStreamFormat();
    }

    int getNumSignificantDigits() const {
        return numSignificantDigits;
    }

    void setNumSignificantDigits(int numSignificantDigits) {
        CsvFile::numSignificantDigits = numSignificantDigits;
        setStreamFormat();
    }

    IntegerWriteScalingStyle getIntegerWriteScalingStyle() const {
        return integerWriteScalingStyle;
    }

    void setIntegerWriteScalingStyle(IntegerWriteScalingStyle integerWriteScalingStyle) {
        CsvFile::integerWriteScalingStyle = integerWriteScalingStyle;
    }

    int getNumChannels() const {
        return numChannels;
    }

    void setNumChannels(int numChannels) {
        CsvFile::numChannels = numChannels;
    }

};


#endif //RESAMPLER_CSV_H
