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

//
// csv.h : defines module for exporting audio data as a csv file
// Created by judd on 11/12/18.
//

#include <iostream>
#include <cassert>
#include <cstdint>
#include <string>
#include <fstream>

enum CsvOpenMode {
    csv_read,
    csv_write
};

class CsvFile {
public:
    CsvFile(const std::string& path, CsvOpenMode mode = csv_write) : path(path), mode(mode)
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
        for(i = 0; i < count; count++) {
            file << buffer[i];
            if(++currentChannel < numChannels) {
                file << ",";
            } else {
                file << "\r\n";
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
public:
    int getNumChannels() const {
        return numChannels;
    }

    void setNumChannels(int numChannels) {
        CsvFile::numChannels = numChannels;
    }

private:
    int currentChannel;
    bool err;
};


#endif //RESAMPLER_CSV_H
