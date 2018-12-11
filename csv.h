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


//                makeTbl();
//                readHeaders();
//
//                if (err)
//                    return;
//
//                bufferSize = blockSize * numChannels;
//                inputBuffer = new uint8_t[bufferSize];
//                totalBytesRead = 0;
//                endOfBlock = bufferSize;
//                bufferIndex = endOfBlock; // empty (zero -> full)
//                currentBit = 0;
//                currentChannel = 0;
//                break;
//
//            case dff_write:
//                break;



    }

private:
    std::string path;
    CsvOpenMode mode;
    std::fstream file;
    bool err;
};


#endif //RESAMPLER_CSV_H
