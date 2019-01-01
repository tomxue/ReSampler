/*
* Copyright (C) 2016 - 2019 Judd Niemann - All Rights Reserved.
* You may use, distribute and modify this code under the
* terms of the GNU Lesser General Public License, version 2.1
*
* You should have received a copy of GNU Lesser General Public License v2.1
* with this file. If not, please refer to: https://github.com/jniemann66/ReSampler
*/

#ifndef OSSPECIFIC_H
#define OSSPECIFIC_H 1

// macros which address differences between operating systems go here. 

#ifdef _WIN32

#if defined(_MSC_VER)
// MSVC specific
#define NOMINMAX // disable min() and max() macros (use std:: library instead)
#pragma warning(disable : 4996) // suppress pointless MS "deprecation" warnings
#pragma warning(disable : 4244) // suppress double-to-float warnings
#endif

#else // Non-Windows:
typedef uint64_t __int64;
#define stricmp strcasecmp
#endif // ends Non-Windows

#endif // OSSPECIFIC_H