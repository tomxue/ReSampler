#ifndef OSSPECIFIC_H
#define OSSPECIFIC_H 1

// macros which address differences between operating systems go here

#ifdef _WIN32
#define ZERO_64 0i64

#else // Non-Windows:
#define ZERO_64 0LL
#define stricmp strcasecmp
#endif // ends Non-Windows

#endif // !OSSPECIFIC_H

