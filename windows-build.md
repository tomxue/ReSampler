## Building ReSampler on Windows

#### using Visual Studio

- Just open the project file (ReSampler.vcxproj) in Visual Studio
- Choose the configuration you want (eg Release x64)
- BUILD !!

#### using minGW-W64 under Windows (to build 64-bit .exe):
- use MinGW-W64 with Posix Threads and SEH
- {...}\MinGW-W64\mingw64\bin must be added to Path
- unlike *nix systems, there is no standard place to put libraries and include files
- locations of include files can be specified by -Idir
- locations of library files can be specified by -Ldir
- if libraries are .lib (instead of .a), use -llib&lt;name&gt; instead of -l&lt;name&gt; 
- launch from git bash
- depends on these 64-bit dlls: **libfftw3-3.dll  libgcc_s_seh-1.dll  libsndfile-1.dll  libstdc++-6.dll  libwinpthread-1.dll**

standard 64-bit build:
~~~
g++ -pthread -std=c++11 ReSampler.cpp -Ilibsndfile/include -Ifftw64 -Lfftw64 -llibfftw3-3 -Llibsndfile/lib -llibsndfile-1 -o x64/minGW-W64/ReSampler.exe -O3
~~~

AVX build:
~~~
g++ -pthread -std=c++11 ReSampler.cpp -Ilibsndfile/include -Ifftw64 -Lfftw64 -llibfftw3-3 -Llibsndfile/lib -llibsndfile-1 -o x64/minGW-W64-AVX/ReSampler.exe -O3 -DUSE_AVX -mavx
~~~