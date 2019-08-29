# building ReSampler for Raspbian on Raspberry Pi

Building on Raspbian consists of the following 4 steps:

- Install development tools
- Build and install fftw library
- Build and install libsndfile library
- Build and install ReSampler

#### build environment
~~~
sudo apt-get install build-essential
~~~

#### fftw build/install

~~~
#sudo apt-get install libfftw3-dev libfftw3-doc
~~~

#### libsndfile build/install

~~~
sudo apt-get install libsndfile-dev
~~~

## building ReSampler

#### using gcc

clone this repository to a local directory, and invoke gcc

~~~
g++ -pthread -std=c++11 main.cpp ReSampler.cpp conversioninfo.cpp -lfftw3 -lsndfile -o ReSampler -O3
~~~

#### using CMake

~~~
mkdir some-directory
cd some-directory
cmake -DCMAKE_BUILD_TYPE=Release path-to-ReSampler-source
make
~~~
