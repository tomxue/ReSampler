# Prerequisites

XCode


# building ReSampler (tested on macOS 10.12 - Sierra)

Building on linux is fairly straightforward and consists of the following 3 steps:

- install fftw library
- install libsndfile library
- Build ReSampler


## install homebrew

~~~
ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)" < /dev/null 2> /dev/null
~~~

## install fftw

~~~
brew install fftw
~~~

## install libsndfile

~~~
brew install libsndfile
~~~

## building ReSampler

clone this repository to a local directory, 

~~~
git clone https://github.com/jniemann66/ReSampler.git
~~~

and use one of the following command lines to compile:

#### using clang:
~~~
clang++ -pthread -std=c++11 ReSampler.cpp -lfftw3 -lsndfile -o ReSampler-clang -O3 -L/usr/local/lib -I/usr/local/include
~~~

#### using gcc:

standard 64-bit build:
~~~
g++ -pthread -std=c++11 ReSampler.cpp -lfftw3 -lsndfile -o ReSampler -O3 -L/usr/local/lib -I/usr/local/include
~~~

AVX Build **(NOT TESTED YET !)**:
~~~
g++ -pthread -std=c++11 ReSampler.cpp -lfftw3 -lsndfile -o ReSampler -O3 -DUSE_AVX -mavx -L/usr/local/lib -I/usr/local/include
~~~