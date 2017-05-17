# building ReSampler (tested on Ubuntu 16.04)

## build environment
~~~
sudo apt-get install build-essential
~~~

## fftw build/install

#sudo apt-get install libfftw3-dev libfftw3-doc

[fftw source](http://www.fftw.org/fftw-3.3.6-pl2.tar.gz)

*extract it to projects folder, and cd to it*
~~~
./configure
make
sudo make install
~~~

## libsndfile build/install

[libsndfile 1.0.28](http://www.mega-nerd.com/libsndfile/files/libsndfile-1.0.28.tar.gz)

*extract it to projects folder, and cd to it*

~~~
sudo apt install autoconf autogen automake build-essential libasound2-dev \
    libflac-dev libogg-dev libtool libvorbis-dev pkg-config python

./configure --prefix=/usr    \
            --disable-static \
            --docdir=/usr/share/doc/libsndfile-1.0.28 &&
make
make check
sudo make install
~~~

## building ReSampler
no SSE, no optimization:
~~~
g++ -pthread -std=c++11 ReSampler.cpp -lfftw3 -lsndfile -o ReSampler
~~~

SSE2, O3 optimization:
~~~
g++ -pthread -std=c++11 ReSampler.cpp -lfftw3 -lsndfile -o ReSampler -D USE_SSE2 -D SSE_CUSTOM_HSUM -O3
~~~


# misc:

## show where gcc is looking for header files:
~~~
`gcc -print-prog-name=cc1plus` -v
`gcc -print-prog-name=cc1` -v
~~~

## show where gcc is looking for libraries:
~~~
gcc -print-search-dirs
~~~

## misc
~~~
export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig
./autogen.sh
./configure --enable-werror
make
make check
~~~