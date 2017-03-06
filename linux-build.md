## fftw build/install
~~~
./configure
make
make install
~~~

## libsndfile build/install
~~~
export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig
./autogen.sh
./configure --enable-werror
make
make check
~~~