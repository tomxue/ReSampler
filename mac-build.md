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



#### using gcc:

standard 64-bit build:
~~~
g++ -pthread -std=c++11 ReSampler.cpp -lfftw3 -lsndfile -o ReSampler -O3
~~~

AVX Build **(NOT TESTED YET !)**:
~~~
g++ -pthread -std=c++11 ReSampler.cpp -lfftw3 -lsndfile -o ReSampler -O3 -DUSE_AVX -mavx
~~~

#### using clang:
~~~
clang++ -pthread -std=c++11 ReSampler.cpp -lfftw3 -lsndfile -o ReSampler-clang -O3
~~~

# misc tasks:

## setting up C++ environment in vscode

ensure that the vscode [c++ tools](https://marketplace.visualstudio.com/items?itemName=ms-vscode.cpptools) are installed.

the g++ commands can be put into the **tasks.json** file. For example.:

~~~
{
    "version": "0.1.0",
    "command": "g++",
    "isShellCommand": true,
    "showOutput": "always",
    "args": [
        "-pthread",
        "-std=c++11",
        "ReSampler.cpp",
        "-l",
        "fftw3",
        "-l",
        "sndfile",
        "-o",
        "ReSampler",
        "-O3",
        "-v"
    ]
}
~~~

[documentation](https://code.visualstudio.com/docs/languages/cpp)

## show where gcc is looking for header files:
~~~
`gcc -print-prog-name=cc1plus` -v
`gcc -print-prog-name=cc1` -v
~~~

## show where gcc is looking for libraries:
~~~
gcc -print-search-dirs
~~~

## display info about the binary you just built:
file ReSampler

*sample output:*

~~~
ReSampler: ELF 64-bit LSB executable, x86-64, version 1 (SYSV), dynamically linked, interpreter /lib64/ld-linux-x86-64.so.2, for GNU/Linux 2.6.32, BuildID[sha1]=2873279a9b0040a268f7c485de9027660ab3617c, not stripped
~~~