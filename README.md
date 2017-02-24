# larlitecv
Analysis processor framework for working with LArLite and LArCV data

`larlitecv` attempts to provide some tools to synchronize input from larlite and larcv data files.  What exists now is a `DataCoordinator` class which allows one to get the right larlite and larcv data entry using run, subrun, and event values.  Hopefully sometime in the near future, support for larlite `ana_processor` and larcv `processmanager` will be added.

## Installation

larlitecv has the following required dependencies:

  * gcc or clang
  * ROOT (5 or 6, with 6 primarily supported at this point)
  * larlite (https://github.com/larlight/larlite)
  * larcv (https://github.com/LArbys/LArCV)

One must download and build the above before trying to compile larlitecv.  Once everything is compiled, one must setup the  environment for each each time a new shell is opened.  There is a script, `example_setup.sh`, in the base folder of larlitecv to assist in this process. 

### Building from scratch

Have compilers and stuff. Instructions on this is left to the internet. But on OS X install XCode and on Ubuntu install `gcc` and one should be pretty much set.

Install ROOT. This README will not cover how to do this. Please refer to ROOT documentation

Make a directory where all the magic will happen

    $ mkdir work

Install larlite

    $ git clone https://github.com/LArLight/larlite larlite
    $ cd larlite
    $ source config/setup.sh
    $ make
    [... if everything builds ...]
    $ cd UserDev/BasicTool
    $ make
    $ cd ../../..

Install larcv. (Note shell must have larlite environment variables set. To set them up, run `source config/setup.sh` from larlite base directory.  once you do this, the environment variable `LARLITE_BASEDIR` is defined.)

    $ git clone https://github.com/LArCV larcv
    $ cd larcv
    $ source configure.sh
    $ make
    cd ..
  
Install larlitecv.

    $ git clone https://github.com/larlitecv larlitecv
    $ cd larlitecv
    $ source configure.sh
    $ make
    cd ..
    
You should have the following directories in `work`

    $ ls
    larlite larcv larlitecv
    
### Resetup the shell environment

    [goto larlite folder]
    $ source config/setup.sh
    [goto larcv folder]
    $ source configure.sh
    [goto larlitecv folder]
    $ source configure.sh
    
Alternatively, make a copy of `example_setup.sh` and change the `MY_LARLITE` and `MY_LARCV` to point to the right folders on your computer. Then run the script via

    source my_example_setup.sh
    
### Quick verification everything is in order

    $ root-config --libs
    -L[path to root directory] -lCore -lRIO ...
    $ larlite-config --includes
    -I[path to larlite diretory]/core
    $ larcv-config --includes
    -I[path to larcv directory]/build/include -I[path to opencv include dir] -I[path to python include dir] ...

### Optional Dependencies

LArCV has an optional dependency on OpenCV, which gives access to some nice image processing functions.  To use this, when compiling LArCV, one needs to define the following environment variables in the shell:

    OPENCVLIBDIR
    OPENCVINCDIR

These variables should point to the directory where the shared object libraries and headers for OpenCV live, respectively.

For example, on Ubuntu 16.10, this is 

    OPENCVLIBDIR=/usr/local/lib
    OPENCVINCDIR=/usr/local/include
    
For those on OS X, one easy way is to install OpenCV using home brew or macports (under no circumstances should you be using both package managers! Pick one!).  For home brew the above, you can find this info using `brew info opencv3`
  
  
## Getting Started

There are two examples available in the app folder: `Example` and `SelectionExample`.

You can also make a minimial application using the command `gen_larlitecv_app.py`. For example,

    $ gen_larlitecv_app.py test


This will make a new folder, `app/test`, which provides

  * an empty class header, `test.h`
  * an empty class source, `test.cxx
  * an executable source file, `run_test.cxx`
  * a GNUmakfile
  * a LinkDef.h file (which helps ROOT build the symbols it needs to have your class be useable in the ROOT interpreter or callable in pyROOT)

You should be able to run `make` in this folder
