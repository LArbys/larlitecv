#!/bin/bash

# clean up previously set env
if [[ -z $FORCE_LARLITECV_BASEDIR ]]; then
    where="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
    export LARLITECV_BASEDIR=${where}
else
    export LARLITECV_BASEDIR=$FORCE_LARLITECV_BASEDIR
fi

if [[ -z $LARLITECV_BUILDDIR ]]; then
    export LARLITECV_BUILDDIR=$LARLITECV_BASEDIR/build
fi

export LARLITECV_COREDIR=$LARLITECV_BASEDIR/core
export LARLITECV_APPDIR=$LARLITECV_BASEDIR/app
export LARLITECV_LIBDIR=$LARLITECV_BUILDDIR/lib
export LARLITECV_INCDIR=$LARLITECV_BUILDDIR/include
export LARLITECV_BINDIR=$LARLITECV_BUILDDIR/bin

# Abort if ROOT not installed. Let's check rootcint for this.
if [ `command -v rootcling` ]; then
    export LARLITECV_ROOT6=1
else 
    if [[ -z `command -v rootcint` ]]; then
	echo
	echo Looks like you do not have ROOT installed.
	echo You cannot use LArLite w/o ROOT!
	echo Aborting.
	echo
	return 1;
    fi
fi

echo
printf "\033[93mLArCV\033[00m FYI shell env. may useful for external packages:\n"
printf "    \033[95mLARLITECV_INCDIR\033[00m   = $LARLITECV_INCDIR\n"
printf "    \033[95mLARLITECV_LIBDIR\033[00m   = $LARLITECV_LIBDIR\n"
printf "    \033[95mLARLITECV_BUILDDIR\033[00m = $LARLITECV_BUILDDIR\n"

# ADD TO PATH AND LD_LIBRARY_PATH LISTS (IF NOT FOUND)
[[ ":$PATH:" != *":${LARLITECV_BINDIR}:"* ]] && PATH="${LARLITECV_BINDIR}:${PATH}"
[[ ":$LD_LIBRARY_PATH:" != *":${LARLITECV_LIBDIR}:"* ]] && LD_LIBRARY_PATH="${LARLITECV_LIBDIR}:${LD_LIBRARY_PATH}"
[[ ":$DYLD_LIBRARY_PATH:" != *":${LARLITECV_LIBDIR}:"* ]] && DYLD_LIBRARY_PATH="${LARLITECV_LIBDIR}:${DYLD_LIBRARY_PATH}"

mkdir -p $LARLITECV_BUILDDIR;
mkdir -p $LARLITECV_LIBDIR;
mkdir -p $LARLITECV_BINDIR;

#export LD_LIBRARY_PATH=$LARLITECV_LIBDIR:$LD_LIBRARY_PATH
#export PYTHONPATH=$LARLITECV_BASEDIR/python:$PYTHONPATH
[[ ":$PYTHONPATH:" != *":${LARLITECV_BASEDIR}/python:"* ]] && PYTHONPATH="${LARLITECV_BASEDIR}/python:${PYTHONPATH}"

if [[ $LARLITE_BASEDIR ]]; then
    printf "\033[93mLArLite\033[00m\n"
    echo "    Found larlite set up @ \$LARLITE_BASEDIR=${LARLITE_BASEDIR}"
else
    printf "\033[93mLArLite\033[00m\n"
    echo "    Missing larlite. Required dependency."
    return 1;
fi

if [[ $LARCV_BASEDIR ]]; then
    printf "\033[93mLArLite\033[00m\n"
    echo "    Found larlite set up @ \$LARLITE_BASEDIR=${LARLITE_BASEDIR}"
else
    printf "\033[93mLArLite\033[00m\n"
    echo "    Missing larlite. Required dependency."
    return 1;
fi

export LARLITECV_CXX=clang++
if [ -z `command -v $LARLITECV_CXX` ]; then
    export LARLITECV_CXX=g++
    if [ -z `command -v $LARLITECV_CXX` ]; then
        echo
        echo Looks like you do not have neither clang or g++!
        echo You need one of those to compile LArCaffe... Abort config...
        echo
        return 1;
    fi
fi

echo
echo "Finish configuration. To build, type:"
echo "> cd \$LARLITECV_BUILDDIR"
echo "> make "
echo
