# LARCV AND LARLITE DIRS
MY_LARLITE=/Users/twongjirad/working/uboone/larlite
MY_LARCV=/Users/twongjirad/working/larbys/LArCV

# Set these to use OpenCV. Needed by StopMu Tagger
#export OPENCV_INCDIR=/usr/local/include
#export OPENCV_LIBDIR=/usr/local/lib

#SETUP LARLITE
cd $MY_LARLITE/config
source setup.sh
cd -

# SETUP LARCV
source $MY_LARCV/configure.sh

# setup LARLITECV
source configure.sh
