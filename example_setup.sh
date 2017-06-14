# LARCV AND LARLITE DIRS
MY_LARLITE=/Users/twongjirad/working/uboone/larlite
MY_LARCV=/Users/twongjirad/working/larbys/LArCV

# OPTIONAL: SPECIFY LOCATION OF OPENCV library and include directories
# This affects LArCV compilation
#export OPENCV_LIBDIR=/usr/local/lib
#export OPENCV_INCDIR=/usr/local/include

#SETUP LARLITE
cd $MY_LARLITE/config
source setup.sh
cd -

# SETUP LARCV
source $MY_LARCV/configure.sh

# setup LARLITECV
source configure.sh
