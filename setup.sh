# LARCV AND LARLITE DIRS
MY_LARLITE=/home/taritree/working/larbys/dev/larlite
MY_LARCV=/home/taritree/working/larbys/dev/larcv

export OPENCV_INCDIR=/usr/local/include
export OPENCV_LIBDIR=/usr/local/lib

#SETUP LARLITE
source $MY_LARLITE/config/setup.sh

# SETUP LARCV
source $MY_LARCV/configure.sh

# setup LARLITECV
source configure.sh

# setup pylard (not required)
#source /home/taritree/working/larbys/dev/pylard/env/bin/activate
