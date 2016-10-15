# WORKS ON NUDOT

# SETUP ROOT
source /home/taritree/setup_root6.sh

# SETUP LARLITE
MY_LARLITE=/home/taritree/working/larlite
cd $MY_LARLITE/config
source setup.sh
cd -

# SETUP LARCV
MY_LARCV=/home/taritree/working/larbys/LArCV
source $MY_LARCV/configure.sh

# SETUP LARLITECV
cd ../..
source configure.sh
cd -

# SETUP PYLARD
source /home/taritree/working/pylard/env/bin/activate
