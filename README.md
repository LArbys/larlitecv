# larlitecv

Analysis processor framework for working with LArLite and LArCV data

`larlitecv` attempts to provide some tools to synchronize input from larlite and larcv data files.  What exists now is a `DataCoordinator` class which allows one to get the right larlite and larcv data entry using run, subrun, and event values.  Hopefully sometime in the near future, support for larlite `ana_processor` and larcv `processmanager` will be added.

### Muon Tagger Branch

This branch contains the muon-tagger code.

* app/ThruMu: through-going muon tagger using pairs of boundary end points
* app/StopMu: tags stopping muon tracks using single boundary end points. relies on output of ThruMu code
* app/ContainedROI: finds contained ROI tracks. uses output of StopMu code

Slides summarizing the algorithms can be found on the [uboone docdb](
http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=7390)

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

Install larlite. Note that in addition to the core of larlite, this branch uses the `BasicTool/GeoAlgo` and `SelectionTool/OpT0Finder` modules in UserDev. We have to build those as well.

    $ git clone https://github.com/LArLight/larlite larlite
    $ cd larlite
    $ source config/setup.sh
    $ make
    [... if everything builds ...]
    $ cd UserDev/BasicTool
    $ make
    $ cd ../SelectionTool/OpT0Finder
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

## Running the code

Each portion of the cosmic muon tagger/contained ROI selection has its own binary

  * ThruMu tagger: `app/ThruMu/bin/thrumu`
  * StopMu tagger: `app/StopMu/bin/run_stopmu2`
  * Contained ROI selection: `app/ContainedROI/bin/run_roi_selection`

You will find READMEs for each (describing how to run them) in their respective `bin` folders.

## Notes: To do/where things can improve

### ThruMu

* first pass for through-going muons should rely on simple straight line (in 3D) fitter. Pairs of end points and the pixels between them should be tagged.  Then A* algorithm can work on images where obvious straight through-going muons have been removed [done.]
* path restriction to avoid wildly deviant paths [done. didn't do much.]
* given path from A* or whatever. tag pixels in image. fill in gaps between nodes. [done.]
* improve post-processing for pass0. check if connected endpoints are really end points and not point in middle of track.  do this by extending past ends and seeing if charge around. [this turned out not to be so useful as endpt often tagged inside track.]
* can we restrict a* to those which have good points near the start/end points
* loosen goodpoint definition in linear3d track to 2 charge only or ( m badch, 3-m charge planes) [done. helped complete long tracks]
* If we can ID 2 more space-points near the line, we can solve for the value of the control points
* Need a post-processor step for linear3d tracks which removes duplicates [done]
* Need post-processor for A-star tracks.
    * Removes tracks which are duplicates of previous linear3d tracks.
    * if post-processor can ID bad tracks, that would be ideal good.
    * Maybe find sharp turns.
    * Fix bulges.

### StopMu

* stop-mu tagger gets stuck while trying to cross thru-mu tagged pixels [addressed]
* an idea is to iteratively find extrema points on the cluster and then use the 3D A* to fit path [done]
* if 3D A* applied, then can have 3D space points to make better flash hypothesis for matching [done]

### Contained ROI Selection

* We want to incorporate 3D constraints on cluster location. This is a bit harder because we do not have information like a good, trustworthy 3D spacepoint or aprior knowledge of the cluster shape.
* Approach is to cluster separately on each plane. Also probably want to form cluster groups as well
* We then match clusters across the planes roughly in time
* We then define the overlap time range.
* Next, we break up the overlap time range into consecutive chunks. For each chunk, we get the wire range of the cluster group.
* For each time chunk, we define the overlap region in (Y,Z).  Since this occurs over some time chunk, we form a polygon in 3D that represents the 3-plane consistent space
* We break up exclude portions of the cluster group inconsistent with the 3-plane overlap polygon.
* Adding up the time chunks, we should have a 3D volume representing the bounding polygon.
* We can do simplistic recon. at this stage or maybe form 'approx' charge distribution by setting charge at centroid of time-chunk volume.
* We want approx. 3D distribution of charge to build flash-hypothesis to test against in-time flash
* We also can use 3D volume to test if near thrumu or stopmu tracks.  This can remove brem photons.


### Shower Reco Variables
A description of the  variables made by the shower reco code
Reco variables:
* shower1_E_U, shower1_E_V, shower1_E_Y = Energy of the first shower in each plane (MeV) (y used in analysis)
* shower2_E_U, shower2_E_V, shower2_E_Y = Energy of the first shower in each plane (MeV) (y used in analysis)
* shower1_gap_U, shower1_gap_V, shower1_gap_Y = size of gap (cm) between the 1st shower start and vertex (NOT FILLED PROPERLY??)
* shower2_gap_U, shower2_gap_V, shower2_gap_Y = size of gap (cm) between the 2nd shower start and vertex (NOT FILLED PROPERLY??)
* shower1_dir_3d_X, shower1_dir_3d_Y, shower1_dir_3d_Z = the components of the 1st 3d shower direction unit vector
* shower2_dir_3d_X, shower2_dir_3d_Y, shower2_dir_3d_Z = the components of the 2nd 3d shower direction unit vector
* shower1_dir_2d_U, shower1_dir_2d_V, shower1_dir_2d_Y = Direction w.r.t wire direction (+x) of first reco shower in each plane
* shower2_dir_2d_U, shower2_dir_2d_V, shower2_dir_2d_Y = Direction w.r.t wire direction (+x) of second reco shower in each plane
* shower1_op_2d_U, shower1_op_2d_V, shower1_op_2d_Y = Opening angle of 2d shower first reco shower in each plane
* shower2_op_2d_U, shower2_op_2d_V, shower2_op_2d_Y = Opening angle of 2d shower second reco shower in each plane
* shower1_start_2d_U_X, shower1_start_2d_U_Y, shower1_start_2d_U_Z = first reco shower, u plane, X = column, Y = row, in standard dl pixels. (Z is empty, as it should be. Accidentally created) 
* shower1_start_2d_V_X, shower1_start_2d_V_Y, shower1_start_2d_V_Z = first reco shower, v plane, X = column, Y = row, in standard dl pixels. (Z is empty, as it should be. Accidentally created) 
* shower1_start_2d_Y_X, shower1_start_2d_Y_Y, shower1_start_2d_Y_Z = first reco shower, y plane, X = column, Y = row, in standard dl pixels. (Z is empty, as it should be. Accidentally created)  
* shower2_start_2d_U_X, shower2_start_2d_U_Y, shower2_start_2d_U_Z = second reco shower, u plane, X = column, Y = row, in standard dl pixels. (Z is empty, as it should be. Accidentally created) 
* shower2_start_2d_V_X, shower2_start_2d_V_Y, shower2_start_2d_V_Z = second reco shower, y plane, X = column, Y = row, in standard dl pixels. (Z is empty, as it should be. Accidentally created)
* shower2_start_2d_Y_X, shower2_start_2d_Y_Y, shower2_start_2d_Y_Z =  second reco shower, y plane, X = column, Y = row, in standard dl pixels. (Z is empty, as it should be. Accidentally created)  
* shower1_impact = from 3d reco, distance of pca line to vertex at closest point, for first reco shower
* shower2_impact = from 3d reco, distance of pca line to vertex at closest point, for second reco shower
* shower1_smallQ_U, shower1_smallQ_V, shower1_smallQ_Y = uses old method of crops to find enclosed ADC (only done for 1 shower case) 
* shower1_sumQ_U, shower1_sumQ_V, shower1_sumQ_Y = Sum of ADC from ssnet shower pixels in first shower in each plane
* shower2_sumQ_U, shower2_sumQ_V, shower2_sumQ_Y = Sum of ADC from ssnet shower pixels in second shower in each plane

Truth Variables:
* haspi0 = 1 if there is a pi0 in the event
* ccnc = 1 if cc,0 if nc
* truefid = 1 if the true vtx is in fiducial
* numshowers = number of true showers
* firstdirection_true_X,firstdirection_true_Y,firstdirection_true_Z = 3d unit vector of the true shower matched to the first reco shower
* seconddirection_true_X,seconddirection_true_Y,seconddirection_true_Z = 3d unit vector of the true shower matched to the second reco shower
* shower_start_2d_true_row = 2d row (time direction) value of the 2d start of the true shower matched to first reco shower
* shower_start_2d_true_colU, shower_start_2d_true_colV, shower_start_2d_true_colY = 2d col (wire direction) value of the 2d start of the true shower matched to first reco shower in each plane
* secondshower_start_2d_true_row = 2d row (time direction) value of the 2d start of the true shower matched to second reco shower
* secondshower_start_2d_true_colU, secondshower_start_2d_true_colV, secondshower_start_2d_true_colY = 2d col (wire direction) value of the 2d start of the true shower matched to second reco shower in each plane
* shower_energy_true = true energy of shower matched to first reco shower
* secondshower_energy_true = true energy of shower matched to second reco shower
* shower_recotrue_dist = 2d dist(yplane) between shower starts of true shower matched to first reco shower and first reco shower
* secondshower_recotrue_dist = 2d dist(yplane) between shower starts of true shower matched to second reco shower and second reco shower
