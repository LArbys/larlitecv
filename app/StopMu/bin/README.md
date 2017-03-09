# Running the Stopping Muon Tagger

## Inputs

One needs a larcv file. It uses no larlite data.

The LArCV file needs to contain:

  * `Image2D` for the TPC event image (often passed on by the thrumu code)
  * `Image2D` whose pixels are tagged by the thrumu tagger (product of `thrumu`)
  * `ChStatus` product which provides the list of badchs
  * Sets of `Pixel2DCluster` which stores the unused Boundary Space Points (product of `thrumu`)

## Configuration File

An example configuration file can be found in `smc.cfg`.

### Parameters to setup before running

  * `InputLArCVImages`: name of tree that contains the TPC images
  * `LArCVInputList`: text file with list of paths to larcv files
  * `LArLiteInputList:`: test file with list of paths to larlite files (not really needed)
  * `EndPointProducers`: List of trees which store the unused boundary end points
  * `EndPointProducers`: List of trees which store the unused boundary end points
  * `ThrumuLabeledImagesProducer`: List of trees that contained images where pixels have been tagged
  * `IOManager:OutFileName`: output larcv file name is set in 
  * `StorageManager:OutFileName`: output larlite file name is set in 

### How do I find the name of the trees?

To find out what values one should use, one can check the ROOT files.

Here is an example file produced by the through-going muon tagger application `thrumu`:

    root -b output_thrumu_larcv.root 
    root [0]
    Attaching file output_thrumu_larcv.root as _file0...
    (TFile *) 0x55fc3d0f4100
    root [1] .ls
    TFile**	output_thrumu_larcv.root
     TFile*		output_thrumu_larcv.root
     KEY: TTree		partroi_tpc_tree;1	tpc tree
     KEY: TTree		chstatus_tpc_tree;1	tpc tree
     KEY: TTree		image2d_modimgs_tree;1	modimgs tree
     KEY: TTree		image2d_marked3d_tree;1	marked3d tree
     KEY: TTree		image2d_realspacehits_tree;1	 realspacehits tree
     KEY: TTree		pixel2d_anodepts_tree;1		 anodepts tree
     KEY: TTree		pixel2d_cathodepts_tree;1	cathodepts tree
     KEY: TTree		pixel2d_imgendpts_tree;1	imgendpts tree
     KEY: TTree		pixel2d_thrumu2d_tree;1	thrumu2d tree
     KEY: TTree		pixel2d_topspacepts_tree;1	topspacepts tree
     KEY: TTree		pixel2d_botspacepts_tree;1	botspacepts tree
     KEY: TTree		pixel2d_upspacepts_tree;1	upspacepts tree
     KEY: TTree		pixel2d_downspacepts_tree;1	downspacepts tree
     KEY: TTree		image2d_badchs_tree;1		badchs tree
     KEY: TTree		image2d_gapchs_tree;1		gapchs tree
     KEY: TTree		image2d_boundarypixels_tree;1	boundarypixels tree
     KEY: TTree		pixel2d_unused_topspacepts_tree;1	       unused_topspacepts tree
     KEY: TTree		pixel2d_unused_botspacepts_tree;1	       unused_botspacepts tree
     KEY: TTree	 	pixel2d_unused_upspacepts_tree;1	       unused_upspacepts tree
     KEY: TTree	 	pixel2d_unused_downspacepts_tree;1	       unused_downspacepts tree
     KEY: TTree	 	pixel2d_unused_anodepts_tree;1		       unused_anodepts tree
     KEY: TTree	 	pixel2d_unused_cathodepts_tree;1	       unused_cathodepts tree
     KEY: TTree		pixel2d_unused_imgendpts_tree;1	unused_imgendpts tree

What we are looking for is the name of the tree that contains our `Image2D` objects.  The patterns for the tree names goes as

    [data product type]_[tree_name_w_underscores]_tree
    
In this case we need images for raw ADC values, which are the `modimgs` tree.
We also need images whose pixels mark where a through-going muon was tagged.
We then will need `marked3d`.

We also need a list of trees that contained the unused boundary points saved as a list of pixels in the images.
Here we will use

    ["unused_topspacepts","unused_botspacepts","unused_upspacepts","unused_downspacepts","unused_anodepts","unused_cathodepts","unused_imgendpts"]

## Running the program

First make a text file containing a list of paths to larcv and larlite files.
For example, if one wants to pass the files `my_larcv_data.root` and `my_larlite_data.root`, one way to make this files is to run

    echo $PWD/my_larcv_data.root > flist_larcv.txt
    echo $PWD/my_larlite_data.root > flist_larlite.txt

Remember to point to these text files in the configuration file.

Then to run the program using some config (e.g. `smc.cfg`) run (from this folder)

    ./run_stopmu2 smc.cfg


Note that the file lists can contain more than one file.

To control which block of events are run you can set the following parameters in the configuration file

    StartEntry: 0
    NumEntries: 10

In this example, entries number 0 to 9 (inclusive) will be run.

## Test data samples

### [MCC7]

  * Sample of MC cosmic + MC nue between [200,600] MeV post ThruMu (on the FNAL system)

      /uboone/data/users/tmw/dl/test_samples/premcc8_nue

  * Sample of MC cosmic + MC numu between [200,600] MeV post ThruMu (on the FNAL system)

      /uboone/data/users/tmw/dl/test_samples/premcc8_numu

## Other Parameter Descriptions

(to do!)

