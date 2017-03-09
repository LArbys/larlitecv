# Running the Through-going Muon Tagger

## Inputs

One needs a larlite and larcv file.

The larlite file needs to contain:

  * Set of opflash products. the code has been run by the output of both SimpleFlash and OpFlash
  * ChStatus product which provides the list of badchs (note: can use larcv version of product instead)

The LArCV file needs to contain:

  * TPC images
  * ChStatus product which provides the list of badchs (note: can use larcv version of product instead)

## Configuration File

An example configuration file can be found in `bmt.cfg`.

### Parameters to setup before running

  * Set `InputLArCVImages` to the name of the TPC `Image2D` branch
  * Set the list of opflash tree names in `BMTFlashTagger:OpflashProducers` (e.g. `["simpleFlashBeam","simpleFlashCosmic"]`)
  * The output larcv file name is set in `IOManager:OutFileName`
  * The output larlite file name is set in `StorageManager:OutFileName`
  * The source of the ChStatus data. Set `ChStatusDataType`. Choose either "LARCV" or "LARLITE".

### How do I find the name of the trees?

To find out what values one should use, one can check the ROOT files.

For LARCV, one might find the following:

    root -b output_nueccqe_larcv.root 
    root [0]
    Attaching file output_nueccqe_larcv.root as _file0...
    (TFile *) 0x55fc3d0f4100
    root [1] .ls
    TFile**	output_nueccqe_larcv.root	
     TFile*		output_nueccqe_larcv.root	
       KEY: TTree	image2d_comb_tpc_tree;8	comb_tpc tree
       KEY: TTree	image2d_comb_tpc_tree;7	comb_tpc tree
       KEY: TTree	image2d_seg_comb_tpc_tree;6	seg_comb_tpc tree
       KEY: TTree	image2d_seg_comb_tpc_tree;5	seg_comb_tpc tree
       KEY: TTree	partroi_comb_tpc_tree;1	comb_tpc tree
       KEY: TTree	chstatus_comb_tpc_tree;1	comb_tpc tree
       KEY: TTree	image2d_comb_pmt_tree;1	comb_pmt tree


What we are looking for is the name of the tree that contains our `Image2D` objects.  The patterns for the tree names goes as

    [data product type]_[tree_name_w_underscores]_tree
    
In this case we want to pass along the images in the `comb_tree` tree (from `image2d_comb_tpc_tree`).

For Larlite, we want to find the names of the opflash and channel status trees.  The trees have a similar naming pattern in larlite.

    root -b output_nueccqe_larlite_full.root 
    root [0] 
    Attaching file output_nueccqe_larlite_full.root as _file0...
    (TFile *) 0x556abb096a40
    root [1] .ls
    TFile**		output_nueccqe_larlite_full.root	
      TFile*		output_nueccqe_larlite_full.root	
        KEY: TTree	wire_caldata_tree;11	wire Tree by caldata
        KEY: TTree	wire_caldata_tree;10	wire Tree by caldata
        KEY: TTree	simphotons_largeant_tree;2	simphotons Tree by largeant
        KEY: TTree	simphotons_largeant_tree;1	simphotons Tree by largeant
        ...
        KEY: TTree	ophit_ophitSat_tree;1	ophit Tree by ophitSat
        KEY: TTree	opflash_opflash_tree;1	opflash Tree by opflash
        KEY: TTree	opflash_opflashSat_tree;1	opflash Tree by opflashSat
        KEY: TTree	trigger_triggersim_tree;1	trigger Tree by triggersim
        ...
        KEY: TTree	chstatus_chstatus_tree;1	chstatus Tree by chstatus
        KEY: TTree	potsummary_generator_tree;1	potsummary Tree by generator

For the opflash trees, we would want `opflashSat`. (Note: `opflash` in this case are a repeat of flashes using the non-saturated corrected PMT waveforms. We ignore this tree.) For the channel status tree, its name is `chstatus`.


## Running the program

First make a text file containing a list of paths to larcv and larlite files.
For example, if one wants to pass the files `my_larcv_data.root` and `my_larlite_data.root`, one way to make this files is to run

    echo $PWD/my_larcv_data.root > flist_larcv.txt
    echo $PWD/my_larlite_data.root > flist_larlite.txt

Then to run the program using some config (e.g. `bmt.cfg`) run (from this folder)

    ./thrumu bmt.cfg flist_larcv.txt flist_larlite.txt


Note that the file lists can contain more than one file.

To control which block of events are run you can set the following parameters in the configuration file

    StartEntry: 0
    NumEntries: 10

## Test data samples

### [MCC7]

  * Sample of MC cosmic + MC nue between [200,600] MeV (on the FNAL system)

      /uboone/data/users/tmw/dl/test_samples/premcc8_nue

  * Sample of MC cosmic + MC numu between [200,600] MeV (on the FNAL system)

      /uboone/data/users/tmw/dl/test_samples/premcc8_numu

## Other Parameter Descriptions

(to do!)

