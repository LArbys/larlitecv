# larlitecv
Analysis processor framework for working with LArLite and LArCV data

This branch contains the muon-tagger code.

* app/ThruMu: through-going muon tagger using pairs of boundary end points
* app/StopMu: tags stopping muon tracks using single boundary end points. relies on output of ThruMu code
* app/ContainedROI: finds contained ROI tracks. uses output of StopMu code



## To Do/where things can improve

### ThruMu

* first pass for through-going muons should rely on simple straight line (in 3D) fitter. Pairs of end points and the pixels between them should be tagged.  Then A* algorithm can work on images where obvious straight through-going muons have been removed [done.]
* path restriction to avoid wildly deviant paths [done. didn't do much.]
* given path from A* or whatever. tag pixels in image. fill in gaps between nodes. [done.]
* improve post-processing for pass0. check if connected endpoints are really end points and not point in middle of track.  do this by extending past ends and seeing if charge around. [this turned out not to be so useful as endpt often tagged inside track.]
* can we restrict a* to those which have good points near the start/end points
* loosen goodpoint definition in linear3d track to 2 charge only or ( m badch, 3-m charge planes) [done. helped complete long tracks]
* If we can ID 2 more space-points near the line, we can solve for the value of the control points
* HIGHEST: need a post-processor step for linear3d tracks which removes duplicates

### StopMu

* stop-mu tagger gets stuck while trying to cross thru-mu tagged pixels
* an idea is to iteratively find extrema points on the cluster and then use the 3D A* to fit path
* if 3D A* applied, then can have 3D space points to make better flash hypothesis for matching

### Contained ROI Selection

* filter out contained clusters that are near muons. tag these as deltas
