# larlitecv
Analysis processor framework for working with LArLite and LArCV data

This branch contains the muon-tagger code.

* app/ThruMu: through-going muon tagger using pairs of boundary end points
* app/StopMu: tags stopping muon tracks using single boundary end points. relies on output of ThruMu code
* app/ContainedROI: finds contained ROI tracks. uses output of StopMu code


## Running the code

Each portion of the cosmic muon tagger/contained ROI selection has its own binary

* ThruMu tagger: `app/ThruMu/bin/thrumu`
* StopMu tagger: `app/StopMu/bin/run_stopmu2`
* Contained ROI selection: `app/ContainedROI/bin/run_roi_selection`

You will find READMEs for each (describing how to run them) in their respective `bin` folders.

## Notes: To Do/where things can improve

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
