# v01-03-01

* 2024-01-17 tmadlener ([PR#19](https://github.com/iLCSoft/Clupatra/pull/19))
  - Migrate to github actions for CI

* 2020-04-13 Frank Gaede ([PR#16](https://github.com/iLCSoft/Clupatra/pull/16))
  - make compatible w/ c++17 for macos/clang
        - patch provided by K.Fujii

# v01-03

* 2018-01-31 Frank Gaede ([PR#13](https://github.com/iLCSoft/Clupatra/pull/13))
  - fix all compiler warnings (gcc54-ub1604)
  - use private inheritance for std::list in NNClusterer 
        - add wrapper functions instead

* 2018-03-13 Marko Petric ([PR#14](https://github.com/iLCSoft/Clupatra/pull/14))
  -  Fix for iLCSoft/LCIO#35

* 2018-03-28 Marko Petric ([PR#15](https://github.com/iLCSoft/Clupatra/pull/15))
  - Fix for the removal of DDSurfaces which have been merged into DDRec 
    -  includes from `DDSurfaces` -> `DDRec`
    - namespace `DDSurfaces` -> `dd4hep::rec`

# v01-02

* 2017-10-12 Frank Gaede ([PR#12](https://github.com/iLCSoft/Clupatra/pull/12))
  - set correct MarlinTrkSystem configuration for every event with `MarlinTrk::TrkSysConfig`

# v01-01

* 2017-06-20 Andre Sailer ([PR#11](https://github.com/iLCSoft/Clupatra/pull/11))
  - Adapt to changes in namespaces and LCDD -->  Detector

* 2017-05-18 Frank Gaede ([PR#10](https://github.com/iLCSoft/Clupatra/pull/10))
  - remove last leftover dependencies to Gear
  - remove any dependency on MarlinKalTest

# v01-00

* 2017-03-30 Andre Sailer ([PR#7](https://github.com/iLCSoft/Clupatra/pull/7))
  - Wrap static ILDDecoder helper object in function to avoid static initialisation at library loading time

* 2017-03-30 Frank Gaede ([PR#6](https://github.com/iLCSoft/Clupatra/pull/6))
  - add units to B-field from DD4hep
        - previously B-field was taken to be zero (10e-13)

* 2017-01-09 Frank Gaede 
  - replaced GEAR with DD4hep
    - use DD4hep::DDRec::FixedPadSizeTPCData for gear::TPCData
    - use DDSurfaces::Vector3D for gear::Vector3D
    - still use gear::Vector3D for MarlinTrk intersection calculations until it is chaned there
  - fix warning from clang ( mostly shadowing variables )  

# v00-14
F. Gaede
- made compatible with c++11
- removed -ansi -pedantic -Wno-long-long
- fixed template parameters in nnclu::(not)inRange<int,int>

# v00-13
- fix for TS @calo: set z0 to 0
- protect against failure in adding a trackerhit in TrackSegmentMerger
- make track state at calo again compatible w/ old sim/rec
- fixed creation of track state at calo
- activate propagete in track fit (-> needed for pull distributions)

# v00-12
- added parameters CaloFaceBarrelID/CaloFaceEndcapID to ClupatraProcessor ( needed for lcgeo based simulations ) default value for both : ILDDetID::Ecal
- adopted LCIOTrackConverter to use these paramters

# v00-11
- no change in CLupatra algorithm ...
- add option to use DDKalTest
- parameter "TrackSystemName"
- modified FixCellIDs: added new parameters
  - FTDCollectionName: set cellIDs for FTD
  - SetSideForVXD: set side to +-1 if true
  - LayerOffset: offset applied to layer number (ZPalanar)
  - needed for compatibility between gear and dd4hep
  - added SET hits ( used for debugging mainly)
  - added debug printout for cellIDs
  - protect against missing gear parameters 
  - added parameter "SetSideForTPC" to allow for setting the side to +/- 1
- minor bug fix in ClupatraProcessor (in DEBUG printout)
- remove FixCellIDs_Errors_TestBeam : 
  - moved to MarlinTPC

# 00-10
- bug fix: remove sqrt from computation of r_phi_res^2 for finding best hit
- included multi-module TPC support via GEAR interface
- ignore hits that are outside of the TPC range (possible in data)
- added new processor FixCellIDs_Error_TestBeam (based on the original FixCellIDs) that sets the correct hit CellIDs and hit uncertainties for test beam data
- added examples for test beam data. Look at examples/README for more details
- current version of algorithm ( May, 2012) ClupatraProcessor.cc
- fill histos for track finding efficiencies TrackEfficiency.cc
- update cellIDs in plder file: FixCellIDs.cc

# v00-09-01
- in TrackCircleDistance::operator()(): adapted to fix in KalTest where curler segments now have correct Z0 (typically close to the IP)
- do not merge segements if the z-positions of both their first hits is within 20 mm of the IP
- changed  LCIOTrackConverter::operator() to use the first constrained fit ( ~= 3rd hit) smoothed to extrapolate back to the calorimeter
- fixed issue in ClupatraProcessor::check() when no  TPCTrackerHits are present ( using LCIterator ) 
- moved local LCIterator to LCIO as UTIL::LCIterator use this in all processors
- updated to use new PlanarDigiProcessor (was: SimplePlanarTestDigiprocessor )
- fixed coverity issues: mostly unitialized c'tor, some unchecked dynamic_casts and division by zero (nhit++)
- deactivate pick'n save feature (#define WRITE_PICKED_DEBUG_TRACKS false)

# v00-09
- apply a refit with larger maximum chi2 increment to tracks that have only a small number of hits used (IMarlinTrkFitter::operator())
- drop poor seed clusters with no hits added except in the very forward (inner) region
- clupatra_new: added optional  argument TrkSystem* to addHitsAndFilter();if given, the smoothed TrackState is used for backward extrapolation; => not currently used, i.e. the last filtered  TrackStae is used ( this should give better results than the previous version that used the filtered TrackState at the 4th hit )
- don't set quality bit on tracker hit: ILDTrkHitQualityBit::USED_IN_FIT anymore, as this won't work on pre-digitized files (and is not unique wrt other pat rec code)
 - added DebugDigiProcessor : takes TrackerHits from a collection of Tracks (e.g. previously written as ClupatraPoorTrackCollection ) and writes them into new Collections ( according to subdetID ) -> see ./example/clupa_debug.xml
-  added debug method printAndSaveTrack() to Clupatra: allows to store a collection of Tracks selected through picking in CED
-  added doxygen documentation of parameters in ClupatraProcessor.h
-  Updated for new lcio TrackState copy constructor taking const reference.

# v00-08
- changed order of hits in track and hits in fit in subdetector hit numbers
- changed way how curler tracks are written out:
  - store only complete track segments - as used in the fit
  - for merged curlers only the primary segment (closest to IP) is stored in track output collection where the other segments are referenced through the track pointers (sorted in z) 
- no change in algorithm 

# v00-07-01
- set radius of innermost hits in Tracks

# v00-07
- added code to merge split tracks seeds - this happens sometimes if one link of a track seed is just a bit further than the current distance cut ( artefact of loops with increasing distance cuts )
- introduced parameter CosAlphaCut and merge hits from consecutive rows if their angle is smaller than acos(CosAlphaCut) -> helps to recover some stiff tracks in forward region   
- added parameter NLoopForSeeding (default 4 ) defining the number of loops for seeding with different distance cuts 
- restrict segment merging to two segments at a time in order to prevent falsly merging almost complete tracks (through a small 'common' stub)
  - repeat this procedure once in order to get also tracks split into three 

- use track state at first hit also for outward extrapolation as this is more accurate when merging segments

# v00-06
- added new TrackSegmentMerger for merging split tracks
- add ppointer to original track segments for merged curling tracks 
- added parameter SegmentCollectionName
- added steering parameters TrackStartsInnerDist, TrackEndsOuterCentralDist, TrackEndsOuterForwardDist, TrackIsCurlerOmega  to define "incomplete" track segments
- adapted pickupSiTracks to use strip hits for the SIT
- fixed issue w/ split seed tracks (occured, e.g. if hit distance is larger than current cut in the loop for one hit pair ...)
- made compatible with llvmgcc4.2 and clang 

# v00-05-01
- patch release: bug fix: protect against tracks w/o TS at firat/last hit

# v00-05
- reduce memory footprint: deal with one KalTest track at a time wherever possible (saves >2MByte/track)
- properly set the bit UTIL::ILDTrkHitQualityBit::USED_IN_FIT for all trackerhits used in fit
- set correct subdetectorHitNumbers() for used in fit
```
 ( trk->subdetectorHitNumbers()[ 2*lcio::ILDDetID::TPC - 1 ] =  usedInFIt ; 
   trk->subdetectorHitNumbers()[ 2*lcio::ILDDetID::TPC - 2 ] =  nHit ;  )
```
- algorithmic improvements:
- introduce a loop over increasing distance cuts for finding the tracks seeds -> fixdes some of the problems seen @ 3 TeV with extremely boosted jets

- use function split_multiplicity() to split up merged cluster seeds if possible ( improves efficiencies in forward region )

- added create_n_clusters and use for splitting up clusters with multiplicities 4 and 5 
- added gear and steering file for CLIC detector
- added some debugging functionality
- new track debugging processor - under developmet: TrackCheckMCTruth
- added debug method printTrackerHit() for picking 
- added optional debug collections: ClupatraPoorQualityTracks, ClupatraOuterSegments, ClupatraInnerSegments, ClupatraMiddleSegments
- improved TrackEfficiencyProcessor
- extended pt range to 500 GeV 
- introduce cut on good hit fraction (96%) (needs MarlinReco v01-01 )

# v00-04-01
 - fixed (small) memory leak in src/ClupatraProcessor.cc

# v00-04
- src/ClupatraProcessor.cc: - improved re-clustering:
  - now also loop in (larger) pad row range
  - increase maxdeltaChi2
- src/FixCellIDs.cc: - reduced verbosity level from WARNING to DEBUG5 for missing collections

# v00-03
- write out canonical track states: AtIP, AtFirst/LastHit, AtCalorimeter
- made hard coded parameters into processor parameters: MaxDeltaChi2, Chi2Cut, MaxStepWithoutHit, PadRowRange, MinLayerFractionWithMultiplicity, MinLayerNumberWithMultiplicity
- added optional picking up of Si-hits (experiemtal)
- add subdetectorHitNumbers:
```
	  2*ILDDetID::TPC-1 : hitsinFit
	  2*ILDDetID::TPC-2 : hitsinTrack
```  
- cleaned up code and removed obsolete processors

# v00-02
- NNClusterer.h : re-factoring of NNCluster code
- clearer type definitions 
- use Element and Cluster (as oposed to GenericHits/GenericClusters)
- clearer memory management
- classes PtrVector and PtrList: user defines ownership
- introduced NNClusterer class to collect cluster functions and conveneient typedef
- clupatra_new.cc: re-factoring of clupatra helper code using new NNClusterer  
- rewrote ClupatraProcessor to use clupatra_new and NNClusterer and MarlinTrk package => KalTest/KalTrack and KalDet like classes no longer needed
- old classs and code still left for comparisons ....

# v00-01
-  first tagged version 
-  first properly working of the algorithm ( ClupatraNew ) 
-  status of ILD Metting, Orsay 2011
 
