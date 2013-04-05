#ifndef ClupatraProcessor_h
#define ClupatraProcessor_h 1

#include "assert.h"

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>

#include "gear/TPCModule.h"

using namespace lcio ;
using namespace marlin ;

namespace MarlinTrk{
  class IMarlinTrkSystem ;
}
namespace gear{
  class TPCParameters ;
  class PadRowLayout2D ;
}

namespace EVENT{ 
  class Track ;
}

/** ClupatraProcessor : nearest neighbour clustering seeded pattern recognition using Kalman-Filter for extrapolation and 
 *  hit search.
 * 
 *   @parameter TPCHitCollection         Name of the tpc hit input collections
 *   @parameter OutputCollection         Name of the output collection with final TPC tracks
 *   @parameter SegmentCollectionName    Name of the output collection that has the individual track segments
 *   @parameter CreateDebugCollections   optionally create some debug collection with intermediate track segments and used and unused hits
 * 
 *   @parameter DistanceCut              Cut for distance between hits in mm for the seed finding
 *   @parameter Chi2Cut                  the maximum chi2-distance for which a hit is considered for merging
 *   @parameter CosAlphaCut              Cut for max.angle between hits in consecutive layers for seed finding - NB value should be smaller than 1 - default is 0.9999999 
 *   @parameter MaxDeltaChi2             the maximum delta Chi2 after filtering for which a hit is added to a track segement
 * 
 *   @parameter DuplicatePadRowFraction  allowed fraction of hits in same pad row per track
 *   @parameter NLoopForSeeding          number of seed finding loops - every loop increases the distance cut by DistanceCut/NLoopForSeeding
 *   @parameter NumberOfZBins            number of bins in z over total length of TPC - hits from different z bins are nver merged
 *   @parameter PadRowRange              number of pad rows used in initial seed clustering
 * 
 *   @parameter MaxStepWithoutHit                 the maximum number of layers without finding a hit before hit search search is stopped 
 *   @parameter MinLayerFractionWithMultiplicity  minimum fraction of layers that have a given multiplicity, when forcing a cluster into sub clusters
 *   @parameter MinLayerNumberWithMultiplicity    minimum number of layers that have a given multiplicity, when forcing a cluster into sub clusters
 *   @parameter MinimumClusterSize                minimum number of hits per cluster
 * 
 *   @parameter TrackEndsOuterCentralDist maximum radial distance [mm] from outer field cage of last hit, such that the track is considered to end at the end 
 *   @parameter TrackEndsOuterForwardDist maximum distance in z [mm] from endplate of last hit, such that the track is considered to end at the end 
 *   @parameter TrackIsCurlerOmega        minimum curvature omega of a track segment for being considered a curler
 *   @parameter TrackStartsInnerDist      maximum radial distance [mm] from inner field cage of first hit, such that the track is considered to start at the beginning 
 * 
 *   @parameter EnergyLossOn             Use Energy Loss in Fit
 *   @parameter MultipleScatteringOn     Use MultipleScattering in Fit
 *   @parameter SmoothOn                 Smooth All Mesurement Sites in Fit
 * 
 *   @parameter pickUpSiHits             try to pick up hits from Si-trackers
 *   @parameter SITHitCollection         name of the SIT hit collections - used to extend TPC tracks if (pickUpSiHits==true)
 *   @parameter VXDHitCollection         name of the VXD hit collections - used to extend TPC tracks if (pickUpSiHits==true)
 * 
 *   @parameter Verbosity               verbosity level of this processor ("DEBUG0-4,MESSAGE0-4,WARNING0-4,ERROR0-4,SILENT")
 * 
 * @author F.Gaede, DESY, 2011/2012
 * @version $Id$ 
 */
class ClupatraProcessor : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new ClupatraProcessor ; }
  
  
  ClupatraProcessor() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ;

  //  void doBLA() ;
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  
  
 protected:

  /** helper method to compute a few track segment parameters (start and end points, z spread,...) 
   */
  void computeTrackInfo(  lcio::Track* lTrk  ) ;


  void pickUpSiTrackerHits( EVENT::LCCollection* trackCol , LCEvent* evt) ;

  /** Input collection name.
   */
  std::string _colName ;
  std::string _vxdColName ;
  std::string _sitColName ;
  std::string _outColName ;
  std::string  _segmentsOutColName ;

  float _distCut ;
  float _cosAlphaCut ;

  float _duplicatePadRowFraction ;
  float _bfield ;
  
  float  _dChi2Max ;
  float  _chi2Cut ;
  int    _maxStep ; 

  float _minLayerFractionWithMultiplicity ;
  int   _minLayerNumberWithMultiplicity ;
  int   _nLoop ;
  
  float _trackStartsInnerDist ;
  float _trackEndsOuterCentralDist ;
  float _trackEndsOuterForwardDist ;
  float _trackIsCurlerOmega ;
  
  int   _minCluSize ;
  int   _padRowRange ; 
  int   _nZBins ;

  bool _MSOn ;
  bool _ElossOn ;
  bool _SmoothOn ;
  bool _pickUpSiHits ;

  bool _createDebugCollections ;

  int _nRun ;
  int _nEvt ;

  MarlinTrk::IMarlinTrkSystem* _trksystem ;

  const gear::TPCParameters*  _gearTPC ;

//   NNClusterer* _clusterer ;

} ;

#endif



