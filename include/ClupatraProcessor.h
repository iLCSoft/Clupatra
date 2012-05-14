#ifndef ClupatraProcessor_h
#define ClupatraProcessor_h 1

#include "assert.h"

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>


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

/** Clustering based pattern recognition for a TPC...
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
  //  float _rCut ;
  float _duplicatePadRowFraction ;
  
  float  _dChi2Max ;
  float  _chi2Cut ;
  int    _maxStep ; 

  float _minLayerFractionWithMultiplicity ;
  int   _minLayerNumberWithMultiplicity ;
  
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
  const gear::PadRowLayout2D* _padLayout ;

//   NNClusterer* _clusterer ;

} ;

#endif



