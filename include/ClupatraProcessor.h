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

  void pickUpSiTrackerHits( EVENT::LCCollection* trackCol , LCEvent* evt) ;

  /** Input collection name.
   */
  std::string _colName ;
  std::string _vxdColName ;
  std::string _sitColName ;
  std::string _outColName ;

  float _distCut ;
  float _rCut ;
  float _duplicatePadRowFraction ;
  int   _minCluSize ;
  int   _padRowRange ; 
  int   _nZBins ;

  bool _MSOn ;
  bool _ElossOn ;
  bool _SmoothOn ;
  bool _pickUpSiHits ;

  int _nRun ;
  int _nEvt ;

  MarlinTrk::IMarlinTrkSystem* _trksystem ;

//   NNClusterer* _clusterer ;

} ;

#endif



