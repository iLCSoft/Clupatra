#ifndef OuterRimSearch_h
#define OuterRimSearch_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>


using namespace lcio ;
using namespace marlin ;


class KalTest ;


/** Clustering based pattern recognition for a TPC.
 *  Start with cluster segments found in the outer rim of the TPC and then use the segements with
 *  with a Kalman fitter (KalTest) to search inwards ...
 */
class OuterRimSearch : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new OuterRimSearch ; }
  
  
  OuterRimSearch() ;
  
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

  /** Input collection name.
   */
  StringVec _colNames ;
  std::string _outColName ;

  float _distCut ;
  float _rCut ;
  float _duplicatePadRowFraction ;
  int _minCluSize ;

  int _nRun ;
  int _nEvt ;

  KalTest* _kalTest ;

//   NNClusterer* _clusterer ;

} ;

#endif



