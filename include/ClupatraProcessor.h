#ifndef ClupatraProcessor_h
#define ClupatraProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>


using namespace lcio ;
using namespace marlin ;


class KalTest ;


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



