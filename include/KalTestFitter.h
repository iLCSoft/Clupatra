#ifndef KalTestFitter_h
#define KalTestFitter_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>


using namespace lcio ;
using namespace marlin ;


class KalTest ;


/** Apply KalTest Kalman Fitter to a given lcio::Track collection 
 */
class KalTestFitter : public Processor {
  
public:
  
  virtual Processor*  newProcessor() { return new KalTestFitter ; }
  
  
  KalTestFitter() ;
  
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

  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  
  
 protected:

  /** Input collection name.
   */
  std::string _inColName ;
  std::string _outColName ;

  int _nRun ;
  int _nEvt ;

  KalTest* _kalTest ;

} ;

#endif



