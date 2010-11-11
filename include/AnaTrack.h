#ifndef AnaTrack_h
#define AnaTrack_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <vector>
#include <string>


using namespace lcio ;
using namespace marlin ;

class TH1 ;

/** Analysis plots for Tracks.
 */
class AnaTrack : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new AnaTrack ; }
  
  
  AnaTrack() ;
  
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

  /** Input collection names.
   */

  std::string _trkColName ;
  std::string _relColName ;

  std::vector< TH1* > _h1 ;

  int _nRun ;
  int _nEvt ;

} ;

#endif



