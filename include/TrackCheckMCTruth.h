#ifndef TrackCheckMCTruth_h
#define TrackCheckMCTruth_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <vector>
#include <string>


using namespace lcio ;
using namespace marlin ;

class TH1 ;

/** Check a track collection based on MCTruth for merged and split tracks.
 */
class TrackCheckMCTruth : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new TrackCheckMCTruth ; }
  
  
  TrackCheckMCTruth() ;
  
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

  std::string _mcpColName ;

  StringVec _sthColNames ;

  std::string _relTrkMCPColName ;
  std::string _relMCPTrkColName ;
  
  std::string _trkColName ;
  
  std::string _mergedColName ;
  std::string _splitColName ;

  float _hitFractionSplit ;
  float _hitFractionMerged ;

  FloatVec _ptRange ;

  std::vector< TH1* > _h1 ;
  
  int _nRun ;
  int _nEvt ;

} ;

#endif



