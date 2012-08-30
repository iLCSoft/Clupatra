#ifndef TrackEfficiency_h
#define TrackEfficiency_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <vector>
#include <string>


using namespace lcio ;
using namespace marlin ;

class TH1 ;

/** Analysis plots for tracking efficiency.
 * 
 *  @author F.Gaede, DESY
 *  @version $Id$ 
 */
class TrackEfficiency : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new TrackEfficiency ; }
  
  
  TrackEfficiency() ;
  
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
  std::string _sthColName ;
  std::string _t2mColName ;
  std::string _m2tColName ;
  std::string _trkColName ;

  std::vector< TH1* > _h1 ;
  
  FloatVec _ptRange ;

  int _nRun ;
  int _nEvt ;

} ;

#endif



