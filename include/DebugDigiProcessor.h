#ifndef DebugDigiProcessor_h
#define DebugDigiProcessor_h 1

#include "marlin/Processor.h"
#include "marlin/EventModifier.h"
#include "lcio.h"
#include <vector>
#include <string>


using namespace lcio ;
using namespace marlin ;


/** DebugDigiProcessor puts all TrackerHits from the tracks in TrackCollection  into one or more TrackerHitCollections.
 *  The Hits are sorted into the output collections based on CellIDDecoder<T>( ILDCellID0::encoder_string )( t )[ ILDCellID0::subdet ] <br/>
 *  Pre-existing collections with the names of the TrackerHitCollections are renamed to NAME_OLD - as is the TrackCollection.
 *
 * @parameter TrackCollection                           Name of the input collections with  Tracks
 * @parameter TrackerHitCollections                     The tracker hit collection names - needs to run parallel to SubDetectorIDs
 * @parameter SubDetectorIDs                            The IDs of the subdetectors corresponding to the hit collections in TrackerHitCollectionNames
 * 
 * @author F.Gaede, DESY
 * @version $Id$ 
 */
class DebugDigiProcessor : public Processor, public EventModifier {
  
 public:
  
  virtual Processor*  newProcessor() { return new DebugDigiProcessor ; }
  
  
  DebugDigiProcessor() ;
  
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

  virtual void modifyEvent( LCEvent * evt ) ;

  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  
  
  virtual const std::string & name() const { return Processor::name() ; }


 protected:

  std::string _trkColName ;

  StringVec _hitColNames ;
  IntVec    _subDetIDs ;

  int _nRun ;
  int _nEvt ;

} ;

#endif



