#ifndef INCLUDE_KalTest
#define INCLUDE_KalTest 1

#include "gear/GearMgr.h"

//LCIO:
#include "UTIL/BitField64.h" 


#include "streamlog/streamlog.h"

#include "TObjArray.h"
#include "TVector3.h"

#include <cmath>
#include <vector>
#include "KalTrack.h"

class TKalDetCradle ;
class TVTrackHit ;

namespace IMPL{
  class TrackImpl ;
}

class ConfigFlags ;


/** Interface to KaltTest Kalman fitter - instantiates and holds the detector geometry.
 */
class KalTest{

  KalTest(){}
  
public:
  
  // define some configuration constants
  static const bool FitBackward   = kIterBackward ;
  static const bool FitForward    = kIterForward ;
  static const bool OrderOutgoing  = true ;
  static const bool OrderIncoming  = false ;
  static const bool PropagateToIP  = true ;
  
  
  /** 'Enums' for identifying detectors
   */
  struct DetID {
    static const int unknown = 0 ;
    static const int VXD =  1 ;
    static const int SIT =  2 ;
    static const int TPC =  3 ;
    static const int SET =  4 ;
    static const int ETD =  5 ;
    static const int FTD =  6 ; // add new detecors here
    static const int Size = 7 ;
    static const int Factor = 10000 ; 
  } ;


  /** 'Enums' for configuration options
   */
  struct CFG {
    static const unsigned  ownsHits = 0 ;
    static const unsigned  useQMS   = 1 ;
    static const unsigned  usedEdx  = 2 ;
    //---
    static const unsigned  size     = 3 ;
  } ;
  
  
  /** Default c'tor, initializes the geometry from GEAR. */
  KalTest( const gear::GearMgr& gearMgr) ;
  
  ~KalTest() ;
  
  KalTrack* createKalTrack() ; // { return new KalTrack( _det ) ; }
  
  /** Create a TVTrackHit from an lcio::TrackerHit - need lcio layer and detector ID.
   */
  TVTrackHit* createHit(EVENT::TrackerHit* hit, int layer, int detectorID) ;
  
  int indexOfFirstLayer( int detectorID ) ;
  
  int maxLayerIndex() ;
  
  void setOption(unsigned CFGOption, bool val) ;  

  void init() ;

protected:

  const gear::GearMgr* _gearMgr ;

  TKalDetCradle* _det ;            // the detector cradle

  std::vector<int> _idOffsets ;

  ConfigFlags* _cfg ;

} ;

#endif
