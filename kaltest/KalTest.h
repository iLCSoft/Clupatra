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

namespace IMPL{
  class TrackImpl ;
}


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
  
/** Enums for identifying detectors
 */
  typedef enum {
    unknown = 0 ,
    VXD =  1 ,
    SIT =  2 ,
    TPC =  3 ,
    SET =  4 ,
    ETD =  5 ,
    FTD =  6 , // add new detecors here
    DetectorID_Size ,
    DetectorID_Factor = 10000 
  } DetectorID ;



  /** Default c'tor, initializes the geometry from GEAR. */
  KalTest( const gear::GearMgr& gearMgr) ;
  
  ~KalTest() ;
  
  KalTrack* createKalTrack()  { return new KalTrack( _det ) ; }

  int indexOfFirstLayer( DetectorID det) ;

protected:
  void init() ;

  const gear::GearMgr* _gearMgr ;

  TKalDetCradle* _det ;            // the detector cradle

  std::vector<int> _idOffsets ;

} ;

#endif
