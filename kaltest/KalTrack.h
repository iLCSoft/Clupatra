#ifndef KalTrack_h
#define KalTrack_h


#include <ostream>
#include <map>
#include <list>
#include <vector>

//#include "TVector3.h"

#include <gearimpl/Vectors.h>
//#include "LCRTRelations.h"
#include "EVENT/LCObject.h"

#include "TVTrackHit.h"
#include "TObjArray.h"

class TVector3 ;
//class TObjArray ;
class TKalTrack ;
class TKalDetCradle ;
class TKalTrackState ;

namespace IMPL{
  class TrackImpl ;
}
namespace EVENT{
  class TrackerHit ;
}

class KalTrack ;

std::ostream& operator<<(std::ostream& o, const KalTrack& trk) ; 


//typedef std::list< std::pair< int, gear::Vector3D> > PointList ;
typedef std::vector< gear::Vector3D* > PointList ;



struct KalTrackLink : EVENT::LCExtension< KalTrackLink, KalTrack> {} ;



/** Track class for pat rec and track fitting with KalTest */
                                                                                
class KalTrack { // : public TKalTrack {

  friend std::ostream& operator<<(std::ostream& o, const KalTrack& trk) ;

private:
  KalTrack() {} // no default c'tor
  
public:

  // /** C'tor - initializing with size*/
  // KalTrack(Int_t n = 1) : TKalTrack(n) {}
  
  /** C'tor - initiale with detector */
  KalTrack(TKalDetCradle* det, bool ownHits=false ) ;

  ~KalTrack() ;


  /**Calculate crossing points with next layers at end of track segment */
  void findXingPoints() ;

  //  const PointList& getXingPoints() { return _xingPts ; }
  PointList& getXingPoints() { return _xingPts ; }
 
  gear::Vector3D* getXingPointForLayer(int l) { return _xingPts.at( l ) ; }
  

  /** Find the next crossing point for the current Track - step: #layers from last hit used in fit */
  void findNextXingPoint(gear::Vector3D& v, int& layer, int step,  bool backward=false)  ;
  

  /** Add a faked hit to get track state at the IP */
  void addIPHit() ;


  /** omega of the track (1./R) */
  double getOmega() ; 

  /** Compute the chi2 ( (a0-a1)^t * cov0^-1 * (a0-a1) ) between two track states - using the covariance matrix of the first */ 
  static double chi2( const KalTrack& t0 , const KalTrack& t1) ;


  /** current track state */
  TKalTrackState& getTrackState() const ;
 
  /** Number of hits on this track */
  unsigned getNHits() const ;


  template <class T>
  void setCluster(T* cluster) {  _cluster = cluster ;  }

  template <class T>
  T* getCluster() { return static_cast<T*>( _cluster ) ; } 

  /** template for adding hits of any type from a container - user needs to provide functor 
      classes for extracting 
      TVector3 position
      int layer.
  */
  template <class In, class Pos, class Layer, class LCIOHit> 
  void addHits( In first, In last, Pos position, Layer layer, LCIOHit lcioHit) {
    

    while( first != last ){
      
      int l = layer( *first ) ;
      TVector3 pos = position( *first ) ;
      
      addHit( pos ,  l , lcioHit( *first ) ) ;

      //      std::cout << " ------- addHit at layer " <<  l << std::endl ; //<< " r: " << pos.Perp() << std::endl ;
      
      ++first ;
    }

    //    addIPHit() ;
  }
  
  /** Add a TVTrackHit - created by the user. To be used with KalTest::CFG::ownsHits = false !
   *  Hits will be used in fit in the order in which they where added with this method. 
   */
  void addHit( TVTrackHit* hit){

    _kalHits->Add( hit ) ;
    //    std::cout <<  "  added hit to KalTrack (" << this << ") - # hits = " << _kalHits->GetEntriesFast() << std::endl ; 
  }

  /** Add the hit and filter/update the track state - will fail if the current track does not reach the surface or 
   *  the deltaChi2 is larger than maxDeltaChi2 (if given).
   */
  bool addAndFilter( TVTrackHit* hit, double maxDeltaChi2 = -1.) ;


  /** Compute the delta Chi2 that would result from adding the hit - the hit is not added and the track state is not changed.
   */
  double testDeltaChi2( TVTrackHit* hit ) ;
    

  /** Fit the hits added with addHits in the specified order (KalTest::FitForward/KalTestFitBackward) */
  void fitTrack( bool fitOrder ) ;

  void toLCIOTrack( IMPL::TrackImpl* trk) ; 

  //void fitTrack(IMPL::TrackImpl* trk) ;
  
  
protected:

  /** add a hit at position in layer */
  void addHit( const TVector3& pos, int layer, EVENT::TrackerHit* hit ) ;
  
  /** clear all hits */
  //  void clear() ; { _kalHits->Clear() ; }

protected:
  void init() ;

  TKalDetCradle* _det ;     // the detector cradle
  TObjArray* _kalHits;      // array to store hits
  TKalTrack* _trk ;         // the KalTest track
  PointList  _xingPts ;     // crossing points with other layers (w/o) hits
  void* _cluster ;          // generic pointer to a cluster of hits

};

#endif
