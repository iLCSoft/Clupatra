#ifndef KalTrack_h
#define KalTrack_h


#include <ostream>
#include <map>
#include <list>
#include <vector>

#include <gearimpl/Vectors.h>

class TVector3 ;
class TObjArray ;
class TKalTrack ;
class TKalDetCradle ;

namespace IMPL{
  class TrackImpl ;
}


class KalTrack ;

std::ostream& operator<<(std::ostream& o, const KalTrack& trk) ; 


//typedef std::list< std::pair< int, gear::Vector3D> > PointList ;
typedef std::vector< gear::Vector3D* > PointList ;


/** Track class for pat rec and track fitting with KalTest */
                                                                                
class KalTrack { // : public TKalTrack {

  friend std::ostream& operator<<(std::ostream& o, const KalTrack& trk) ;

private:
  KalTrack() {} // no default c'tor
  
public:

  // /** C'tor - initializing with size*/
  // KalTrack(Int_t n = 1) : TKalTrack(n) {}
  
  /** C'tor - initiale with detector */
  KalTrack(TKalDetCradle* det) ;

  ~KalTrack() ;


  /**Calculate crossing points with next layers at end of track segment */
  void findXingPoints() ;

  const PointList& getXingPoints() { return _xingPts ; }
 
  gear::Vector3D* getXingPointForLayer(int l) { return _xingPts.at( l ) ; }
  
  /** Add a faked hit to get track state at the IP */
  void addIPHit() ;
  
  
  template <class T>
  void setCluster(T* cluster) {  _cluster = cluster ;  }

  template <class T>
  T* getCluster() { return static_cast<T*>( _cluster ) ; } 

  /** template for adding hits of any type from a container - user needs to provide functor 
      classes for extracting 
      TVector3 position
      int layer.
  */
  template <class In, class Pos, class Layer > 
  void addHits( In first, In last, Pos position, Layer layer) {
    
    while( first != last ){
      
      int l = layer( *first ) ;
      TVector3 pos = position( *first ) ;
      
      addHit( pos ,  l ) ;
      
      ++first ;
    }

    //    addIPHit() ;
  }
  
  void fitTrack() ;

  void toLCIOTrack( IMPL::TrackImpl* trk) ; 

  //void fitTrack(IMPL::TrackImpl* trk) ;
  
  
protected:

  /** add a hit at position in layer */
  void addHit( const TVector3& pos, int layer ) ;
  
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
