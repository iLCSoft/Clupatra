#include "OuterRimSearch.h"

#include <time.h>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

#include <math.h>
#include <cmath>

//---- MarlinUtil 
#include "NNClusters.h"
#include "ClusterShapes.h"
#include "MarlinCED.h"

//---- LCIO ---
#include "IMPL/LCCollectionVec.h"
#include "IMPL/TrackImpl.h"
#include "IMPL/TrackerHitImpl.h"
#include "EVENT/SimTrackerHit.h"
#include "IMPL/LCFlagImpl.h"
#include "UTIL/Operators.h"
#include "UTIL/LCTOOLS.h"


//---- GEAR ----
#include "marlin/Global.h"
#include "gear/GEAR.h"
#include "gear/TPCParameters.h"
#include "gear/PadRowLayout2D.h"
#include "gear/BField.h"

#include "LCIterator.h"

#include "KalTest.h"

using namespace lcio ;
using namespace marlin ;

namespace OuterRimHelper{

  typedef GenericCluster<TrackerHit> HitCluster ;
  typedef GenericClusterVec<TrackerHit> HitClusterVec ;
  typedef GenericHit<TrackerHit>     Hit ;
  typedef GenericHitVec<TrackerHit>  HitVec ;


  //--- delete helper
  template<class P>  void delete_ptr(P* p) { delete p;}

  /** helper class that maps array to gear::Vector3D */
  struct VecFromArray{
    gear::Vector3D _v ;
    VecFromArray( const double* v) : _v( v[0] , v[1] , v[2] ) {}
    VecFromArray( const float* v) : _v( v[0] , v[1] , v[2] ) {}
    const gear::Vector3D& v(){ return _v ; }
  } ;

  /** helper function to restrict the range of the azimuthal angle to ]-pi,pi]*/
  inline double toBaseRange( double phi){
    while( phi <= -M_PI ){  phi += 2. * M_PI ; }
    while( phi >   M_PI ){  phi -= 2. * M_PI ; }
    return phi ;
  }

  /** helper class to compute the chisquared of two points in rho and z coordinate */
  class Chi2_RPhi_Z{
    double _sigsr, _sigsz ;
  public :
    Chi2_RPhi_Z(double sigr, double sigz) : _sigsr( sigr * sigr ) , _sigsz( sigz * sigz ){}
    double operator()( const gear::Vector3D& v0, const gear::Vector3D& v1) {

      //    return (v0 - v1 ).r() ;

      //double dRPhi = v0.rho() * v0.phi() - v1.rho() * v1.phi() ;

      double dPhi = std::abs(  v0.phi() - v1.phi() )  ;
      if( dPhi > M_PI )
	dPhi = 2.* M_PI - dPhi ;

      double dRPhi =  dPhi *  v0.rho() ; 

      double dZ = v0.z() - v1.z() ;

      return  dRPhi * dRPhi / _sigsr + dZ * dZ / _sigsz  ;
    }
  };

  /** helper class to compute the chisquared of two points in rho and z coordinate */
  struct Chi2_RPhi_Z_Hit{
    double operator()( const TrackerHit* h, const gear::Vector3D& v1) {


      gear::Vector3D v0( h->getPosition()[0] ,  h->getPosition()[1] ,  h->getPosition()[2] ) ;

      double sigsr =  sqrt( h->getCovMatrix()[0] + h->getCovMatrix()[2] ) ;
      double sigsz =  h->getCovMatrix()[5] ;
      // double sigsr =  0.01 ; 
      // double sigsz =  0.1 ;
    

      double dPhi = std::abs(  v0.phi() - v1.phi() )  ;
      if( dPhi > M_PI )
	dPhi = 2.* M_PI - dPhi ;

      double dRPhi =  dPhi *  v0.rho() ; 

      double dZ = v0.z() - v1.z() ;

      return  dRPhi * dRPhi / sigsr + dZ * dZ / sigsz  ;
    }
  };
  //------------------------------------------------------

  struct HitInfoStruct{
    HitInfoStruct() :layerID(-1), usedInTrack(false) , chi2Residual(-1.) , deltaChi2(-1.), pNNHit(0) , nNNHit(0), pDist(-1.), nDist(-1.)  {} //, hit(0) {}
    int layerID ;
    bool usedInTrack ;
    double chi2Residual ;  
    double deltaChi2 ;
    Hit* pNNHit ;
    Hit* nNNHit ;
    double pDist ;
    double nDist ;
    //    TVTrackHit* hit ;
  } ;
  struct HitInfo : LCOwnedExtension<HitInfo, HitInfoStruct> {} ;

  // create an owned extension pointer to TVTrackHit - will be deleted when lcio::TrackerHit is deleted
  struct KTHit :  LCOwnedExtension<KTHit,TVTrackHit> {} ;

  //------------------------------------------------------
  // helper class to assign additional parameters to TrackerHit clusters
  struct ClusterInfoStruct{
    ClusterInfoStruct() : track(0){ }
    KalTrack* track ;
    gear::Vector3D nextXPoint ;
    int nextLayer ;
  } ;
  struct ClusterInfo : LCOwnedExtension<ClusterInfo, ClusterInfoStruct> {} ;


  struct TrackInfoStruct{
    TrackInfoStruct() {}
    gear::Vector3D nextXPoint ;
    int nextLayer ;
  } ;
  struct TrackInfo : LCOwnedExtension<TrackInfo, TrackInfoStruct> {} ;


  //------------------------------------------------------
  // function to extract position for Kaltest:
  TVector3 hitPosition( Hit* h)  { 
    return TVector3( h->first->getPosition()[0],   
		     h->first->getPosition()[1],
		     h->first->getPosition()[2]  ) ; 
  }   

  //-----------------------------------------------------------------------------------
  // function to extract layerID from generic Hit:
  int hitLayerID( const Hit* h, int offset=0) { return  h->first->ext<HitInfo>()->layerID + offset  ; } 

  // functor for layer ID
  class HitLayerID{
    int _off ;
    HitLayerID(){}
  public:
    HitLayerID( int off) : _off(off) {}
    int operator()(const Hit* h){ return hitLayerID( h, _off) ; } 
  } ;


  // same for lcio::TrackerHit
  int lcioLayerID( TrackerHit* h, int offset=0) { return  h->ext<HitInfo>()->layerID + offset  ; } 

  // functor for layer ID
  class LCIOLayerID{
    int _off ;
    LCIOLayerID(){}
  public:
    LCIOLayerID( int off) : _off(off) {}
    int operator()(TrackerHit* h){ return lcioLayerID( h, _off) ; } 
  } ;
  
  //-------------------


  struct LCIOTrackerHit{ EVENT::TrackerHit* operator()( Hit* h) { return h->first ; }   } ;

  //---------------------------------------------------
  // helper for sorting cluster wrt layerID
  template <bool SortDirection>
  struct LayerSort{
    bool operator()( const Hit* l, const Hit* r) {
      return hitLayerID( l ) < hitLayerID( r ) ; 
    }
  } ;
  template<>
  struct LayerSort<KalTest::OrderIncoming>{
    bool operator()( const Hit* l, const Hit* r) {
      return hitLayerID( r ) < hitLayerID( l ) ; 
    }
  } ;

  //--- order clusters wrt |omega| - smallest (highest pt) first 
  struct OmegaSort{
    bool operator()( HitCluster* l, HitCluster* r) {
      return std::abs( l->ext<ClusterInfo>()->track->getOmega() )   <  std::abs( r->ext<ClusterInfo>()->track->getOmega()  ) ; 
    }
  };
  
  
  //------ ordering of KalTracks 
  struct KalTrackLengthSort {
    bool operator()( const KalTrack* t0, const KalTrack* t1) {
      return ( t0->getNHits() >= t1->getNHits() );
    }
  };



  //------------------------------
  //helpers for z ordering of hits
  struct TrackerHitCast{
    TrackerHit* operator()(LCObject* o) { return (TrackerHit*) o ; }
  };

  struct ZSort {
    bool operator()( const TrackerHit* l, const TrackerHit* r) {
      return ( l->getPosition()[2] < r->getPosition()[2] );
    }
  };



  void printZ(TrackerHit* h) { 
    std::cout << h->getPosition()[2] << ", " ;
    if(!( h->id() % 30 )) std::cout << std::endl ;
  }



  //-------------------------------
  template <class T>
  void delete_elements(T* t) { delete t ; }

  //-------------------------------

  //-------------------------------------------------------------------------
  template <bool HitOrder, bool FitOrder, bool PropagateIP=false>

  struct KalTestFitter{

    KalTest* _kt ; 
  
    KalTestFitter(KalTest* k) : _kt( k ) {}
  
    KalTrack* operator() (HitCluster* clu) {  
    
      static HitLayerID tpcLayerID( _kt->indexOfFirstLayer( KalTest::DetID::TPC )  )  ;
    
      clu->sort( LayerSort<HitOrder>() ) ;
    
    
      // need to reverse the order for incomming track segments (curlers)
      // assume particle comes from IP
      Hit* hf = clu->front() ;
      Hit* hb = clu->back() ;

      bool reverse_order = ( ( HitOrder ==  KalTest::OrderOutgoing ) ?    
			     ( std::abs( hf->first->getPosition()[2] ) > std::abs( hb->first->getPosition()[2]) + 3. )   :   
			     ( std::abs( hf->first->getPosition()[2] ) < std::abs( hb->first->getPosition()[2]) + 3. )   ) ;
    
      // reverse_order = false ;


      KalTrack* trk = _kt->createKalTrack() ;

      // store mutual pointers between tracks and 'clusters'
      trk->setCluster<HitCluster>( clu ) ;
      //      clu->ext<ClusterInfo>() = new ClusterInfoStruct ;
      clu->ext<ClusterInfo>()->track = trk ;
    

      if( PropagateIP  && HitOrder == KalTest::OrderOutgoing ) {
      
	trk->addIPHit() ;
      }  
    
      // // ----- debug ----
      // std::set<int> layers ;
      // for( HitCluster::iterator it=clu->begin() ; it != clu->end() ; ++it){
      //   if( layers.find( tpcLayerID( *it ) ) != layers.end()  )
      // 	std::cout << " +++++++++++++++++++ duplicate layerID in addHits : " <<  tpcLayerID( *it ) << std::endl ;
      //   layers.insert( tpcLayerID( *it ) ) ;
      // }
      // // ---- end debug ----------
    
      // if( reverse_order )
      // 	trk->addHits( clu->rbegin() , clu->rend() , hitPosition, tpcLayerID , LCIOTrackerHit() ) ; 
      // else
      // 	trk->addHits( clu->begin() , clu->end() , hitPosition, tpcLayerID , LCIOTrackerHit() ) ; 

      if( reverse_order ){
	
	for( HitCluster::reverse_iterator it=clu->rbegin() ; it != clu->rend() ; ++it){   

	  TVTrackHit* h = (*it)->first->ext<KTHit>() ;
	  if( h != 0 )
	    trk->addHit( h ) ; 
	  else
	    streamlog_out( DEBUG ) <<  "   hit not added ;-( " << std::endl ;
	}
	
      } else {
	
	for( HitCluster::iterator it=clu->begin() ; it != clu->end() ; ++it){   

	  TVTrackHit* h = (*it)->first->ext<KTHit>() ;
	  if( h != 0 )
	    trk->addHit( h ) ; 
	  else
	    streamlog_out( DEBUG ) <<  "   hit not added ;-( " << std::endl ;
	}
      }


      if( PropagateIP  && HitOrder == KalTest::OrderIncoming ) {
      
	trk->addIPHit() ;
      }  

      if( reverse_order )
	trk->fitTrack( ! FitOrder  ) ;
      else
 	trk->fitTrack(   FitOrder  ) ;
     
      return trk;
    }
  };


  struct KalTrack2LCIO{
    TrackImpl* operator() (KalTrack* trk) {  
      TrackImpl* lTrk = new TrackImpl ;
      trk->toLCIOTrack( lTrk  ) ;

      if( streamlog_level( DEBUG4 ) ){
	lTrk->ext<TrackInfo>() = new TrackInfoStruct ;       

	streamlog_out( DEBUG4 ) <<  "   ---- KalTrack2LCIO : "  <<  trk->getCluster<HitCluster>()			     
			       <<  " next xing point at layer: "   <<  trk->getCluster<HitCluster>()->ext<ClusterInfo>()->nextLayer
				<< " : " <<  trk->getCluster<HitCluster>()->ext<ClusterInfo>()->nextXPoint ;

	lTrk->ext<TrackInfo>()->nextLayer   =   trk->getCluster<HitCluster>()->ext<ClusterInfo>()->nextLayer ;
	lTrk->ext<TrackInfo>()->nextXPoint  =  trk->getCluster<HitCluster>()->ext<ClusterInfo>()->nextXPoint ;
      }
      return lTrk ;
    }
  };

  //-------------------------------------------------------------------------
  template <class T>
  class RCut {
  public:
    RCut( double rcut ) : _rcut( rcut ) {}  
  
    // bool operator() (T* hit) {  // DEBUG ....
    //   return   std::abs( hit->getPosition()[2] ) > 2000. ;
    bool operator() (T* hit) {  
      return  ( (std::sqrt( hit->getPosition()[0]*hit->getPosition()[0] +
			    hit->getPosition()[1]*hit->getPosition()[1] )   > _rcut )   ||
		( std::abs( hit->getPosition()[2] ) > (500. + _rcut ) )
		); 
    }
  protected:
    RCut() {} ;
    double _rcut ;
  } ;

  template <class T>
  class RCutInverse {
  public:
    RCutInverse( double rcut ) : _rcut( rcut ) {}  
  
    bool operator() (T* hit) {  
      return (  ( std::sqrt( hit->getPosition()[0]*hit->getPosition()[0] +
			     hit->getPosition()[1]*hit->getPosition()[1] )   <= _rcut )   &&
		(  std::abs( hit->getPosition()[2] ) <= (500. + _rcut ) )
		) ;

    }
  protected:
    RCutInverse() {} ;
    double _rcut ;
  } ;


  //---------------------------------------------------------------------------------

  /** Predicate class for identifying clusters with duplicate pad rows - returns true
   * if the fraction of duplicate hits is larger than 'fraction'.
   */
  struct DuplicatePadRows{

    unsigned _N ;
    float _f ; 
    DuplicatePadRows(unsigned nPadRows, float fraction) : _N( nPadRows), _f( fraction )  {}

    bool operator()(const HitCluster* cl) const {
 
      // check for duplicate layer numbers
      std::vector<int> hLayer( _N )  ; 
      typedef HitCluster::const_iterator IT ;

      unsigned nHit = 0 ;
      for(IT it=cl->begin() ; it != cl->end() ; ++it ) {
	TrackerHit* th = (*it)->first ;
	++ hLayer[ th->ext<HitInfo>()->layerID ]   ;
	++ nHit ;
      } 
      unsigned nDuplicate = 0 ;
      for(unsigned i=0 ; i < _N ; ++i ) {
	if( hLayer[i] > 1 ) 
	  nDuplicate += hLayer[i] ;
      }
      return double(nDuplicate)/nHit > _f ;
    }
  };

  //---------------------------------------------------------------------------------

  /** Predicate class for 'distance' of NN clustering.
   */
  //template <class HitClass, typename PosType > 
  class HitDistance{
    typedef TrackerHit HitClass ;
    typedef double PosType ;
  public:

    /** Required typedef for cluster algorithm 
     */
    typedef HitClass hit_type ;

    /** C'tor takes merge distance */
    HitDistance(float dCut) : _dCutSquared( dCut*dCut ) , _dCut(dCut)  {} 


    /** Merge condition: true if distance  is less than dCut given in the C'tor.*/ 
    inline bool mergeHits( GenericHit<HitClass>* h0, GenericHit<HitClass>* h1){
    
      if( std::abs( h0->Index0 - h1->Index0 ) > 1 ) return false ;
    
      //     int l0 =  h0->first->ext<HitInfo>()->layerID ;
      //     int l1 =  h1->first->ext<HitInfo>()->layerID ;

      //     //------- don't merge hits from same layer !
      //     if( l0 == l1 )
      //       return false ;

      if( h0->first->ext<HitInfo>()->layerID == h1->first->ext<HitInfo>()->layerID )
	return false ;

      const PosType* pos0 =  h0->first->getPosition() ;
      const PosType* pos1 =  h1->first->getPosition() ;
    
      return 
	( pos0[0] - pos1[0] ) * ( pos0[0] - pos1[0] ) +
	( pos0[1] - pos1[1] ) * ( pos0[1] - pos1[1] ) +
	( pos0[2] - pos1[2] ) * ( pos0[2] - pos1[2] )   
	< _dCutSquared ;
    }
  
  protected:
    HitDistance() ;
    float _dCutSquared ;
    float _dCut ;
  } ;
  //---------------------------------------------------------------------------------
  
  /** Predicate class for 'distance' of NN clustering based on previous nearest neighbor search.
   */
  class NearestHitDistance{
    typedef TrackerHit HitClass ;
    typedef double PosType ;
  public:
    
    /** Required typedef for cluster algorithm 
     */
    typedef HitClass hit_type ;
    
    /** C'tor takes merge distance */
    NearestHitDistance(float dCut) : _dCutSquared( dCut*dCut ) , _dCut(dCut)  {} 
    
    /** Merge condition: true if distance  is less than dCut given in the C'tor.*/ 
    inline bool mergeHits( GenericHit<HitClass>* h0, GenericHit<HitClass>* h1){
      
      if( std::abs( h0->Index0 - h1->Index0 ) > 1 ) return false ;
      
      int l0 =  h0->first->ext<HitInfo>()->layerID ;
      int l1 =  h1->first->ext<HitInfo>()->layerID ;
      
      if(  std::abs( l0 - l1 ) != 1 ) 
	return false ;
      

      bool ret = ( l0 < l1 ? 
		   h0->first->ext<HitInfo>()->nNNHit ==  h1  &&  h1->first->ext<HitInfo>()->pNNHit ==  h0  :
		   h0->first->ext<HitInfo>()->pNNHit ==  h1  &&  h1->first->ext<HitInfo>()->nNNHit ==  h0  ) ;

      streamlog_out( DEBUG4 ) << "--- NearestHitDistance::mergeHits: "  << std::endl 
			      << " h0: " << h0 
			      << " h0.pNNHit " << h0->first->ext<HitInfo>()->pNNHit
			      << " h0.nNNHit " << h0->first->ext<HitInfo>()->nNNHit << std::endl 
			      << " h1: " << h1 
			      << " h1.pNNHit " << h1->first->ext<HitInfo>()->pNNHit
			      << " h1.nNNHit " << h1->first->ext<HitInfo>()->nNNHit 
			      << "  ---> " <<  ( h0->first->ext<HitInfo>()->nNNHit ==  h1  &&  h1->first->ext<HitInfo>()->pNNHit ==  h0 )
			      << "  , " <<     ( h0->first->ext<HitInfo>()->pNNHit ==  h1  &&  h1->first->ext<HitInfo>()->nNNHit ==  h0 )
			      << std::endl ;

      return ret ;
    }
    
  protected:
    NearestHitDistance() ;
    float _dCutSquared ;
    float _dCut ;
  } ;
  
  //-------------------------------------------------------------------------------------------------------
  class HitDistance_2{
    typedef TrackerHit HitClass ;
    typedef double PosType ;
  public:

    /** Required typedef for cluster algorithm 
     */
    typedef HitClass hit_type ;

    /** C'tor takes merge distance */
    HitDistance_2(float dCut) : _dCutSquared( dCut*dCut ) , _dCut(dCut)  {} 


    /** Merge condition: true if distance  is less than dCut given in the C'tor.*/ 
    inline bool mergeHits( GenericHit<HitClass>* h0, GenericHit<HitClass>* h1){
    
      //if( std::abs( h0->Index0 - h1->Index0 ) > 1 ) return false ;

      int l0 =  h0->first->ext<HitInfo>()->layerID ;
      int l1 =  h1->first->ext<HitInfo>()->layerID ;


      //------- don't merge hits from same layer !
      if( l0 == l1 )
	return false ;


      const PosType* pos0 =  h0->first->getPosition() ;
      const PosType* pos1 =  h1->first->getPosition() ;
    
      return inRange<-2,2>(  l0 - l1 )  &&  
	( pos0[0] - pos1[0] ) * ( pos0[0] - pos1[0] ) +
	( pos0[1] - pos1[1] ) * ( pos0[1] - pos1[1] ) +
	( pos0[2] - pos1[2] ) * ( pos0[2] - pos1[2] )   
	< _dCutSquared ;
    }
  
  protected:
    HitDistance_2() ;
    float _dCutSquared ;
    float _dCut ;
  } ;


  template <class T>
  struct LCIOTrack{
  
    lcio::Track* operator() (GenericCluster<T>* c) {  
    
      TrackImpl* trk = new TrackImpl ;
    
      double e = 0.0 ;
      int nHit = 0 ;
      for( typename GenericCluster<T>::iterator hi = c->begin(); hi != c->end() ; hi++) {
      
	trk->addHit(  (*hi)->first ) ;
	e += (*hi)->first->getEDep() ;
	nHit++ ;
      }

   
      trk->setdEdx( e/nHit ) ;
      trk->subdetectorHitNumbers().push_back( 1 ) ;  // workaround for bug in lcio::operator<<( Tracks ) - used for picking ....
 
      // FIXME - these are no meaningfull tracks - just a test for clustering tracker hits
      return trk ;
    }

  } ;

  // helper for creating lcio header for short printout
  template <class T> 
  const std::string & myheader(){return header(*(T*)(0)); }


  //------------------------------------------------------------------------------------
  void printTrackShort(const LCObject* o){
  
    Track* trk = const_cast<Track*> ( dynamic_cast<const Track*> (o) ) ; 
  
    if( o == 0 ) {
    
      streamlog_out( ERROR ) << "  printTrackShort : dynamic_cast<Track*> failed for LCObject : " << o << std::endl ;
      return  ;
    }
  
    streamlog_out( MESSAGE ) << myheader<Track>()  
			     << lcshort( trk ) << std::endl  ;
  
  
    double r0 = 1. / trk->getOmega() ;
    double d0 = trk->getD0() ;
    double p0 = trk->getPhi() ;
  
    double x0 = ( r0 - d0 ) * sin( p0 ) ;
    double y0 = ( d0 - r0 ) * cos( p0 ) ;
  
    streamlog_out( MESSAGE ) << " circle: r = " << r0 << ", xc = " << x0 << " , yc = " << y0 << std::endl ;
    
    if( streamlog_level( DEBUG ) ){
      if( trk->ext<TrackInfo>() )
	streamlog_out( MESSAGE ) << "  next xing point at layer  "  << trk->ext<TrackInfo>()->nextLayer << " : " <<  trk->ext<TrackInfo>()->nextXPoint ;
    }

  }

  //------------------------------------------------------------------------------------

  void printTrackerHit(const LCObject* o){
  
    TrackerHit* trk = const_cast<TrackerHit*> ( dynamic_cast<const TrackerHit*> (o) ) ; 
  
    if( o == 0 ) {
    
      streamlog_out( ERROR ) << "  printTrackerHit : dynamic_cast<TrackerHit*> failed for LCObject : " << o << std::endl ;
      return  ;
    }
  
    streamlog_out( MESSAGE ) << *trk << std::endl 
			     << " err: rPhi" <<  sqrt( trk->getCovMatrix()[0] + trk->getCovMatrix()[2] ) 
			     << " z :  " <<   trk->getCovMatrix()[5] << std::endl 
			     << " chi2 residual to best matching track : " << trk->ext<HitInfo>()->chi2Residual
			     << " delta chi2 to best matching track : " << trk->ext<HitInfo>()->deltaChi2
			     << std::endl ;

    
  }
  //------------------------------------------------------------------------------------


  /** Predicate class for track merging with NN clustering.
   */
  class TrackStateDistance{
    typedef Track HitClass ;

  public:
  
    /** Required typedef for cluster algorithm 
     */
    typedef HitClass hit_type ;
  
    /** C'tor takes merge distance */
    TrackStateDistance(float dCut) : _dCutSquared( dCut*dCut ) , _dCut(dCut)  {} 


    /** Merge condition: true if distance  is less than dCut given in the C'tor.*/ 
    inline bool mergeHits( GenericHit<HitClass>* h0, GenericHit<HitClass>* h1){
    
      if( std::abs( h0->Index0 - h1->Index0 ) > 1 ) return false ;
    
      Track* trk0 = h0->first ;
      Track* trk1 = h1->first ;

      //------- dont' merge complete tracks: ------------------------
      unsigned nh0 = trk0->getTrackerHits().size() ;
      unsigned nh1 = trk1->getTrackerHits().size() ;
      if( nh0 > 220 ) return false ;
      if( nh1 > 220 ) return false ;
      //------------------------------------------------------


      KalTrack* ktrk0 = h0->first->ext<KalTrackLink>() ; 
      KalTrack* ktrk1 = h1->first->ext<KalTrackLink>() ; 

      //---- sanity check on radii-------------------
      double r0 = 1. / trk0->getOmega() ;
      double r1 = 1. / trk1->getOmega() ;

      if( r0 < 300. || r1 < 300. )
	return false ;

      if( std::abs( r0 - r1 ) / std::abs( r0 + r1 )  > 0.02 )  // relative difference larger than 1%
	return false ;
      //---------------------------------------------


      double chi2  = KalTrack::chi2( *ktrk0 , *ktrk1 ) ; 	   

      return chi2  < _dCut ;
      
    }
  
  protected:
    TrackStateDistance() ;
    float _dCutSquared ;
    float _dCut ;
  } ;

  //------------------------------------------------------------------------------------


  /** helper class for merging track segments, based on circle (and tan lambda) */
  class TrackCircleDistance{
    typedef Track HitClass ;

  public:
  
    /** Required typedef for cluster algorithm
     */
    typedef HitClass hit_type ;
  
  
    /** C'tor takes merge distance */
    TrackCircleDistance(float dCut) : _dCutSquared( dCut*dCut ) , _dCut(dCut){}
  
    /** Merge condition: ... */
    inline bool mergeHits( GenericHit<HitClass>* h0, GenericHit<HitClass>* h1){
    
      static const double DRMAX = 0.1 ; // make parameter
      static const double DTANLMAX = 0.2 ; //   " 

      if( std::abs( h0->Index0 - h1->Index0 ) > 1 ) return false ;
    
      Track* trk0 = h0->first ;
      Track* trk1 = h1->first ;

      // //------- dont' merge complete tracks: ------------------------
      // unsigned nh0 = trk0->getTrackerHits().size() ;
      // unsigned nh1 = trk1->getTrackerHits().size() ;
      // if( nh0 > 220 ) return false ;
      // if( nh1 > 220 ) return false ;
      // //------------------------------------------------------

      // KalTrack* ktrk0 = h0->first->ext<KalTrackLink>() ; 
      // KalTrack* ktrk1 = h1->first->ext<KalTrackLink>() ; 

      double tl0 = trk0->getTanLambda() ;
      double tl1 = trk1->getTanLambda() ;

      double dtl = 2. * ( tl0 - tl1 ) / ( tl0 + tl1 ) ;
      dtl *= dtl ;
      if(  dtl >  DTANLMAX * DTANLMAX ) 
	return false ;

      double r0 = 1. / trk0->getOmega()  ;
      double r1 = 1. / trk1->getOmega()  ;


      double d0 = trk0->getD0() ;
      double d1 = trk1->getD0() ;

      double p0 = trk0->getPhi() ;
      double p1 = trk1->getPhi() ;

      double x0 = ( r0 - d0 ) * sin( p0 ) ;
      double x1 = ( r1 - d1 ) * sin( p1 ) ;

      double y0 = ( d0 - r0 ) * cos( p0 ) ;
      double y1 = ( d1 - r1 ) * cos( p1 ) ;
    
      double dr = 2. * std::abs( ( r0 -r1 ) / (r0 + r1 ) ) ;

      double distMS = sqrt ( ( x0 - x1 ) * ( x0 - x1 ) + ( y0 - y1 ) * ( y0 - y1 )  ) ;
    
    
      return ( dr < DRMAX && distMS < _dCut*r0 ) ;

    }
  
  protected:
    float _dCutSquared ;
    float _dCut ;
  } ; 

  //-----------------------------------------------------------------

  /** Find the nearest hits in the previous and next layers - (at most maxStep layers appart - NOT YET).
   */
  void findNearestHits( HitCluster& clu, KalTest* kt, int maxStep=1){

    typedef std::list<Hit*> HitList ;
    typedef std::vector< HitList > HitListVector ;
    HitListVector hitsInLayer( kt->maxLayerIndex() ) ;
    HitLayerID  tpcLayerID( kt->indexOfFirstLayer( KalTest::DetID::TPC )  ) ;
    
    for(HitCluster::iterator it= clu.begin() ; it != clu.end() ; ++it ){
      Hit* hit = *it ; 
      hitsInLayer[ tpcLayerID( hit ) ].push_back( hit )  ;
    }
    //-----------------------------

    for(HitCluster::iterator it= clu.begin() ; it != clu.end() ; ++it ){
      
      Hit* hit0 = *it ; 
      gear::Vector3D pos0( hit0->first->getPosition()[0] , hit0->first->getPosition()[1] , hit0->first->getPosition()[2] ) ;
      int layer = tpcLayerID( hit0 ) ;
      
      int l=0 ;
      if( (l=layer+1) <  kt->maxLayerIndex() ) {
	
	HitList& hL = hitsInLayer[ l ] ;
	
	double minDist2 = 1.e99 ;
	Hit*   bestHit = 0 ;

	for( HitList::iterator iH = hL.begin() ; iH != hL.end() ; ++iH ) {
	  
	  Hit* hit1 = *iH ; 
	  gear::Vector3D pos1( hit1->first->getPosition()[0] , hit1->first->getPosition()[1] , hit1->first->getPosition()[2] ) ;
	  
	  double dist2 = ( pos0 - pos1 ).r2() ;

	  if( dist2 < minDist2 ){
	    minDist2 = dist2 ;
	    bestHit = hit1 ;
	  }
	}
	if( bestHit ){

	  hit0->first->ext<HitInfo>()->nNNHit =  bestHit ;
	  hit0->first->ext<HitInfo>()->nDist  =  minDist2 ;
	}
      }

      if( (l=layer-1) > 0  ) {
	
	HitList& hL = hitsInLayer[ l ] ;
	
	double minDist2 = 1.e99 ;
	Hit* bestHit = 0 ;
	
	for( HitList::iterator iH = hL.begin() ; iH != hL.end() ; ++iH ) {
	  
	  Hit* hit1 = *iH ; 
	  gear::Vector3D pos1( hit1->first->getPosition()[0] , hit1->first->getPosition()[1] , hit1->first->getPosition()[2] ) ;
	  
	  double dist2 = ( pos0 - pos1 ).r2() ;
	  
	  if( dist2 < minDist2 ){
	    minDist2 = dist2 ;
	    bestHit = hit1 ;
	  }
	}
	if( bestHit ){

	  hit0->first->ext<HitInfo>()->pNNHit =  bestHit ;
	  hit0->first->ext<HitInfo>()->pDist  =  minDist2 ;
	}
      }

    }
  }




}  // namespace OuterRimhelper

using namespace OuterRimHelper ;


//========================================================================================================

OuterRimSearch aOuterRimSearch ;


//*********************************************************************************************************
OuterRimSearch::OuterRimSearch() : Processor("OuterRimSearch") {
//*********************************************************************************************************
  

  // modify processor description
  _description = "OuterRimSearch : simple nearest neighbour clustering" ;
  
  
  StringVec colDefault ;
  colDefault.push_back("AllTPCTrackerHits" ) ;

  registerInputCollections( LCIO::TRACKERHIT,
			    "HitCollections" , 
			    "Name of the input collections"  ,
			    _colNames ,
			    colDefault ) ;
  
  registerOutputCollection( LCIO::TRACK,
			    "OutputCollection" , 
			    "Name of the output collections"  ,
			    _outColName ,
			    std::string("CluTracks" ) ) ;
  
  
  registerProcessorParameter( "DistanceCut" , 
			      "Cut for distance between hits in mm"  ,
			      _distCut ,
			      (float) 40.0 ) ;
  
  registerProcessorParameter( "MinimumClusterSize" , 
			      "minimum number of hits per cluster"  ,
			      _minCluSize ,
			      (int) 3) ;
  

  registerProcessorParameter( "DuplicatePadRowFraction" , 
			      "allowed fraction of hits in same pad row per track"  ,
			      _duplicatePadRowFraction,
			      (float) 0.01 ) ;

  registerProcessorParameter( "RCut" , 
 			      "Cut for r_min in mm"  ,
 			      _rCut ,
 			      (float) 0.0 ) ;
  
}


void OuterRimSearch::init() { 

  // usually a good idea to
  printParameters() ;


  _kalTest = new KalTest( *marlin::Global::GEAR ) ;
  _kalTest->setOption( KalTest::CFG::ownsHits, false ) ;
  _kalTest->init() ;

  _nRun = 0 ;
  _nEvt = 0 ;
}

void OuterRimSearch::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 

void OuterRimSearch::processEvent( LCEvent * evt ) { 

  clock_t start =  clock() ; 

  GenericHitVec<TrackerHit> h ;
  GenericHitVec<TrackerHit> hSmallR ; 
  
  GenericClusterVec<TrackerHit> cluList ;
  
  RCut<TrackerHit> rCut( _rCut ) ;
  RCutInverse<TrackerHit> rCutInverse( _rCut ) ;
  
  ZIndex<TrackerHit,200> zIndex( -2750. , 2750. ) ; 
  
  //  NNDistance< TrackerHit, double> dist( _distCut )  ;

  HitDistance dist0( _distCut ) ;
  HitDistance dist( 20. ) ;
  
  LCIOTrack<TrackerHit> converter ;
  

  // -------------------- assign layer number to Trackerhits ---------------------------------
  
  const gear::TPCParameters& gearTPC = Global::GEAR->getTPCParameters() ;
  const gear::PadRowLayout2D& padLayout = gearTPC.getPadLayout() ;
  unsigned nPadRows = padLayout.getNRows() ;
  
  // create a vector of generic hits from the collection applying a cut on r_min
  for( StringVec::iterator it = _colNames.begin() ; it !=  _colNames.end() ; it++ ){  
    
    LCCollectionVec* col =  dynamic_cast<LCCollectionVec*> (evt->getCollection( *it )  ); 
    
    int nHit = col->getNumberOfElements() ;
    for(int i=0 ; i < nHit ; ++i ) {
      
      TrackerHitImpl* th = (TrackerHitImpl*) col->getElementAt(i) ;
      
      gear::Vector3D v( th->getPosition()[0],th->getPosition()[1], 0 ) ; 
      
      int padIndex = padLayout.getNearestPad( v.rho() , v.phi() ) ;
      
      th->ext<HitInfo>() = new HitInfoStruct ;
      
      th->ext<HitInfo>()->layerID = padLayout.getRowNumber( padIndex ) ;

      // add a TVTrackHit:
      th->ext<KTHit>() = _kalTest->createHit( th , lcioLayerID( th )  , KalTest::DetID::TPC ) ;
      
      
    } //-------------------- end assign layernumber -----------------------------------------
    
    
    std::list< TrackerHit*> hitList ;
    TrackerHitCast cast ;
    ZSort zsort ;
    std::transform(  col->begin(), col->end(), std::back_inserter( hitList ), cast ) ;

    hitList.sort( zsort ) ;
    //    std::for_each( hitList.begin() , hitList.end() , printZ ) ;

    addToGenericHitVec( h, hitList.begin() , hitList.end() , rCut ,  zIndex ) ;
    
    // create a vector with the hits at smaller R
    addToGenericHitVec( hSmallR, hitList.begin() , hitList.end() , rCutInverse ,  zIndex ) ;
  } 
  
  // cluster the sorted hits  ( if |diff(z_index)|>1 the loop is stopped)
  cluster_sorted( h.begin() , h.end() , std::back_inserter( cluList )  , &dist0 , _minCluSize ) ;
  //cluster( h.begin() , h.end() , std::back_inserter( cluList )  , &dist , _minCluSize ) ;
  
  streamlog_out( DEBUG ) << "   ***** clusters: " << cluList.size() << std::endl ; 
  
  LCCollectionVec* allClu = new LCCollectionVec( LCIO::TRACK ) ;
  std::transform(cluList.begin(), cluList.end(), std::back_inserter( *allClu ) , converter ) ;
  evt->addCollection( allClu , "AllTrackClusters" ) ;
  

  //---------------------------------------------------------------------------------------------
  // find 'odd' clusters that have duplicate hits in pad rows

  HitClusterVec ocs ;

  split_list( cluList, std::back_inserter(ocs),  DuplicatePadRows( nPadRows, _duplicatePadRowFraction  ) ) ;

  LCCollectionVec* oddCol = new LCCollectionVec( LCIO::TRACK ) ;
  std::transform( ocs.begin(), ocs.end(), std::back_inserter( *oddCol ) , converter ) ;
  evt->addCollection( oddCol , "OddClu_1" ) ;

  streamlog_out( DEBUG ) << "   ***** clusters: " << cluList.size() 
			 << "   ****** oddClusters " << ocs.size() 
			 << std::endl ; 


  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // clean up the clusters with duplicate pad rows ....
  //  
  HitClusterVec cclu ; // clean clusters
 
  NearestHitDistance nnDist(0.) ;
 
  for( HitClusterVec::iterator iC = ocs.begin() ; iC != ocs.end() ; ++iC ){

    HitCluster* clu = *iC ;

    findNearestHits( *clu , _kalTest ) ;

    //---- get hits from cluster into a vector
    std::vector< Hit* > v ;
    v.reserve( clu->size() ) ;
    clu->takeHits( std::back_inserter( v )  ) ;
    delete clu ;

    //-- recluster with nearest neighbours

    cluster( v.begin(), v.end() , std::back_inserter( cclu ),  &nnDist  , _minCluSize ) ;

  }
  ocs.clear() ;

  LCCollectionVec* ccol = new LCCollectionVec( LCIO::TRACK ) ;
  std::transform( cclu.begin(), cclu.end(), std::back_inserter( *ccol ) , converter ) ;
  evt->addCollection( ccol , "cleaned_clusters" ) ;


  cluList.merge( cclu ) ;
  
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  LCCollectionVec* cluCol = new LCCollectionVec( LCIO::TRACK ) ;
  std::transform( cluList.begin(), cluList.end(), std::back_inserter( *cluCol ) , converter ) ;
  evt->addCollection( cluCol , "CluTrackSegments" ) ;



  //---------------------------------------------------------------------------------------------
  // create vector with left over hits

  std::vector< Hit* > leftOverHits ;
  leftOverHits.reserve(  h.size() + hSmallR.size() ) ;

  typedef HitVec::const_iterator GHVI ;

  for( GHVI it = h.begin(); it != h.end() ; ++it ){

    if ( (*it)->second == 0 ) leftOverHits.push_back( *it ) ;
  }

  // add all hits that failed the rcut 
  std::copy( hSmallR.begin() , hSmallR.end() , std::back_inserter( leftOverHits )  ) ;
  
  
  
  //*********************************************************
  //   run KalTest on good track segments 
  //*********************************************************

  for( HitClusterVec::iterator icv = cluList.begin() ; icv != cluList.end() ; ++ icv ) {
    HitCluster* clu  = *icv ;
    clu->ext<ClusterInfo>() = new ClusterInfoStruct ;
  }

  streamlog_out( DEBUG ) <<  "************* fitted segments and KalTest tracks : **********************************" 
			 << std::endl ;
  
  
  std::list< KalTrack* > ktracks ;
  
  KalTestFitter<KalTest::OrderOutgoing, KalTest::FitBackward > fitter( _kalTest ) ;
  
  std::transform( cluList.begin(), cluList.end(), std::back_inserter( ktracks ) , fitter ) ;
  

  //  std::for_each( ktracks.begin(), ktracks.end(), std::mem_fun( &KalTrack::findXingPoints ) ) ;
  
  
  LCCollectionVec* trksegs = new LCCollectionVec( LCIO::TRACK ) ;
  std::transform( ktracks.begin(), ktracks.end(), std::back_inserter( *trksegs ) , KalTrack2LCIO() ) ;
  evt->addCollection( trksegs , "KalTrackSegments" ) ;

  
  //  //=========== assign left over hits ... ==================================================================

  //------------- create vector of left over hits per layer
  typedef std::list<Hit*> HitList ;
  typedef std::vector< HitList > HitListVector ;
  HitListVector hitsInLayer( _kalTest->maxLayerIndex() ) ;

  HitLayerID  tpcLayerID( _kalTest->indexOfFirstLayer( KalTest::DetID::TPC )  ) ;
 
  
  for( GHVI ih = leftOverHits.begin() ; ih != leftOverHits.end() ; ++ih ) {
    
    Hit* hit = *ih ;
    //      std::cout << " ++++++  layerId: " << tpcLayerID( hit ) << " max layer index : " <<  _kalTest->maxLayerIndex() << std::endl  ;
    hitsInLayer[ tpcLayerID( hit ) ].push_back( hit )  ;
  }
  //-----------------------------


  OmegaSort omegaSort ;
  cluList.sort( omegaSort ) ;

  Chi2_RPhi_Z_Hit ch2rzh ;
  static const double chi2Cut = 100 ; //15. ; // FIXME: make parameter

  
  if( streamlog_level( DEBUG4 ) ) {
    for( HitClusterVec::iterator icv = cluList.begin() ; icv != cluList.end() ; ++ icv ) {
      HitCluster* clu  = *icv ;
      KalTrack* trk =  clu->ext<ClusterInfo>()->track ;
      gear::Vector3D xv ;
      int  layer ;
      trk->findNextXingPoint(  xv , layer ) ;
      
      clu->ext<ClusterInfo>()->nextXPoint = xv ;
      clu->ext<ClusterInfo>()->nextLayer = layer ;

      streamlog_out( DEBUG4 ) <<  "   ----  FINDNEXTXINGPOINT: "  <<  clu
			      <<  " next xing point at layer: "   <<  clu->ext<ClusterInfo>()->nextLayer
			      << " : " <<  clu->ext<ClusterInfo>()->nextXPoint ;
    }
  }
  
  for( HitClusterVec::iterator icv = cluList.begin() ; icv != cluList.end() ; ++ icv ) {
    
    HitCluster* clu  = *icv ;
    KalTrack* trk =  clu->ext<ClusterInfo>()->track ;
    
    static const int maxStep =  3 ;   //FIXME: make parameter - max step size  (#layers) w/o hit 
    int step = 0 ;
    
    while( step < maxStep + 1 ) {
      
      gear::Vector3D xv ;
      int layer ;
      
      bool hitAdded = false ;

      trk->findNextXingPoint(  xv , layer , step ) ;
      
      streamlog_out( DEBUG4 ) <<  "  -- searching in leftover hits for cluster : " << std::endl 
			     <<  "  omega : " <<  trk->getOmega()  
			     <<  "  Step : " << step 
			     <<  "  at layer: "   << layer         
			     <<  "  next xing point : " <<  xv  ;
      
      if( layer > -1 ) { // found a crossing point 
	
 	HitList& hLL = hitsInLayer.at( layer ) ;
	
	double ch2Min = 1.e99 ;
	Hit* bestHit = 0 ;
	
 	for( HitList::const_iterator ih = hLL.begin() ; ih != hLL.end() ; ++ih ){
	  
 	  Hit* hit = *ih ;
	  
 	  double ch2 = ch2rzh( hit->first , xv )  ;
	  
 	  if( ch2 < ch2Min ){
	    
 	    ch2Min = ch2 ;
 	    bestHit = hit ;
 	  }
 	}
	
 	if( bestHit != 0 ){
	  
 	  VecFromArray hPos(  bestHit->first->getPosition() ) ;
	  
 	  bestHit->first->ext<HitInfo>()->chi2Residual = ch2Min ;
	  
	  if( ch2Min  < chi2Cut ) { 
	  
	    streamlog_out( DEBUG ) <<   " ---- assigning left over hit : " << hPos.v() << " <-> " << xv
				   <<   " dist: " <<  (  hPos.v() - xv ).r()
				   <<   " chi2: " <<  ch2Min 
				   <<   "  hit errors :  rphi=" <<  sqrt( bestHit->first->getCovMatrix()[0] + bestHit->first->getCovMatrix()[2] ) 
				   <<	 "  z= " <<  sqrt( bestHit->first->getCovMatrix()[5] )
				   << std::endl ;
	    
	    
	    
	    //	    if( trk->addAndFilter( bestHit->first->ext<KTHit>() ) ) {
	      
	    double deltaChi2 =  trk->testDeltaChi2(  bestHit->first->ext<KTHit>() ) ;

	    bestHit->first->ext<HitInfo>()->deltaChi2 = deltaChi2 ;

	    if(  deltaChi2 < 25.   &&
		 trk->addAndFilter( bestHit->first->ext<KTHit>() )  ){

	      hitAdded = true ;
	      
	      hLL.remove(  bestHit ) ;
	      clu->addHit( bestHit ) ;
	      
	      streamlog_out( DEBUG ) <<   " ---- track state filtered with new hit ! ------- " << std::endl ;

	    }

	  } // chi2Cut 

	} // bestHit

      } // layer > -1 	

      if( hitAdded ){     step = 1 ;
      } else {  	++step  ;    }
      
    } // while step < maxStep
  }

  //  static const bool use_best_track = false ;

  //  if( use_best_track ) {

  //    streamlog_out( DEBUG ) << "  ------ assign left over hits - best matching track for every hit ..."  << std::endl ;

  //    Chi2_RPhi_Z ch2rz( 0.1 , 1. ) ; // fixme - need proper errors ....
  //    Chi2_RPhi_Z_Hit ch2rzh ;


  //    HitLayerID  tpcLayerID( _kalTest->indexOfFirstLayer( KalTest::DetID::TPC )  ) ;
    
  //    for( GHVI ih = leftOverHits.begin() ; ih != leftOverHits.end() ; ++ih ){
      
  //      Hit* hit = *ih ;
  //      VecFromArray hPos(  hit->first->getPosition() ) ;
      
  //      double ch2Min = 999999999999999. ;
  //      KalTrack* bestTrk = 0 ;
      
  //      for( std::list< KalTrack* >::iterator it = ktracks.begin() ; it != ktracks.end() ; ++it ){
	
  // 	const gear::Vector3D* kPos = (*it)->getXingPointForLayer( tpcLayerID( hit ) ) ;
	
  // 	// double rh  =  hPos.v().rho() ;
  // 	// double rk  =  kPos->rho() ;
  // 	// if( std::abs( rh - rk ) > 0.1 ) {
  // 	// 	streamlog_out( WARNING ) << " --- different radii for hit and crossing point : " <<  tpcLayerID( hit ) << ": " << rh << " - " << rk 
  // 	// 				 <<  *kPos  << std::endl ;
  // 	// } 
	
  // 	if( kPos != 0 ){
	  
  // 	  //	  double ch2 = ch2rz( hPos.v() , *kPos )  ;
  // 	  double ch2 = ch2rzh( hit->first , *kPos )  ;
	  
  // 	  if( ch2 < ch2Min ){
	    
  // 	    ch2Min = ch2 ;
  // 	    bestTrk = *it ;
  // 	  }
	  
  // 	}
	
  // 	// else {
  // 	// 	streamlog_out( MESSAGE ) << " --- no crossing point found for layer : " <<  tpcLayerID( hit ) << ": " << hPos.v() << std::endl ;
  // 	// }
	
  //      }
  //      if( bestTrk ) {
	
  // 	const gear::Vector3D* kPos = bestTrk->getXingPointForLayer( tpcLayerID( hit ) ) ;
	
  // 	// double rh  =  hPos.v().rho() ;
  // 	// double rk  =  kPos->rho() ;
  // 	// if( std::abs( rh - rk ) > 0.1 ) {
  // 	// 	streamlog_out( WARNING ) << "  different radii for hit and crossing point : " << rh << " - " << rk << std::endl ;
  // 	// } 
	
  // 	//      if( std::abs( hPos.v().rho() - kPos->rho() ) < 0.5 &&   std::abs( hPos.v().z() - kPos->z() ) < 5. ) {
	
  // 	if(  (  hPos.v() - *kPos ).r()  < 3. ) {   // check for bad outliers... FIXME: need proper criterion here .....
	  
	  
  // 	  HitCluster* clu = bestTrk->getCluster< HitCluster >() ;
	  
  // 	  streamlog_out( DEBUG ) << " ---- assigning left over hit : " << hPos.v() << " <-> " << *kPos  
  // 				 <<   " dist: " <<  (  hPos.v() - *kPos ).r()  << std::endl ;
	  
  // 	  clu->addHit( hit ) ;
  // 	}	
  // 	else 
  // 	  streamlog_out( DEBUG ) << " ---- NOT assigning left over hit : " << hPos.v() << " <-> " << *kPos << std::endl ;
  //      }
  //      else
  // 	streamlog_out( DEBUG ) << " ---- NO best track found ??? ---- " << std::endl ;
      
  //    }
    

  //    //        ==========================================================================================
  //  } else { // ================== use best matching hit for every track segment =========================
  //    //        ==========================================================================================
    


  //    streamlog_out( DEBUG1 ) << "  ------ assign left over hits - best matching hit for every track ..."  << std::endl ;
    
  //    HitLayerID  tpcLayerID( _kalTest->indexOfFirstLayer( KalTest::DetID::TPC )  ) ;
    

  //    //------------- create vector of left over hits per layer
  //    typedef std::list<Hit*> HitList ;
  //    typedef std::vector< HitList > HitListVector ;
  //    HitListVector hitsInLayer( _kalTest->maxLayerIndex() ) ;
    
    
  //    for( GHVI ih = leftOverHits.begin() ; ih != leftOverHits.end() ; ++ih ) {
      
  //      Hit* hit = *ih ;
  //      //      std::cout << " ++++++  layerId: " << tpcLayerID( hit ) << " max layer index : " <<  _kalTest->maxLayerIndex() << std::endl  ;
  //      hitsInLayer[ tpcLayerID( hit ) ].push_back( hit )  ;
  //    }
  //    //-----------------------------
    
  //    std::map< HitCluster* , KalTrack* > clu2trkMap ;

  //    const bool use_segment_hits = false ; //true ;
    
  //    if( use_segment_hits  ){
      
  //      // store first and last hit of every segment in map with leftover hits in this layer
      
  //      for( GenericClusterVec<TrackerHit>::iterator icv = cluList.begin() ; icv != cluList.end() ; ++ icv ) {
	
  // 	Hit* h0 = (*icv)->front() ;
  // 	Hit* h1 = (*icv)->back() ;
	
  // 	hitsInLayer[ tpcLayerID( h0 ) ].push_back( h0 )  ;
  // 	hitsInLayer[ tpcLayerID( h1 ) ].push_back( h1 )  ;
  //      }
      
  //      // sort the tracks wrt. lenghts (#hits)
  //      ktracks.sort( KalTrackLengthSort() ) ;

  //      // store assoaciation between cluster and track 
  //      for( std::list< KalTrack* >::iterator it = ktracks.begin() ; it != ktracks.end() ; ++it ){
  // 	HitCluster* c = (*it)->getCluster< HitCluster >() ;
  // 	clu2trkMap[ c ] = *it ;
  //      }	   
  //    }
  //    //-------------------------------
    

  //    //    Chi2_RPhi_Z ch2rz( 0.1 , 1. ) ; // fixme - need proper errors 
  //    Chi2_RPhi_Z_Hit  ch2rzh ;
    

  //    for( std::list< KalTrack* >::iterator it = ktracks.begin() ; it != ktracks.end() ; ++it ){
      
  //      KalTrack* theTrack = *it ;
  //      if( theTrack == 0 ) 
  // 	continue ;
      
      
  //      // ----- define chi2 cut    ~15 for 1 GeV pt 
  //      double chi2Cut = 100000. / ( std::log(1.) - std::log( std::abs(theTrack->getOmega()) ) ) ;


  //      streamlog_out( DEBUG3 ) << " ------- searching for leftover hits for track : " << theTrack 
  // 			      << "   chi2 cut : " << chi2Cut  << " -  omega : " << theTrack->getOmega() <<  std::endl ;
      
  //      int xpLayer = 0 ;
      
  //      // const PointList& xptList = theTrack->getXingPoints() ;
  //      // for(PointList::const_iterator itXP = xptList.begin() ; itXP != xptList.end() ; ++itXP , xpLayer++ ) {
  //      // 	const gear::Vector3D* kPos =  *itXP ;
      
  //      PointList& xpVec = theTrack->getXingPoints() ;
  //      for( unsigned ixp=0 ; ixp < xpVec.size() ; ++ixp, xpLayer++  ) {
  //      	const gear::Vector3D* kPos =  xpVec[ixp]  ;
	
  // 	if( kPos == 0 ) {   // we don't have a xing point
  // 	  continue ;
  // 	}
	
  //       	double ch2Min = 10e99 ;
  // 	Hit* bestHit = 0 ;
	
  // 	HitList& hLL = hitsInLayer.at( xpLayer ) ;
	
  // 	for( HitList::const_iterator ih = hLL.begin() ; ih != hLL.end() ; ++ih ){
	  
  // 	  Hit* hit = *ih ;
	  
  // 	  //VecFromArray hPos(  hit->first->getPosition() ) ;
  // 	  //double ch2 = ch2rz( hPos.v() , *kPos )  ;
  // 	  double ch2 = ch2rzh( hit->first , *kPos )  ;

  // 	  if( ch2 < ch2Min ){
	    
  // 	    ch2Min = ch2 ;
  // 	    bestHit = hit ;
  // 	  }
  // 	}
	
  // 	if( bestHit != 0 ) {
	  
  // 	  VecFromArray hPos(  bestHit->first->getPosition() ) ;
	  
  // 	  //	  if( ch2Min  <  6. ) { // Sum( pdf(ch2,ndf==2) )_0^6 ~ 95% )
  // 	  //	  if( ch2Min  <  20. ) { // Sum( pdf(ch2,ndf==2) )_0^20 ~ 99.x% ?? ) // FIXME: need steering parameter and optimize value
	  
	  
  // 	  bestHit->first->ext<HitInfo>()->chi2Residual = ch2Min ;


  // 	  if( ch2Min  < chi2Cut ) { 
	    
  // 	    streamlog_out( DEBUG1 ) <<   " ---- assigning left over hit : " << hPos.v() << " <-> " << *kPos
  // 				    <<   " dist: " <<  (  hPos.v() - *kPos ).r()
  // 				    <<   " chi2: " <<  ch2Min 
  // 				    <<   "  hit errors :  rphi=" <<  sqrt( bestHit->first->getCovMatrix()[0] + bestHit->first->getCovMatrix()[2] ) 
  // 				    <<	 "  z= " <<  sqrt( bestHit->first->getCovMatrix()[5] )
  // 				    << std::endl ;
	    
	    
  // 	    if( bestHit->second != 0 ) { //--------------------------------------------------------------------------------
	      
  // 	      // hit is already part of a track segment 
	      
	      
  // 	      HitCluster* c = bestHit->second  ;
  // 	      KalTrack* trk = clu2trkMap[ c ] ;
	      

  // 	      if( trk == theTrack ) {
  // 		streamlog_out( ERROR ) << " =======================  found best matching hit from track itself: " << *bestHit->first 
  // 				       <<     std::endl  
  // 				       <<  "      track has  " << trk->getNHits()  << " hits " << std::endl ;

  // 		for( unsigned ii=0 ; ii < xpVec.size() ; ++ii) {
  // 		  if( xpVec[ii] ) 
  // 		    streamlog_out( ERROR ) << "  xing pt : "  << ii << " - " << *xpVec[ii]  ;
  // 		}
		
		
  // 		for( HitCluster::iterator its = c->begin(); its != c->end() ; ++its ){
  // 		  Hit* hit = *its ;
  // 		  VecFromArray hPos(  hit->first->getPosition() ) ;
  // 		  streamlog_out( ERROR ) << "  hit  : layer: "  <<   tpcLayerID( hit )   << " - " << hPos.v()  ;
  // 		}
		

  // 	      } else {

		
  // 		streamlog_out( DEBUG3 ) << " +++++++++ found best hit already part of a track segment !!!!!! " 
  // 					<< " trk : " << trk  << " #hits: " << trk->getNHits() 
  // 					<< " cluster " << c << c->size() 
  // 					<< std::endl ;   
		
		
  // 		unsigned goodHits(0), allHits(0) ;
		
  // 		double chi2Max = 10. ; // fixme parameter
		
  // 		for( HitCluster::iterator its = c->begin(); its != c->end() ; ++its ){
		  
  // 		  ++allHits ;
		  
  // 		  Hit* hit = *its ;
  // 		  VecFromArray hPos(  hit->first->getPosition() ) ;
		  
  // 		  const gear::Vector3D* kPos = theTrack->getXingPointForLayer( tpcLayerID( hit ) ) ;
		  
  // 		  if( kPos != 0 ) {
		    
  // 		    //double chi2 = ch2rz( hPos.v() , *kPos )  ;
  // 		    double chi2 = ch2rzh( hit->first , *kPos )  ;

  // 		    streamlog_out( DEBUG3 ) << " +++++++++ chi2 : " << chi2 << hPos.v() 
  // 					    << " +++++++++                  " << *kPos 
  // 					    << " +++++++++  hit id " << std::hex << hit->first->id() << std::dec 
  // 					    << std::endl ;
		    
  // 		    if( chi2 < chi2Max ){
		      
  // 		      ++goodHits ;
  // 		    }
  // 		  }
  // 		}
		
  // 		double goodFraction = double( goodHits ) / double(  allHits ) ;
		
  // 		streamlog_out( DEBUG3 ) << " +++++++++ fraction of matching hits : " << goodFraction 
  // 					<< std::endl ;   
		
		
  // 		// ---------------------  merge the track segements -------------------------
		
  // 		if( goodFraction > 0.5  ) { // fixme: what is reasonable here - make parameter ...
		  
		  
  // 		  for( HitCluster::iterator its = c->begin(); its != c->end() ; ++its ){

  // 		    delete  xpVec[  tpcLayerID( *its ) ] ; // erase crossing points for these hit
  // 		    xpVec[  tpcLayerID( *its ) ]  = 0 ;   
  // 		  }
  // 		  HitCluster* clu = theTrack->getCluster< HitCluster >() ;
		  
  // 		  // merge the cluster into the larger one and delete it - remove the hits from the hitsinlayer vector first
		  
  // 		  Hit* h0 = c->front() ;
  // 		  Hit* h1 = c->back() ;
		  
  // 		  hitsInLayer[ tpcLayerID( h0 ) ].remove( h0 )  ;
  // 		  hitsInLayer[ tpcLayerID( h1 ) ].remove( h1 )  ;
		  
  // 		  clu->mergeClusters( c ) ;
		  
  // 		  cluList.remove( c  ) ;
		  
		  
  // 		  streamlog_out( DEBUG3) << " ************ found matching segment, merged all hits: delete cluster : " << c 
  // 					 << " and track : " << trk << std::endl ;
		  
  // 		  delete c ;
		  
  // 		  ktracks.remove( trk ) ;
		  
  // 		} //-------------------------------------------------------------

  // 	      }

		
  // 	    }  else  {  //--------------------------------------------------------------------------------
	      
  // 	      hLL.remove(  bestHit ) ;
	      
  // 	      HitCluster* clu = theTrack->getCluster< HitCluster >() ;
	      
  // 	      streamlog_out( DEBUG3) << "    ************ found matching hit, add to  cluster : " << clu  << std::endl ;
	      
  // 	      clu->addHit( bestHit ) ;
  // 	    }


  // 	  }
  // 	} 
  //          // else {
  // 	  //	  streamlog_out( DEBUG1 ) << "????????????????????? no best Hit found xing pnt  : chi2  " << *xpVec[ixp]  << " : " << ch2Min << std::endl ;
  // 	  //	}

  //      }
  //    }
    
  //  }

  //   //================================================================================================================ 
  


  //  //std::transform( ktracks.begin(), ktracks.end(), std::back_inserter( *kaltracks ) , KalTrack2LCIO() ) ;

  std::list< KalTrack* > newKTracks ;

  //  //KalTestFitter<KalTest::OrderIncoming, KalTest::FitForward, KalTest::PropagateToIP > ipFitter( _kalTest ) ;
  //  //  KalTestFitter<KalTest::OrderOutgoing, KalTest::FitBackward, KalTest::PropagateToIP > ipFitter( _kalTest ) ;

  //  //FIXME: DEBUG - non ip fitter
  KalTestFitter<KalTest::OrderOutgoing, KalTest::FitBackward > ipFitter( _kalTest ) ;
  
  std::transform( cluList.begin(), cluList.end(), std::back_inserter( newKTracks ) , ipFitter  ) ;


  LCCollectionVec* kaltracks = new LCCollectionVec( LCIO::TRACK ) ;
  LCFlagImpl trkFlag(0) ;
  trkFlag.setBit( LCIO::TRBIT_HITS ) ;
  kaltracks->setFlag( trkFlag.getFlag()  ) ;
  
  std::transform( newKTracks.begin(), newKTracks.end(), std::back_inserter( *kaltracks ) , KalTrack2LCIO() ) ;

  //  //  std::transform( cluList.begin(), cluList.end(), std::back_inserter( *kaltracks ) , converter ) ;

  evt->addCollection( kaltracks , _outColName ) ;
  
 

  // //================================================================================================================ 
  //  //   merge track segments based on track parameters and errors ...
  //  //
  //  static const int merge_track_segments = true ;

  //  if( merge_track_segments ) {

  //    GenericHitVec<Track> trkVec ;
  //    GenericClusterVec<Track> trkCluVec ;
  //    LCCollectionVec* mergedTracks = new LCCollectionVec( LCIO::TRACK ) ;
  
  //    addToGenericHitVec( trkVec , kaltracks ,  AllwaysTrue()  ) ;

  //    //    TrackStateDistance trkMerge( 50. ) ;
  //    TrackCircleDistance trkMerge( 0.1 ) ; 

  //    cluster( trkVec.begin() , trkVec.end() , std::back_inserter( trkCluVec ), &trkMerge  , 2 ) ;


  //    streamlog_out( DEBUG4 ) << " ===== merged tracks - # cluster: " << trkCluVec.size()   << "  ============================== " << std::endl ;
    
  //    for( GenericClusterVec<Track>::iterator it= trkCluVec.begin() ; it != trkCluVec.end() ; ++it) {
      
  //      streamlog_out( DEBUG2 ) <<  myheader<Track>() << std::endl ;
      
  //      GenericCluster<Track>* trkClu = *it ;
      
  //      std::list<Track*> mergedTrk ;
  //      for( GenericCluster<Track>::iterator itC = trkClu->begin() ; itC != trkClu->end() ; ++ itC ){
	
  // 	streamlog_out( DEBUG2 ) << lcshort(  (*itC)->first ) << std::endl ; 
	
  // 	mergedTrk.push_back( (*itC)->first ) ; 
  //      }

  //      TrackImpl* trk = new TrackImpl ;
  //      Track* bestTrk = 0 ;
  //      double chi2Min = 99999999999999999. ;
  //      for( std::list<Track*>::iterator itML = mergedTrk.begin() ; itML != mergedTrk.end() ; ++ itML ){
	
  // 	const TrackerHitVec& hV = (*itML)->getTrackerHits() ;

  // 	for(unsigned i=0 ; i < hV.size() ; ++i){

  // 	  trk->addHit( hV[i] ) ;
  // 	  double chi2ndf = (*itML)->getChi2() / (*itML)->getNdf() ;

  // 	  if( chi2ndf < chi2Min ){
  // 	    bestTrk = (*itML) ;
  // 	    chi2Min = chi2ndf ;
  // 	  }
  // 	}
  //      }
  //      if( bestTrk != 0 ){ 

  // 	trk->setD0( bestTrk->getD0() ) ;
  // 	trk->setOmega( bestTrk->getOmega() ) ;
  // 	trk->setPhi( bestTrk->getPhi() ) ;
  // 	trk->setZ0( bestTrk->getZ0() ) ;
  // 	trk->setTanLambda( bestTrk->getTanLambda() ) ;
  // 	trk->setCovMatrix( bestTrk->getCovMatrix()  ) ;
  // 	// ...
	
  //      }
  //      else{
  // 	streamlog_out( ERROR ) << "   no best track found for merged tracks ... !? " << std::endl ; 
  //      }
  //      mergedTracks->addElement( trk )  ;

  //    }

  //    // add all tracks that have not been merged :
  //    for( GenericHitVec<Track>::iterator it = trkVec.begin(); it != trkVec.end() ;++it){

  //      if( (*it)->second == 0 ){

  // 	mergedTracks->addElement(  new TrackImpl( *dynamic_cast<TrackImpl*>( (*it)->first ) ) ) ;
  //      }
  //    }


  //    evt->addCollection( mergedTracks , "MergedKalTracks" ) ;
  //  }


  //------ register some debugging print funtctions for picking in CED :

  CEDPickingHandler::getInstance().registerFunction( LCIO::TRACKERHIT , &printTrackerHit ) ; 
  CEDPickingHandler::getInstance().registerFunction( LCIO::TRACK , &printTrackShort ) ; 
  //======================================================================================================

 
  //========== cleanup KalTracks ========
  std::for_each( ktracks.begin() , ktracks.end() , delete_ptr<KalTrack> ) ;

  //FIXME: memory leak - need for debugging....

  //  std::for_each( newKTracks.begin() , newKTracks.end() , delete_ptr<KalTrack> ) ;
  //=====================================



  //========  create collections of used and unused TPC hits ===========================================

  LCCollectionVec* usedHits = new LCCollectionVec( LCIO::TRACKERHIT ) ;   ;
  LCCollectionVec* unUsedHits = new LCCollectionVec( LCIO::TRACKERHIT ) ;   ;
  usedHits->setSubset() ;
  unUsedHits->setSubset() ;
  usedHits->reserve( h.size() ) ;
  unUsedHits->reserve( h.size() ) ;
  //  typedef GenericHitVec<TrackerHit>::iterator GHVI ;
  for( GHVI it = h.begin(); it != h.end() ;++it){
    if( (*it)->second != 0 ){
      usedHits->push_back( (*it)->first ) ;
    } else {
      unUsedHits->push_back( (*it)->first ) ;          
    }
  }
  for( GHVI it = hSmallR.begin(); it != hSmallR.end() ;++it){
    if( (*it)->second != 0 ){
      usedHits->push_back( (*it)->first ) ;
    } else {
      unUsedHits->push_back( (*it)->first ) ;          
    }
  }
  evt->addCollection( usedHits ,   "UsedTPCCluTrackerHits" ) ;
  evt->addCollection( unUsedHits , "UnUsedTPCCluTrackerHits" ) ;
  
  //========================================================================================================
  
  _nEvt ++ ;

  clock_t end = clock () ; 
  
  streamlog_out( DEBUG )  << "---  clustering time: " 
 			  <<  double( end - start ) / double(CLOCKS_PER_SEC) << std::endl  ;
  
}


/*************************************************************************************************/
void OuterRimSearch::check( LCEvent * evt ) { 
  /*************************************************************************************************/

  //std::string colName( "MergedKalTracks"  ) ;
  std::string colName(  _outColName ) ;


  bool checkForDuplicatePadRows =  false ; //true ;
  bool checkForMCTruth =  true ;

  bool checkForSplitTracks =  true ;   // WARNING: DEBUG only - this requires the kaltracks to not be deleted in processEvent !!!!!!!!! 


  streamlog_out( MESSAGE ) <<  " check called.... " << std::endl ;

  const gear::TPCParameters& gearTPC = Global::GEAR->getTPCParameters() ;
  const gear::PadRowLayout2D& pL = gearTPC.getPadLayout() ;


  //====================================================================================
  // check for duplicate padRows 
  //====================================================================================

  if( checkForDuplicatePadRows ) {

    LCCollectionVec* oddCol = new LCCollectionVec( LCIO::TRACK ) ;
    oddCol->setSubset( true ) ;
    // try iterator class ...

    LCIterator<Track> trIt( evt, colName ) ;
    while( Track* tr = trIt.next()  ){

      
      // check for duplicate layer numbers
      std::vector<int> hitsInLayer( pL.getNRows() ) ; 
      const TrackerHitVec& thv = tr->getTrackerHits() ;
      typedef TrackerHitVec::const_iterator THI ;
      for(THI it = thv.begin() ; it  != thv.end() ; ++it ) {
	TrackerHit* th = *it ;
	++ hitsInLayer.at( th->ext<HitInfo>()->layerID )   ;
      } 
      unsigned nHit = thv.size() ;
      unsigned nDouble = 0 ;
      for(unsigned i=0 ; i < hitsInLayer.size() ; ++i ) {
	if( hitsInLayer[i] > 1 ){
	  ++nDouble ;
	  streamlog_out( DEBUG4 ) << " &&&&&&&&&&&&&&&&&&&&&&&&&& duplicate hit in layer : " << i << std::endl ;
	}
      }
      if( double(nDouble) / nHit > _duplicatePadRowFraction ){
	//if( nDouble  > 0){
	streamlog_out( DEBUG4 ) << " oddTrackCluster found with "<< 100. * double(nDouble) / nHit 
				<< "% of double hits " << std::endl ;
	oddCol->addElement( tr ) ;
      }
    }
    evt->addCollection( oddCol , "OddCluTracks" ) ;
  }
  //====================================================================================
  // check Monte Carlo Truth via SimTrackerHits 
  //====================================================================================

  if( checkForMCTruth ) {
 

    LCCollectionVec* oddCol = new LCCollectionVec( LCIO::TRACK ) ;
    oddCol->setSubset( true ) ;

    LCCollectionVec* splitCol = new LCCollectionVec( LCIO::TRACK ) ;
    splitCol->setSubset( true ) ;
    
    typedef std::map<Track* , unsigned > TRKMAP ; 
    
    typedef std::map< MCParticle* , TRKMAP > MCPTRKMAP ;
    MCPTRKMAP mcpTrkMap ;
    
    typedef std::map< MCParticle* , unsigned > MCPMAP ;
    MCPMAP hitMap ;
    
    
    // if( streamlog_level( DEBUG4) )
    //   LCTOOLS::printTracks( evt->getCollection("KalTestTracks") ) ;


    LCIterator<Track> trIt( evt, colName  ) ;  
    //    "KalTestTracks" ) ;
    //    LCIterator<Track> trIt( evt, _outColName ) ;
    //    LCIterator<Track> trIt( evt, "TPCTracks" ) ;

    while( Track* tr = trIt.next()  ){
      
      MCPMAP mcpMap ;

      const TrackerHitVec& thv = tr->getTrackerHits() ;
      typedef TrackerHitVec::const_iterator THI ;

      // get relation between mcparticles and tracks
      for(THI it = thv.begin() ; it  != thv.end() ; ++it ) {

	TrackerHit* th = *it ;
	// FIXME:
	// we know that the digitizer puts the sim hit into the raw hit pointer
	// but of course the proper way is to go through the LCRelation ...
	SimTrackerHit* sh = (SimTrackerHit*) th->getRawHits()[0] ;
	MCParticle* mcp = sh->getMCParticle() ;

	
	hitMap[ mcp ] ++ ;   // count all hits from this mcp
	
	mcpMap[ mcp ]++ ;    // count hits from this mcp for this track
	
	mcpTrkMap[ mcp ][ tr ]++ ;  // map between mcp, tracks and hits
	
      } 

      // check for tracks with hits from several mcparticles
      unsigned nHit = thv.size() ;
      unsigned maxHit = 0 ; 
      for( MCPMAP::iterator it= mcpMap.begin() ;
	   it != mcpMap.end() ; ++it ){
	if( it->second  > maxHit ){
	  maxHit = it->second ;
	}
      }

      if( double(maxHit) / nHit < 0.99 ){ // What is acceptable here ???
	//if( nDouble  > 0){
	streamlog_out( MESSAGE ) << " oddTrackCluster found with only "
				 << 100.*double(maxHit)/nHit 
				 << "% of hits  form one MCParticle " << std::endl ;
	oddCol->addElement( tr ) ;
      }
    }
    evt->addCollection( oddCol , "OddMCPTracks" ) ;
    
    
    if( checkForSplitTracks ) {
      
      streamlog_out( DEBUG ) << " checking for split tracks - mcptrkmap size : " <<  mcpTrkMap.size() << std::endl ;
      
      // check for split tracks 
      for( MCPTRKMAP::iterator it0 = mcpTrkMap.begin() ; it0 != mcpTrkMap.end() ; ++it0){
	
	streamlog_out( DEBUG ) << " checking for split tracks - map size : " <<  it0->second.size() << std::endl ;
	
	
	if( it0->second.size() > 1 ) {
	  
	  
	  typedef std::list< EVENT::Track* > TL ;
	  TL trkList ;
	  
	  for( TRKMAP::iterator it1 = it0->second.begin() ; it1 != it0->second.end() ; ++it1){
	    
	    double totalHits = hitMap[ it0->first ]  ; // total hits for this track 
	    
	    double thisMCPHits = it1->second ;     //  hits from this mcp
	    
	    double ratio =  thisMCPHits / totalHits  ;
	    
	    streamlog_out( DEBUG ) << " checking for split tracks - ratio : " 
				   << thisMCPHits << " / " << totalHits << " = " << ratio << std::endl ;
	    
	    if( ratio > 0.03 && ratio < 0.95 ){
	      // split track
	      
	      splitCol->addElement( it1->first ) ; 
	      
	      trkList.push_back( it1->first ) ;
	    } 
	  }
	  // chi2 between split track segments :
	  // for( TRKMAP::iterator ist0 = it0->second.begin() ; ist0 != it0->second.end() ; ++ist0){
	    
	  //   KalTrack* sptrk0 = ist0->first->ext<KalTrackLink>() ; 
	    
	  //   TRKMAP::iterator ist0_pp = ist0 ;
	  //   ++ist0_pp ;

	  //   for( TRKMAP::iterator ist1 = ist0_pp ; ist1 != it0->second.end() ; ++ist1){
	  
	  //     KalTrack* sptrk1 = ist1->first->ext<KalTrackLink>() ; 
	      
	  //     double chi2 =  KalTrack::chi2( *sptrk0 ,  *sptrk1 ) ;
	      
	  //     streamlog_out( DEBUG4 ) << " *********************  chi2 between split tracks : "  << chi2 << std::endl 
	  // 			      << myheader< Track >() << std::endl 
	  // 			      << lcshort( ist0->first )  << std::endl 
	  // 			      << lcshort( ist1->first )	 << std::endl ; 
	      
	  //   }
	  // }



	  streamlog_out( DEBUG2 ) << " ------------------------------------------------------ " << std::endl ;
	  
	  for( TL::iterator it0 = trkList.begin() ; it0 != trkList.end() ; ++it0 ){
	    
	    
	    //	    KalTrack* trk0 = (*it0)->ext<KalTrackLink>() ; 
	    
	    HelixClass hel ;
	    hel.Initialize_Canonical( (*it0)->getPhi(),
				      (*it0)->getD0(),
				      (*it0)->getZ0(),
				      (*it0)->getOmega(),
				      (*it0)->getTanLambda(),
				      3.50 ) ;
	    
	    streamlog_out( DEBUG1 ) << hel.getXC() << "\t"
				    << hel.getYC() << "\t"
				    << hel.getRadius() << "\t" 
				    << hel.getTanLambda() << std::endl ; 
	    
	    
	    // streamlog_out( DEBUG1 ) << (*it0)->getPhi() << "\t"
	    // 			  << (*it0)->getD0()  << "\t"
	    // 			  << (*it0)->getOmega()  << "\t"
	    // 			  << (*it0)->getZ0()  << "\t"
	    // 			  << (*it0)->getTanLambda()  << "\t"
	    // 			  << std::endl ;
	    
	    //	    streamlog_out( DEBUG1 ) << " trk0 : " << *trk0 << std::endl ;
	    
	    // TL::iterator its = it0 ;
	    // ++its ;
	    
	    // for( TL::iterator it1 =  its ; it1 != trkList.end() ; ++it1 ){
	      
	    //   KalTrack* trk1 = (*it1)->ext<KalTrackLink>() ; 
	      
	    //   streamlog_out( DEBUG1 ) << "    - trk0 : " << *trk0 << std::endl ;
	    //   streamlog_out( DEBUG1 ) << "    - trk1 : " << *trk1 << std::endl ;
	      
	    //   double chi2 =  KalTrack::chi2( *trk0 ,  *trk1 ) ;
	      
	    //   streamlog_out( DEBUG1 ) << " +++++++++++++++++  chi2 between split tracks : " 
	    // 			      << trk0 << " - " << trk1 << " : " << chi2 << std::endl ; 
	      
	      
	    // }
	  }
	  
	}
      }
      evt->addCollection( splitCol , "SplitTracks" ) ;
    }

  }
  //====================================================================================

}


void OuterRimSearch::end(){ 
  
  //   std::cout << "OuterRimSearch::end()  " << name() 
  // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
  // 	    << std::endl ;
  
  delete _kalTest ;
}


//====================================================================================================
