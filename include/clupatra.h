#ifndef clupatra_h
#define clupatra_h

#include "assert.h"

#include "NNClusters.h"
#include "lcio.h"
#include "EVENT/TrackerHit.h"
#include "IMPL/TrackImpl.h"
#include "UTIL/Operators.h"
#include "LCRTRelations.h"
#include "TVTrackHit.h"
#include "KalTest.h"

#include "GEAR.h"

// copy_if (missing from STL)
template <class In, class Out, class Pred> Out copy_if(In first, In last, Out res, Pred p){
  
  while( first != last ){
    if( p( *first) ){
 
     *res++ = *first ;
    }
    ++first ;
  }
  return res ;
}


namespace clupatra{
  
  typedef TVTrackHit  FitHit ;
  

  /** Extend the gear::Vector3D w/ a c'tor that can be initiallized with an array such as positions from the lcio objects.
   */
  struct Vector3D : public gear::Vector3D{
    Vector3D( const double* v) : gear::Vector3D( v[0] , v[1] , v[2] ) {}
    Vector3D( const float* v) :  gear::Vector3D( v[0] , v[1] , v[2] ) {}
  };

  //-----------------------------------------------------------------------

  /** Wrapper hit class that holds hit related information needed for pat rec as well as pointers to the original lcio::TrackerHit 
   *  and a pointer to the hit used in the Fitter.
   */
  struct ClupaHit {
    
    ClupaHit() :layerID(-1), 
		zIndex(-1), 
		phiIndex(-1), 
		usedInTrack(false) , 
		chi2Residual(-1.) , 
		deltaChi2(-1.), 
		pNNHit(0) , 
		nNNHit(0), 
		pDist(-1.), 
		nDist(-1.),
		fitHit(0),
		lcioHit(0), 
		pos(0.,0.,0.) {}
    ~ClupaHit() { 
      delete fitHit ; 
    } 

    int layerID ;
    int zIndex ;
    int phiIndex ;
    bool usedInTrack ;
    double chi2Residual ;  
    double deltaChi2 ;
    ClupaHit* pNNHit ;
    ClupaHit* nNNHit ;
    double pDist ;
    double nDist ;
    FitHit* fitHit ;
    lcio::TrackerHit* lcioHit ;
    gear::Vector3D pos ;
  };
  inline lcio::TrackerHit* lcioHit( const ClupaHit* h) { return h->lcioHit ; }

  // allow assignment to lcio::TrackerHit (memory mgmt)
  struct HitInfo : lcrtrel::LCOwnedExtension<HitInfo, ClupaHit > {} ;
  

  inline int z_index( ClupaHit* h) { return h->zIndex ; } 



  // helper class to assign additional parameters to TrackerHit clusters
  struct ClusterInfoStruct{
    ClusterInfoStruct() : track(0){ }
    KalTrack* track ;
    gear::Vector3D nextXPoint ;
    int nextLayer ;
  } ;
  struct ClusterInfo : lcrtrel::LCOwnedExtension<ClusterInfo, ClusterInfoStruct> {} ;


  struct TrackInfoStruct{
    TrackInfoStruct() : nextXPoint(0.,0.,0.), nextLayer(-1) {}
    gear::Vector3D nextXPoint ;
    int nextLayer ;
    float zMin ;
    float zAvg ;
    float zMax ;

  } ;
  struct TrackInfo : lcrtrel::LCOwnedExtension<TrackInfo, TrackInfoStruct> {} ;

  //------------------ typedefs for generic hits and clusters ---------

  typedef GenericHit<ClupaHit>        GHit ;
  typedef GenericHitVec<ClupaHit>     GHitVec ;
  typedef GenericCluster<ClupaHit>    GCluster ;
  typedef GenericClusterVec<ClupaHit> GClusterVec ;

  //-------------  vector of left over hits per layer
  typedef std::list<GHit*> GHitList ;
  typedef std::vector< GHitList > GHitListVector ;

  //-------------------------------------------------------------------------------------
  
  // helper struct: returns true if |z| index is in ]z0,z1]
  struct ZIndexInRange{
    int _z0, _z1 ;
    ZIndexInRange(int z0, int z1 ) : _z0(z0), _z1(z1) {}
    bool operator()( GHit* hit ){
      int z = std::abs( hit->first->zIndex ) ;
      return  ( _z0 <= z && z < _z1 ) ;
    }
  } ;

  //-----------------------------------

  /** Add the generic hits from (First,Last) to the GHitListVector - vector needs to be initialized, e.g. with
   *     hLV.resize( KalTest::getMaxLayerIndex() ) 
   */
  //HitLayerID  tpcLayerID( _kalTest->indexOfFirstLayer( KalTest::DetID::TPC )  ) ;
  template <class First, class Last>
  void addToHitListVector( First first, Last last, GHitListVector& hLV, unsigned offset=0 ){
    
    while( first != last ){
      hLV[ (*first)->first->layerID + offset  ].push_back( *first )  ;
      ++first ;
    }
  }

  //-------------------------------------------------------------------------------------

  /** Try to add hits from hLV (hit lists per layer) to the cluster. The cluster needs to have a fitted KalTrack associated to it.
   *  Hits are added if the resulting delta Chi2 is less than dChiMax - a maxStep is the maximum number of steps (layers) w/o 
   *  successfully merging a hit.
   */
  void addHitsAndFilter( GCluster* clu, GHitListVector& hLV , double dChiMax, double chi2Cut, unsigned maxStep, 
			 bool backward=false) ; 

  //-------------------------------------------------------------------------------------
  /** Find the nearest neighbor hit in the next and previous layer(s) 
   */
  void  findNearestHits( GCluster& clu,  int maxLayerID, int maxStep=1) ;

  //-------------------------------------------------------------------------------------

  /** Force the hits in cluster clu into two clusters - only hits from layers with exactly two clusters are used.
   *  The algorithm assumes that the cluster consists of effectively two tracks - hits are split according to the 'hemisphere' 
   *  defined by the difference vector of the hit pair in the previous layer. 
   */
  //  std::pair<GCluster*, GCluster* > create_two_clusters( GCluster& clu,  GHitVec& hV,  int maxLayerID) ;
  void create_two_clusters( const GHitVec& hV, GClusterVec& cluVec,  int maxLayerID) ;
  
  
  /** Split the hits in this cluster up into three clusters. Hits are merged according to the smallest angle between the last hit in 
   *  the cluster and the new hits (as seen from the IP). Works for three prong tau decays in fwd. direction - might not for curved tracks 
   *  in central region ...
   */
  void create_three_clusters( const GHitVec& hV, GClusterVec& cluVec,  int maxLayerID) ;
  //-------------------------------------------------------------------------------------
  
  
  //----------------  delete helper
  template<class P>  void delete_ptr(P* p) { delete p;}
  //-------------------
  

  /** helper function to restrict the range of the azimuthal angle to ]-pi,pi]
   */
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
    //    double operator()( const TrackerHit* h, const gear::Vector3D& v1) {
    double operator()( const ClupaHit* h, const gear::Vector3D& v1) {


      //      gear::Vector3D v0( h->getPosition()[0] ,  h->getPosition()[1] ,  h->getPosition()[2] ) ;
      const gear::Vector3D& v0 = h->pos ;

      double sigsr =  sqrt( h->lcioHit->getCovMatrix()[0] + h->lcioHit->getCovMatrix()[2] ) ;
      double sigsz =  h->lcioHit->getCovMatrix()[5] ;
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
  // function to extract position for Kaltest:
  inline TVector3 hitPosition( GHit* h)  { 
    return h->first->pos ;
    // return TVector3( h->first->getPosition()[0],   
    // 		     h->first->getPosition()[1],
    // 		     h->first->getPosition()[2]  ) ; 
  }   

  // function to extract layerID from generic GHit:
  inline int hitLayerID( const GHit* h, int offset=0) { return  h->first->layerID + offset  ; } 

  // functor for layer ID
  class HitLayerID{
    int _off ;
    HitLayerID(){}
  public:
    HitLayerID( int off) : _off(off) {}
    int operator()(const GHit* h){ return hitLayerID( h, _off) ; } 
  } ;


  // // same for lcio::TrackerHit
  // int lcioLayerID( TrackerHit* h, int offset=0) { return  h->ext<HitInfo>()->layerID + offset  ; } 

  // // functor for layer ID
  // class LCIOLayerID{
  //   int _off ;
  //   LCIOLayerID(){}
  // public:
  //   LCIOLayerID( int off) : _off(off) {}
  //   int operator()(TrackerHit* h){ return lcioLayerID( h, _off) ; } 
  // } ;
  
  //-----------------

  struct LCIOTrackerHit{ lcio::TrackerHit* operator()( GHit* h) { return h->first->lcioHit  ; }   } ;

  //---------------------------------------------------
  // helper for sorting cluster wrt layerID
  template <bool SortDirection>
  struct LayerSort{
    bool operator()( const GHit* l, const GHit* r) {
      return hitLayerID( l ) < hitLayerID( r ) ; 
    }
  } ;
  template<>
  struct LayerSort<KalTest::OrderIncoming>{
    bool operator()( const GHit* l, const GHit* r) {
      return hitLayerID( r ) < hitLayerID( l ) ; 
    }
  } ;

  //------ ordering of KalTracks 
  struct KalTrackLengthSort {
    bool operator()( const KalTrack* t0, const KalTrack* t1) {
      return ( t0->getNHits() >= t1->getNHits() );
    }
  };



  //------------------------------
  //helpers for z ordering of hits
  struct TrackerHitCast{
    lcio::TrackerHit* operator()(lcio::LCObject* o) { return (lcio::TrackerHit*) o ; }
  };


  struct ZSort { 
    bool operator()( const ClupaHit* l, const ClupaHit* r) {
      return ( l->pos.z() < r->pos.z() ); 
    }
    // bool operator()( const TrackerHit* l, const TrackerHit* r) {
    //   return ( l->getPosition()[2] < r->getPosition()[2] );
    //}
  };
  


  //-------------------------------
  template <class T>
  void delete_elements(T* t) { delete t ; }

  //-------------------------------

  //-------------------------------------------------------------------------
  template <bool HitOrder, bool FitOrder, bool PropagateIP=false>

  struct KalTestFitter{

    KalTest* _kt ; 
  
    KalTestFitter(KalTest* k) : _kt( k ) {}
  
    KalTrack* operator() (GCluster* clu) {  
    
      KalTrack* trk = _kt->createKalTrack() ;


      // store mutual pointers between tracks and 'clusters'
      trk->setCluster<GCluster>( clu ) ;
      
      if( !  clu->ext<ClusterInfo>() )
       	clu->ext<ClusterInfo>() = new ClusterInfoStruct ;
      
      clu->ext<ClusterInfo>()->track = trk ;
      

      if( clu->size() < 3 ) 
	return trk ;
      
      static HitLayerID tpcLayerID( _kt->indexOfFirstLayer( KalTest::DetID::TPC )  )  ;
    
      clu->sort( LayerSort<HitOrder>() ) ;
    
    
      // need to reverse the order for incomming track segments (curlers)
      // assume particle comes from IP
      GHit* hf = clu->front() ;
      GHit* hb = clu->back() ;

      bool reverse_order = ( ( HitOrder ==  KalTest::OrderOutgoing ) ?    
			     ( std::abs( hf->first->pos.z() ) > std::abs( hb->first->pos.z()) + 3. )   :   
			     ( std::abs( hf->first->pos.z() ) < std::abs( hb->first->pos.z()) + 3. )   ) ;
    
      // reverse_order = false ;


   
      if( PropagateIP  && HitOrder == KalTest::OrderOutgoing ) {
      
	trk->addIPHit() ;
      }  
    
      // // ----- debug ----
      // std::set<int> layers ;
      // for( GCluster::iterator it=clu->begin() ; it != clu->end() ; ++it){
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
	
	for( GCluster::reverse_iterator it=clu->rbegin() ; it != clu->rend() ; ++it){   

	  TVTrackHit* h = (*it)->first->fitHit ;
	  if( h != 0 )
	    trk->addHit( h ) ; 
	  else
	    streamlog_out( DEBUG ) <<  "   hit not added ;-( " << std::endl ;
	}
	
      } else {
	
	for( GCluster::iterator it=clu->begin() ; it != clu->end() ; ++it){   

	  TVTrackHit* h = (*it)->first->fitHit ;
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

  // template <bool HitOrder, bool FitOrder, bool PropagateIP=false>

  // struct KalTestFitter{

  //   KalTest* _kt ; 
  
  //   KalTestFitter(KalTest* k) : _kt( k ) {}
  
  //   KalTrack* operator() (GCluster* clu) {  
    
  //     static HitLayerID tpcLayerID( _kt->indexOfFirstLayer( KalTest::DetID::TPC )  )  ;
    
  //     clu->sort( LayerSort<HitOrder>() ) ;
    
    
  //     // need to reverse the order for incomming track segments (curlers)
  //     // assume particle comes from IP
  //     Hit* hf = clu->front() ;
  //     Hit* hb = clu->back() ;

  //     bool reverse_order = ( ( HitOrder ==  KalTest::OrderOutgoing ) ?    
  // 			     ( std::abs( hf->first->pos.z() ) > std::abs( hb->first->pos.z() ) + 3  )  :
  // 			     ( std::abs( hf->first->pos.z() ) < std::abs( hb->first->pos.z() ) + 3  )  ) ;
    
  //     KalTrack* trk = _kt->createKalTrack() ;

  //     trk->setCluster<GCluster>( clu ) ;
    

  //     if( PropagateIP  && HitOrder == KalTest::OrderOutgoing ) {
      
  // 	trk->addIPHit() ;
  //     }  
    
  //     // // ----- debug ----
  //     // std::set<int> layers ;
  //     // for( GCluster::iterator it=clu->begin() ; it != clu->end() ; ++it){
  //     //   if( layers.find( tpcLayerID( *it ) ) != layers.end()  )
  //     // 	std::cout << " +++++++++++++++++++ duplicate layerID in addHits : " <<  tpcLayerID( *it ) << std::endl ;
  //     //   layers.insert( tpcLayerID( *it ) ) ;
  //     // }
  //     // // ---- end debug ----------
    
  //     if( reverse_order )
  // 	trk->addHits( clu->rbegin() , clu->rend() , hitPosition, tpcLayerID , LCIOTrackerHit() ) ; 
  //     else
  // 	trk->addHits( clu->begin() , clu->end() , hitPosition, tpcLayerID , LCIOTrackerHit() ) ; 
    

  //     if( PropagateIP  && HitOrder == KalTest::OrderIncoming ) {
      
  // 	trk->addIPHit() ;
  //     }  

  //     trk->fitTrack( FitOrder  ) ;
    
  //     return trk;
  //   }
  // };


  // struct KalTrack2LCIO{
  //   lcio::TrackImpl* operator() (KalTrack* trk) {  
  //     lcio::TrackImpl* lTrk = new lcio::TrackImpl ;
  //     trk->toLCIOTrack( lTrk  ) ;
  //     return lTrk ;
  //   }
  // };

   struct KalTrack2LCIO{
     lcio::TrackImpl* operator() (KalTrack* trk) {  

       lcio::TrackImpl* lTrk = new lcio::TrackImpl ;
       trk->toLCIOTrack( lTrk  ) ;
       
       lTrk->ext<TrackInfo>() = new TrackInfoStruct ;       
       
       // streamlog_out( DEBUG4 ) <<  "   ---- KalTrack2LCIO : "  <<  trk->getCluster<GCluster>()			     
       // 			 <<  " next xing point at layer: "   <<  trk->getCluster<GCluster>()->ext<ClusterInfo>()->nextLayer
       // 			 << " : " <<  trk->getCluster<GCluster>()->ext<ClusterInfo>()->nextXPoint ;
       
       if( streamlog_level( DEBUG ) ){
	 lTrk->ext<TrackInfo>()->nextLayer   =  trk->getCluster<GCluster>()->ext<ClusterInfo>()->nextLayer ;
	 lTrk->ext<TrackInfo>()->nextXPoint  =  trk->getCluster<GCluster>()->ext<ClusterInfo>()->nextXPoint ;
       }
       
       
       // compute z-extend of this track segment
       const lcio::TrackerHitVec& hv = lTrk->getTrackerHits() ;
       float zMin =  1e99 ;
       float zMax = -1e99 ;
       float zAvg =  0. ;
       for(unsigned i=0; i < hv.size() ; ++i ){
	 
	 float z = hv[i]->getPosition()[2] ;
	 
	 if( z < zMin )   zMin = z ;
	 if( z > zMax )   zMax = z ;
	 
	 zAvg += z ;
       }
       zAvg /= hv.size() ;
       
       // if( hv.size() >  1 ) {
       // 	 zMin = hv[            0  ]->getPosition()[2] ;
       // 	 zMax = hv[ hv.size() -1  ]->getPosition()[2] ;
       // }

       if( zMin > zMax ){ // swap 
	 float d = zMax ;
	 zMax = zMin ;
	 zMin = d  ;
       }

       lTrk->ext<TrackInfo>()->zMin = zMin ;
       lTrk->ext<TrackInfo>()->zMax = zMax ;
       lTrk->ext<TrackInfo>()->zAvg = zAvg ;
       
       
       return lTrk ;
     }
   };
  //-------------------------------------------------------------------------
  template <class T>
  class RCut {
  public:
    RCut( double rcut ) : _rcut( rcut ) {}  
    
    bool operator() (T* hit) {  
      
      return ( ( hit->pos.rho() > _rcut )   ||
	       ( hit->pos.z()   > (500. + _rcut ) ) )  ; 
	       
	       // return  ( (std::sqrt( hit->getPosition()[0]*hit->getPosition()[0] +
      // 			    hit->getPosition()[1]*hit->getPosition()[1] )   > _rcut )   ||
      // 		( std::abs( hit->getPosition()[2] ) > (500. + _rcut ) )
      // 		); 
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
      return ( ( hit->pos.rho() <= _rcut )  &&
	       ( hit->pos.z()   <= (500. + _rcut ) ) ) ; 

      // return (  ( std::sqrt( hit->getPosition()[0]*hit->getPosition()[0] +
      // 			     hit->getPosition()[1]*hit->getPosition()[1] )   <= _rcut )   &&
      // 		(  std::abs( hit->getPosition()[2] ) <= (500. + _rcut ) )
      // 		) ;

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

    bool operator()(const GCluster* cl) const {
 
      // check for duplicate layer numbers
      std::vector<int> hLayer( _N )  ; 
      typedef GCluster::const_iterator IT ;

      unsigned nHit = 0 ;
      for(IT it=cl->begin() ; it != cl->end() ; ++it ) {
	++ hLayer[ (*it)->first->layerID ]   ;
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
  //TODO: create a faster predicate for no duplicate pad rows ....

  //---------------------------------------------------------------------------------
  
  /** Predicate class for 'distance' of NN clustering based on previous nearest neighbor search.
   */
  class NearestHitDistance{
    typedef ClupaHit HitClass ;
    typedef double PosType ;
  public:
    
    /** Required typedef for cluster algorithm 
     */
    typedef HitClass hit_type ;
    
    /** C'tor takes merge distance */
    NearestHitDistance(float dCut) : _dCutSquared( dCut*dCut ) , _dCut(dCut)  {} 
    
    /** Merge condition: true if distance  is less than dCut given in the C'tor.*/ 
    inline bool mergeHits( GenericHit<HitClass>* h0, GenericHit<HitClass>* h1){
      
      //      if( std::abs( h0->Index0 - h1->Index0 ) > 1 ) return false ;
      
      int l0 =  h0->first->layerID ;
      int l1 =  h1->first->layerID ;
      
      if(  std::abs( l0 - l1 ) != 1 ) 
	return false ;
      

      bool ret = ( l0 < l1 ? 
		   h0->first->nNNHit ==  h1->first  &&  h1->first->pNNHit ==  h0->first  :
		   h0->first->pNNHit ==  h1->first  &&  h1->first->nNNHit ==  h0->first  ) ;

      streamlog_out( DEBUG4 ) << "--- NearestHitDistance::mergeHits: "  << std::endl 
			      << " h0: " << h0->first 
			      << " h0.pNNHit " << h0->first->pNNHit
			      << " h0.nNNHit " << h0->first->nNNHit << std::endl 
			      << " h1: " << h1->first 
			      << " h1.pNNHit " << h1->first->pNNHit
			      << " h1.nNNHit " << h1->first->nNNHit 
			      << "  ---> " <<  ( h0->first->nNNHit ==  h1->first  &&  h1->first->pNNHit ==  h0->first )
			      << "  , " <<     ( h0->first->pNNHit ==  h1->first  &&  h1->first->nNNHit ==  h0->first )
			      << std::endl ;

      return ret ;
    }
    
  protected:
    NearestHitDistance() ;
    float _dCutSquared ;
    float _dCut ;
  } ;
  
  //---------------------------------------------------------------------------------

  /** Predicate class for 'distance' of NN clustering.
   */
  //template <class HitClass, typename PosType > 
  class HitDistance{
    typedef ClupaHit HitClass ;
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
    
      //     int l0 =  h0->first->layerID ;
      //     int l1 =  h1->first->layerID ;

      //     //------- don't merge hits from same layer !
      //     if( l0 == l1 )
      //       return false ;

      if( h0->first->layerID == h1->first->layerID )
	return false ;

      // const PosType* pos0 =  h0->first->getPosition() ;
      // const PosType* pos1 =  h1->first->getPosition() ;
    
      // return 
      // 	( pos0[0] - pos1[0] ) * ( pos0[0] - pos1[0] ) +
      // 	( pos0[1] - pos1[1] ) * ( pos0[1] - pos1[1] ) +
      // 	( pos0[2] - pos1[2] ) * ( pos0[2] - pos1[2] )   
      // 	< _dCutSquared ;

      return ( h0->first->pos - h1->first->pos).r2()  < _dCutSquared ;

    }
  
  protected:
    HitDistance() ;
    float _dCutSquared ;
    float _dCut ;
  } ;

  class HitDistance_2{
    typedef ClupaHit HitClass ;
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

      int l0 =  h0->first->layerID ;
      int l1 =  h1->first->layerID ;
      
      
      //------- don't merge hits from same layer !
      if( l0 == l1 )
	return false ;


      // const PosType* pos0 =  h0->first->getPosition() ;
      // const PosType* pos1 =  h1->first->getPosition() ;
    
      return inRange<-2,2>(  l0 - l1 )  &&  
	( h0->first->pos - h1->first->pos ).r2()  
	< _dCutSquared ;

      // const PosType* pos0 =  h0->first->getPosition() ;
      // const PosType* pos1 =  h1->first->getPosition() ;
    
      // return inRange<-2,2>(  l0 - l1 )  &&  
      // 	( pos0[0] - pos1[0] ) * ( pos0[0] - pos1[0] ) +
      // 	( pos0[1] - pos1[1] ) * ( pos0[1] - pos1[1] ) +
      // 	( pos0[2] - pos1[2] ) * ( pos0[2] - pos1[2] )   
      // 	< _dCutSquared ;
    }
  
  protected:
    HitDistance_2() ;
    float _dCutSquared ;
    float _dCut ;
  } ;


  template <class T>
  struct LCIOTrack{
  
    lcio::Track* operator() (GenericCluster<T>* c) {  
    
      lcio::TrackImpl* trk = new lcio::TrackImpl ;
    
      double e = 0.0 ;
      int nHit = 0 ;
      for( typename GenericCluster<T>::iterator hi = c->begin(); hi != c->end() ; hi++) {
      
	trk->addHit(  (*hi)->first->lcioHit ) ;
	e += (*hi)->first->lcioHit->getEDep() ;
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
  const std::string & myheader(){    return lcio::header(*(T*)(0));  }


  inline void printTrackShort(const lcio::LCObject* o){

    lcio::Track* trk = const_cast<lcio::Track*> ( dynamic_cast<const lcio::Track*> (o) ) ; 
    
    if( trk == 0 ) {
      
      streamlog_out( ERROR ) << "  printTrackShort : dynamic_cast<Track*> failed for LCObject : " << o << std::endl ;
      return  ;
    }
    
    streamlog_out( MESSAGE ) << myheader<lcio::Track>()  
			     << lcio::lcshort( trk ) << std::endl  ;
    
    
    double r0 = 1. / trk->getOmega() ;
    double d0 = trk->getD0() ;
    double p0 = trk->getPhi() ;
  
    double x0 = ( r0 - d0 ) * sin( p0 ) ;
    double y0 = ( d0 - r0 ) * cos( p0 ) ;
  
    streamlog_out( MESSAGE ) << " circle: r = " << r0 << ", xc = " << x0 << " , yc = " << y0 << std::endl ;
    //    if( streamlog_level( DEBUG ) ){
    if( trk->ext<TrackInfo>() )
      streamlog_out( MESSAGE ) << "  next xing point at layer  "  << trk->ext<TrackInfo>()->nextLayer << " : " 
			       <<  trk->ext<TrackInfo>()->nextXPoint 
			       << " zMin " << trk->ext<TrackInfo>()->zMin  
			       << " zMax " << trk->ext<TrackInfo>()->zMax 
			       << " zAvg " << trk->ext<TrackInfo>()->zAvg 
			       <<  std::endl ;
      // }
   
  }

  inline void printTrackerHit(const lcio::LCObject* o){
  
    using namespace UTIL ;

    lcio::TrackerHit* hit = const_cast<lcio::TrackerHit*> ( dynamic_cast<const lcio::TrackerHit*> (o) ) ; 
  
    if( hit == 0 ) {
    
      streamlog_out( ERROR ) << "  printTrackerHit : dynamic_cast<TrackerHit*> failed for LCObject : " << o << std::endl ;
      return  ;
    }
  
    // streamlog_out( MESSAGE ) << *hit << std::endl 
    // 			     << " err: rPhi" <<  sqrt( hit->getCovMatrix()[0] + hit->getCovMatrix()[2] ) 
    // 			     << " z :  " <<   hit->getCovMatrix()[5] << std::endl 
    // 			     << " chi2 residual to best matching track : " << hit->ext<HitInfo>()->chi2Residual << std::endl ;
    streamlog_out( MESSAGE ) << *hit << std::endl 
			     << " err: rPhi" <<  sqrt( hit->getCovMatrix()[0] + hit->getCovMatrix()[2] ) 
			     << " z :  " <<   hit->getCovMatrix()[5] << std::endl 
			     << " chi2 residual to best matching track : " << hit->ext<HitInfo>()->chi2Residual
			     << " delta chi2 to best matching track : " << hit->ext<HitInfo>()->deltaChi2
			     << " zIndex: " << hit->ext<HitInfo>()->zIndex 
			     << " layerID: " << hit->ext<HitInfo>()->layerID 
			     << std::endl ;
    

  }


  /** Predicate class for track merging with NN clustering.
   */
  class TrackStateDistance{
    typedef lcio::Track HitClass ;

  public:
  
    /** Required typedef for cluster algorithm 
     */
    typedef HitClass hit_type ;
  
    /** C'tor takes merge distance */
    TrackStateDistance(float dCut) : _dCutSquared( dCut*dCut ) , _dCut(dCut)  {} 


    /** Merge condition: true if distance  is less than dCut given in the C'tor.*/ 
    inline bool mergeHits( GenericHit<HitClass>* h0, GenericHit<HitClass>* h1){
    
      if( std::abs( h0->Index0 - h1->Index0 ) > 1 ) return false ;
    
      lcio::Track* trk0 = h0->first ;
      lcio::Track* trk1 = h1->first ;

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


  //=======================================================================================

  /** helper class for merging track segments, based on circle (and tan lambda) */

  class TrackCircleDistance{

    typedef lcio::Track HitClass ;

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
      
      // if( std::abs( h0->Index0 - h1->Index0 ) > 1 ) return false ;
      
      lcio::Track* trk0 = h0->first ;
      lcio::Track* trk1 = h1->first ;
      
      const TrackInfoStruct* ti0 =  trk0->ext<TrackInfo>() ;
      const TrackInfoStruct* ti1 =  trk1->ext<TrackInfo>() ;
      
      
      // don't allow  overlaps in z
      float epsilon = 0. ; 
      if(  ti0->zAvg > ti1->zAvg ){
	
      	if( ti1->zMax > ( ti0->zMin + epsilon )  )
      	  return false ;
	
      } else {
	
      	if( ti0->zMax > ( ti1->zMin + epsilon ) )
      	  return false ;
	
      }
      // //------- dont' merge complete tracks: ------------------------
      // unsigned nh0 = trk0->getTrackerHits().size() ;
      // unsigned nh1 = trk1->getTrackerHits().size() ;
      // if( nh0 > 220 ) return false ;
      // if( nh1 > 220 ) return false ;
      // //------------------------------------------------------
      
      // KalTrack* ktrk0 = h0->first->ext<KalTrackLink>() ; 
      // KalTrack* ktrk1 = h1->first->ext<KalTrackLink>() ; 
      
      double tl0 = std::abs( trk0->getTanLambda() ) ;
      double tl1 = std::abs( trk1->getTanLambda() ) ;

      double dtl = 2. * std::abs( tl0 - tl1 ) / ( tl0 + tl1 ) ;
      if(  dtl > 2.  * _dCut  &&  std::abs( tl0 + tl1 ) > 1.e-2  )
      	return false ;
      // for very steep tracks (tanL < 0.001 ) tanL might differ largely for curlers due to multiple scattering

      double r0 = 1. / trk0->getOmega()   ;
      double r1 = 1. / trk1->getOmega()  ;

      double r0abs = std::abs( r0 ) ; 
      double r1abs = std::abs( r1 ) ; 


      double d0 = trk0->getD0() ;
      double d1 = trk1->getD0() ;

      double z0 = trk0->getZ0() ;
      double z1 = trk1->getZ0() ;

      double rIP = 20. ; 
      // don't merge tracks that come from an area of 20 mm around the IP
      if(  std::abs( d0 ) < rIP &&  std::abs( z0 ) < rIP  &&
      	   std::abs( d1 ) < rIP &&  std::abs( z1 ) < rIP     )
      	return false ;

      double p0 = trk0->getPhi() ;
      double p1 = trk1->getPhi() ;

      double x0 = ( r0 - d0 ) * sin( p0 ) ;
      double x1 = ( r1 - d1 ) * sin( p1 ) ;

      double y0 = ( d0 - r0 ) * cos( p0 ) ;
      double y1 = ( d1 - r1 ) * cos( p1 ) ;
    
      double dr = 2. * std::abs( r0abs - r1abs )  / (r0abs + r1abs )  ;

      double distMS = sqrt ( ( x0 - x1 ) * ( x0 - x1 ) + ( y0 - y1 ) * ( y0 - y1 )  ) ;
    
    
      //      return ( dr < DRMAX && distMS < _dCut * std::abs( r0 )  ) ;
      return ( dr < _dCut * std::abs( r0 )  &&  distMS < _dCut * std::abs( r0 )  ) ;

    }
  
  protected:
    float _dCutSquared ;
    float _dCut ;
  } ; 









}
#endif
