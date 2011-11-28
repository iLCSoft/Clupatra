#ifndef clupatra_new_h
#define clupatra_new_h

#include <cmath>
#include <time.h>
#include <math.h>
#include <sstream>
#include "assert.h"

#include "NNClusterer.h"

#include "lcio.h"
#include "EVENT/TrackerHit.h"
#include "IMPL/TrackImpl.h"
#include "UTIL/Operators.h"
#include "UTIL/ILDConf.h"

#include "GEAR.h"

#include "MarlinTrk/IMarlinTrack.h"
#include "MarlinTrk/IMarlinTrkSystem.h"

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"


namespace clupatra_new{
  
  /** Small wrapper extension of the LCIO Hit
   */
  struct ClupaHit {
    
    ClupaHit() :layer(-1), 
		zIndex(-1), 
		phiIndex(-1), 
		lcioHit(0), 
		pos(0.,0.,0.) {}
    int layer ;
    int zIndex ;
    int phiIndex ;
    lcio::TrackerHit* lcioHit ;
    gear::Vector3D pos ;

  };
  
  //  inline lcio::TrackerHit* lcioHit( const ClupaHit* h) { return h->lcioHit ; }

  
//------------------ typedefs for elements and clusters ---------

  typedef nnclu::NNClusterer< ClupaHit > Clusterer ;
  
  typedef Clusterer::element_type Hit ;
  typedef Clusterer::cluster_type CluTrack ;

  typedef Clusterer::element_vector HitVec ;
  typedef Clusterer::cluster_vector CluTrackVec ;
  
  typedef std::list<Hit*>        HitList ;
  typedef std::vector< HitList > HitListVector ;
  

  // typedef GenericHitVec<ClupaHit>      GHitVec ;
  // typedef GenericClusterVec<ClupaHit>  GClusterVec ;
  
  //------------------------------------------------------------------------------------------

  struct MarTrk : lcrtrel::LCExtension<MarTrk, MarlinTrk::IMarlinTrack> {} ;

  //------------------------------------------------------------------------------------------
  
  /** Simple predicate class for computing an index from N bins of the z-coordinate of LCObjects
   *  that have a float/double* getPostion() method.
   */
  class ZIndex{
  public:
    /** C'tor takes zmin and zmax - NB index can be negative and larger than N */
    ZIndex( float zmin , float zmax , int n ) : _zmin( zmin ), _zmax( zmax ) , _N (n) {}  

    template <class T>
    inline int operator() (T* hit) {  
      
      return (int) std::floor( ( hit->getPosition()[2] - _zmin ) / ( _zmax - _zmin ) * _N ) ; 
    }
    
    inline int index( double z) {  return  (int) std::floor( ( z - _zmin ) / ( _zmax - _zmin ) * _N ) ;  } 

  protected:
    ZIndex() {} ;
    float _zmin ;
    float _zmax ;
    int _N ;
  } ;
  
  //------------------------------------------------------------------------------------------

  struct ZSort { 
    inline bool operator()( const Hit* l, const Hit* r) {      
      return ( l->first->pos.z() < r->first->pos.z() ); 
    }
  };
  

  //------------------------------------------------------------------------------------------

  struct PtSort {  // sort tracks wtr to pt - largest first
    inline bool operator()( const lcio::LCObject* l, const lcio::LCObject* r) {      
      return ( std::abs( ( (const lcio::Track*) l )->getOmega() ) < std::abs( ( (const lcio::Track*) r )->getOmega() )  );  // pt ~ 1./omega  
    }
  };

  //------------------------------------------------------------------------------------------


  /** Add the elements (Hits) from (First,Last) to the HitListVector - vector needs to be initialized, e.g. with
   *  hLV.resize( MAX_LAYER_INDEX ) 
   */
  template <class First, class Last>
  void addToHitListVector( First first, Last last, HitListVector& hLV ){
    
    while( first != last ){
      hLV[ (*first)->first->layer ].push_back( *first )  ;
      ++first ;
    }
  }

  //------------------------------------------------------------------------------------------

  /** Predicate class for 'distance' of NN clustering. */
  class HitDistance{
  public:

    HitDistance(float dCut) : _dCutSquared( dCut*dCut ) {} 

    /** Merge condition: true if distance  is less than dCut */ 
    inline bool operator()( Hit* h0, Hit* h1){
    
      if( std::abs( h0->Index0 - h1->Index0 ) > 1 ) return false ;
    
      if( h0->first->layer == h1->first->layer )
	return false ;

      return ( h0->first->pos - h1->first->pos).r2()  < _dCutSquared ;
    }
  protected:
    HitDistance() ;
    float _dCutSquared ;
  } ;
  
  
  // /** Predicate class for 'distance' of NN clustering. */

  // struct HitDistance{  float _dCutSquared ;
  //   HitDistance(float dCut) : _dCutSquared( dCut*dCut ) {} 
    
  //   inline bool operator()( Hit* h0, Hit* h1){
      
  //     if( h0->first->layer == h1->first->layer )
  // 	return false ;
      
  //     return ( h0->first->pos - h1->first->pos).r2()  
  // 	< _dCutSquared ;
  //   }
  // } ;
  


  //------------------------------------------------------------------------------------------
  
  struct LCIOTrackConverter{
    
    lcio::Track* operator() (CluTrack* c) {  
    
      lcio::TrackImpl* trk = new lcio::TrackImpl ;
 
   
      double e = 0.0 ;
      int nHit = 0 ;
      for( CluTrack::iterator hi = c->begin(); hi != c->end() ; hi++) {
      
	trk->addHit(  (*hi)->first->lcioHit ) ;
	e += (*hi)->first->lcioHit->getEDep() ;
	nHit++ ;
      }

      MarlinTrk::IMarlinTrack* mtrk = c->ext<MarTrk>()  ;

      trk->ext<MarTrk>()  = mtrk ;

      trk->setdEdx( e/nHit ) ;
      trk->subdetectorHitNumbers().resize( 2 * lcio::ILDDetID::ETD ) ;
      trk->subdetectorHitNumbers()[ 2*lcio::ILDDetID::TPC - 1 ] =  nHit ;  // ??used in fit?? 
      trk->subdetectorHitNumbers()[ 2*lcio::ILDDetID::TPC - 2 ] =  nHit ;  
      
      if( mtrk != 0 ){
	
	lcio::TrackStateImpl* tsIP =  new lcio::TrackStateImpl ;
	lcio::TrackStateImpl* tsFH =  new lcio::TrackStateImpl ;
	lcio::TrackStateImpl* tsLH =  new lcio::TrackStateImpl ;
	//	lcio::TrackStateImpl* tsCA =  new lcio::TrackStateImpl ;
	
	tsIP->setLocation(  lcio::TrackState::AtIP ) ;
	tsFH->setLocation(  lcio::TrackState::AtFirstHit ) ;
	tsLH->setLocation(  lcio::TrackState::AtLastHit) ;
	//	tsCA->setLocation(  lcio::TrackState::AtCalo ) ;
	
	double chi2 ;
	int ndf  ;
	int code ;
	
	Hit* hf = c->front() ;
	Hit* hb = c->back() ;
	bool reverse_order =   ( std::abs( hf->first->pos.z() ) > std::abs( hb->first->pos.z()) + 3. ) ;

	lcio::TrackerHit* fHit =  ( reverse_order ?  hb->first->lcioHit  :  hf->first->lcioHit ) ;
	lcio::TrackerHit* lHit =  ( reverse_order ?  hf->first->lcioHit  :  hb->first->lcioHit ) ;
	
	// ======= get TrackState at first hit  ========================
	
	code = mtrk->getTrackState( fHit, *tsFH, chi2, ndf ) ;
	
	if( code != MarlinTrk::IMarlinTrack::success ){
	  
 	  streamlog_out( ERROR ) << "  >>>>>>>>>>> LCIOTrackConverter :  could not get TrackState at first Hit !!?? " << std::endl ; 
	}
	
	// ======= get TrackState at last hit  ========================
	code = mtrk->getTrackState( lHit, *tsLH, chi2, ndf ) ;
	
	if( code != MarlinTrk::IMarlinTrack::success ){
	  
 	  streamlog_out( ERROR ) << "  >>>>>>>>>>> LCIOTrackConverter :  could not get TrackState at last Hit !!?? " << std::endl ; 
	}
	
	// ======= get TrackState at calo face  ========================
	//
	//     need ecal face as material layers in KalDet ....
	//     or define planes and use helix utilitites ....
	//
	// ======= get TrackState at IP ========================
	
	const gear::Vector3D ipv( 0.,0.,0. );
	
	// fg: propagate is quite slow  and not really needed for the TPC
	// int ret = mtrk->propagate( ipv, *tsIP, chi2, ndf ) ;
	code = mtrk->extrapolate( ipv, *tsIP, chi2, ndf ) ;
	
	if( code != MarlinTrk::IMarlinTrack::success ){
	  
 	  streamlog_out( ERROR ) << "  >>>>>>>>>>> LCIOTrackConverter :  could not extrapolate TrackState to IP !!?? " << std::endl ; 
	}
	
	// trk->trackStates().push_back( tsIP ) ;
	// trk->trackStates().push_back( tsFH ) ;
	// trk->trackStates().push_back( tsLH ) ;
	// //	trk->trackStates().push_back( tsCA ) ;

	trk->addTrackState( tsIP ) ;
	trk->addTrackState( tsFH ) ;
	trk->addTrackState( tsLH ) ;
	//	trk->addTrackState( tsCA ) ;
	
	trk->setChi2( chi2 ) ;
	trk->setNdf( ndf ) ;
	

      } else {
	  
	//	streamlog_out( ERROR ) << "  >>>>>>>>>>> LCIOTrackConverter::operator() (CluTrack* c)  :  no MarlinTrk::IMarlinTrack* found for cluster !!?? " << std::endl ; 
      }


      return trk ;
    }

  } ;

  //------------------------------------------------------------------------------------------
  
  
  /** Predicate class for identifying clusters with duplicate pad rows - returns true
   *  if the fraction of duplicate hits is larger than 'fraction'.
   */
  struct DuplicatePadRows{

    unsigned _N ;
    float _f ; 
    DuplicatePadRows(unsigned nPadRows, float fraction) : _N( nPadRows), _f( fraction )  {}

    bool operator()(const CluTrack* cl) const {
 
      // check for duplicate layer numbers
      std::vector<int> hLayer( _N )  ; 
      typedef CluTrack::const_iterator IT ;

      unsigned nHit = 0 ;
      for(IT it=cl->begin() ; it != cl->end() ; ++it ) {
	++ hLayer[ (*it)->first->layer]   ;
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

  //------------------------------------------------------------------------------------------

  // helper for sorting cluster wrt layerID
  struct LayerSortOut{
    bool operator()( const Hit* l, const Hit* r) { return l->first->layer < r->first->layer ; }
  } ;
  struct LayerSortIn{
    bool operator()( const Hit* l, const Hit* r) { return l->first->layer > r->first->layer ; }
  } ;
  
  


  //------------------------------------------------------------------------------------------

  struct IMarlinTrkFitter{
    
    MarlinTrk::IMarlinTrkSystem* _ts ; 
    
    IMarlinTrkFitter(MarlinTrk::IMarlinTrkSystem* ts) : _ts( ts ) {}
    
    MarlinTrk::IMarlinTrack* operator() (CluTrack* clu) {  
      
      MarlinTrk::IMarlinTrack* trk = _ts->createTrack();
      
      if( clu->empty()  ){

	streamlog_out( ERROR ) << " IMarlinTrkFitter::operator() : cannot fit empty cluster track ! " << std::endl ;

	return trk ;
      }

      clu->ext<MarTrk>() = trk ;
      
      clu->sort( LayerSortOut() ) ;
            
      // need to reverse the order for incomming track segments (curlers)
      // assume particle comes from IP
      Hit* hf = clu->front() ;
      Hit* hb = clu->back() ;
      
      bool reverse_order =   ( std::abs( hf->first->pos.z() ) > std::abs( hb->first->pos.z()) + 3. ) ;
      
      if( reverse_order ){
	
	for( CluTrack::reverse_iterator it=clu->rbegin() ; it != clu->rend() ; ++it){   
	  
	  trk->addHit( (*it)->first->lcioHit  ) ; 
	  
	  //	  streamlog_out( DEBUG ) <<  "   hit  added  " <<  *(*it)->first->lcioHit   << std::endl ;
	}
	
	trk->initialise( MarlinTrk::IMarlinTrack::forward ) ;

      } else {
	
	for( CluTrack::iterator it=clu->begin() ; it != clu->end() ; ++it){   
	  
	  trk->addHit( (*it)->first->lcioHit   ) ; 
	  
	  //streamlog_out( DEBUG ) <<  "   hit  added  "<<  *(*it)->first->lcioHit   << std::endl ;
	}

	trk->initialise( MarlinTrk::IMarlinTrack::backward ) ;
      }


      trk->fit() ;

      return trk;
    }
  };

  //-------------------------------------------------------------------------------------

  /** Try to add hits from hLV (hit lists per layer) to the cluster. The cluster needs to have a fitted KalTrack associated to it.
   *  Hits are added if the resulting delta Chi2 is less than dChiMax - a maxStep is the maximum number of steps (layers) w/o 
   *  successfully merging a hit.
   */
  void addHitsAndFilter( CluTrack* clu, HitListVector& hLV , double dChiMax, double chi2Cut, unsigned maxStep, ZIndex& zIndex,  bool backward=false) ; 
  //------------------------------------------------------------------------------------------
  
  /** Try to add a hit from the given HitList in layer of subdetector to the track.
   *  A hit is added if the resulting delta Chi2 is less than dChiMax.
   */
  bool addHitAndFilter( int detectorID, int layer, CluTrack* clu, HitListVector& hLV , double dChiMax, double chi2Cut) ; 
  
  //------------------------------------------------------------------------------------------

  /** Returns the number of rows where cluster clu has i hits in mult[i] for i=1,2,3,4,.... -
   *  mult[0] counts all rows that have hits
   */
  void getHitMultiplicities( CluTrack* clu, std::vector<int>& mult ) ;

  //------------------------------------------------------------------------------------------

  /** Split the cluster into two clusters.
   */
  void create_two_clusters( Clusterer::cluster_type& clu, Clusterer::cluster_list& cluVec ) ;

  //------------------------------------------------------------------------------------------

  /** Split the cluster into three clusters.
   */
  void create_three_clusters( Clusterer::cluster_type& clu, Clusterer::cluster_list& cluVec ) ;




  //=======================================================================================

  struct TrackInfoStruct{  TrackInfoStruct() : zMin(0.), zAvg(0.), zMax(0.) {}
    float zMin ;
    float zAvg ;
    float zMax ;
  } ;
  struct TrackInfo : lcrtrel::LCOwnedExtension<TrackInfo, TrackInfoStruct> {} ;

  /** Helper class to compute track segment properties.
   */
  struct ComputeTrackerInfo{ 

    void operator()( lcio::LCObject* o ){ 

      lcio::Track* lTrk  = (lcio::Track*) o ; 

      // compute z-extend of this track segment
      const lcio::TrackerHitVec& hv = lTrk->getTrackerHits() ;

      float zMin =  1e99 ;
      float zMax = -1e99 ;
      float zAvg =  0. ;
      // for(unsigned i=0; i < hv.size() ; ++i ){
	
      // 	float z = hv[i]->getPosition()[2] ;
	
      // 	if( z < zMin )   zMin = z ;
      // 	if( z > zMax )   zMax = z ;
	
      // 	zAvg += z ;
      // }
      // zAvg /= hv.size() ;
      
      if( hv.size() >  1 ) {
	zMin = hv[            0  ]->getPosition()[2] ;
	zMax = hv[ hv.size() -1  ]->getPosition()[2] ;
	zAvg = ( zMax + zMin ) / 2. ;
      }
      
      if( zMin > zMax ){ // swap 
	float d = zMax ;
	zMax = zMin ;
	zMin = d  ;
      }

      if( ! lTrk->ext<TrackInfo>() )
	lTrk->ext<TrackInfo>() =  new TrackInfoStruct ;

      lTrk->ext<TrackInfo>()->zMin = zMin ;
      lTrk->ext<TrackInfo>()->zMax = zMax ;
      lTrk->ext<TrackInfo>()->zAvg = zAvg ;
    }
  } ;
  //=======================================================================================
  
  /** helper class for merging track segments, based on circle (and tan lambda) */
  
  class TrackCircleDistance{
    
  public:
    /** C'tor takes merge distance */
    TrackCircleDistance(float dCut) : _dCutSquared( dCut*dCut ) , _dCut(dCut){}
    
    /** Merge condition: ... */
    inline bool operator()( nnclu::Element<lcio::Track>* h0, nnclu::Element<lcio::Track>* h1){
      
     

      lcio::Track* trk0 = h0->first ;
      lcio::Track* trk1 = h1->first ;
      
      const TrackInfoStruct* ti0 =  trk0->ext<TrackInfo>() ;
      const TrackInfoStruct* ti1 =  trk1->ext<TrackInfo>() ;


      streamlog_out( DEBUG0 ) << "TrackCircleDistance:: operator() : " <<  trk0->id() << " <-> "  << trk1->id() 
			      << "  (  ti0->zAvg > ti1->zAvg ) " << (  ti0->zAvg > ti1->zAvg )
			      << std::endl ;


      // don't allow  overlaps in z !!!!
      float epsilon = 0. ; 
      if(  ti0->zAvg > ti1->zAvg ){
	
      	if( ti1->zMax > ( ti0->zMin + epsilon )  )
      	  return false ;
	
      } else {
	
      	if( ti0->zMax > ( ti1->zMin + epsilon ) )
      	  return false ;
	
      }
      
      double tl0 = std::abs( trk0->getTanLambda() ) ;
      double tl1 = std::abs( trk1->getTanLambda() ) ;
      
      double dtl = 2. * std::abs( tl0 - tl1 ) / ( tl0 + tl1 ) ;

      streamlog_out( DEBUG ) << "TrackCircleDistance:: operator() : " <<  trk0->id() << " <-> "  << trk1->id() 
			     << "(  dtl > 2.  * _dCut  &&  std::abs( tl0 + tl1 ) > 1.e-2  ) " << (  dtl > 2.  * _dCut  &&  std::abs( tl0 + tl1 ) > 1.e-2  )
			     << std::endl ;
      
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

      // streamlog_out( DEBUG3 ) << "TrackCircleDistance:: operator() : " <<  trk0->id() << " <-> "  << trk1->id() 
      // 			      << " (  std::abs( d0 ) < rIP &&  std::abs( z0 ) < rIP  && std::abs( d1 ) < rIP &&  std::abs( z1 ) < rIP     ) " 
      // 			      << (  std::abs( d0 ) < rIP &&  std::abs( z0 ) < rIP  && std::abs( d1 ) < rIP &&  std::abs( z1 ) < rIP     )
      // 			      << std::endl ;
      

      // // don't merge tracks that come from an area of 20 mm around the IP
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
    
      streamlog_out( DEBUG ) << "TrackCircleDistance:: operator() : " <<  trk0->id() << " <-> "  << trk1->id() 
			     << "( dr < _dCut * std::abs( r0 )  &&  distMS < _dCut * std::abs( r0 )  ) " 
			     << ( dr < _dCut * std::abs( r0 )  &&  distMS < _dCut * std::abs( r0 )  ) 
			     <<  " dr : " << dr 
			     <<  " _dCut * std::abs( r0 ) " << _dCut * std::abs( r0 )
			     << " distMS :" << distMS 
			     << std::endl ;
      
      //      return ( dr < DRMAX && distMS < _dCut * std::abs( r0 )  ) ;
      return ( dr < _dCut * std::abs( r0 )  &&  distMS < _dCut * std::abs( r0 )  ) ;

    }
  
  protected:
    float _dCutSquared ;
    float _dCut ;
  } ; 


  //=======================================================================================
  
  /** Helper class that creates an Elements for an LCOjects of type T.
   */
  template <class T>
  struct MakeLCIOElement{  
    nnclu::Element<T>*  operator()( lcio::LCObject* o) { return new nnclu::Element<T>( (T*) o) ;  }    
  } ;
  
  
  //=======================================================================================

  class Timer{
  public:
    Timer(){
      _clocks.reserve( 100 ) ;
      _names.reserve( 100 ) ;
      
      _clocks.push_back(0) ;
      _names.push_back(  "start"  ) ;
    }
    unsigned registerTimer( const std::string& name ){
      _clocks.push_back(0) ;
      _names.push_back( name ) ;
      return _clocks.size() - 1 ;
    }

    void time(unsigned index){
      _clocks[ index ] = clock() ;
    }
    void start() { time(0) ; }


    std::string toString(){
      
      std::stringstream s ;

      s << " ============= Timer ================================ "  << std::endl ;
      unsigned N=_clocks.size() ;
      for( unsigned i=1 ;  i < N ; ++i){
	s << "    " << _names[i] << " : " <<  double(_clocks[i] - _clocks[i-1] ) / double(CLOCKS_PER_SEC) << std::endl ;
      } 
      s << "         Total  : " <<  double(_clocks[N-1] - _clocks[0] ) / double(CLOCKS_PER_SEC) << std::endl ;
      s << " ==================================================== "  << std::endl ;

      return s.str() ;
    }
  protected:  
    std::vector< clock_t> _clocks ;
    std::vector< std::string > _names ;
  };


}
#endif
