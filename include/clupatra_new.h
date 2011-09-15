#ifndef clupatra_new_h
#define clupatra_new_h

#include <cmath>
#include <math.h>
#include "assert.h"

#include "NNClusterer.h"

#include "lcio.h"
#include "EVENT/TrackerHit.h"
#include "IMPL/TrackImpl.h"
#include "UTIL/Operators.h"

#include "GEAR.h"

#include "MarlinTrk/IMarlinTrack.h"
#include "MarlinTrk/IMarlinTrkSystem.h"


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

  /** Predicate class for 'distance' of NN clustering.
   */
  class HitDistance{
    typedef double PosType ;
  public:

    /** Required typedef for cluster algorithm 
     */
    typedef ClupaHit element_type ;

    /** C'tor takes merge distance */
    HitDistance(float dCut) : _dCutSquared( dCut*dCut ) , _dCut(dCut)  {} 


    /** Merge condition: true if distance  is less than dCut given in the C'tor.*/ 
    inline bool operator()( Hit* h0, Hit* h1){
    
      if( std::abs( h0->Index0 - h1->Index0 ) > 1 ) return false ;
    
      //     //------- don't merge hits from same layer !
      //     if( l0 == l1 )
      //       return false ;

      if( h0->first->layer == h1->first->layer )
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

   
      trk->setdEdx( e/nHit ) ;
      trk->subdetectorHitNumbers().push_back( 1 ) ;  // workaround for bug in lcio::operator<<( Tracks ) - used for picking ....
 
      // FIXME - these are no meaningfull tracks - just a test for clustering tracker hits
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

  struct MarTrk : lcrtrel::LCExtension<MarTrk, MarlinTrk::IMarlinTrack> {} ;

  //------------------------------------------------------------------------------------------

  struct IMarlinTrkFitter{
    
    MarlinTrk::IMarlinTrkSystem* _ts ; 
    
    IMarlinTrkFitter(MarlinTrk::IMarlinTrkSystem* ts) : _ts( ts ) {}
    
    MarlinTrk::IMarlinTrack* operator() (CluTrack* clu) {  
      
      MarlinTrk::IMarlinTrack* trk = _ts->createTrack();
      
      // // store mutual pointers between tracks and 'clusters'
      // trk->setCluster<GCluster>( clu ) ;
      // if( !  clu->ext<ClusterInfo>() )
      //  	clu->ext<ClusterInfo>() = new ClusterInfoStruct ;
      //      clu->ext<ClusterInfo>()->track = trk ;
      // if( clu->size() < 3 ) 
      // 	return trk ;
      //      static HitLayerID tpcLayerID( _ts->indexOfFirstLayer( KalTest::DetID::TPC )  )  ;
      

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
  void addHitsAndFilter( CluTrack* clu, HitListVector& hLV , double dChiMax, double chi2Cut, unsigned maxStep, 
			 bool backward=false) ; 
  //------------------------------------------------------------------------------------------



}
#endif
