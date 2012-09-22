#ifndef clupatra_new_h
#define clupatra_new_h

#include <cmath>
#include <algorithm>
#include <time.h>
#include <math.h>
#include <sstream>
#include <memory>
#include "assert.h"

#include "NNClusterer.h"

#include "lcio.h"
#include "EVENT/TrackerHit.h"
#include "IMPL/TrackImpl.h"
#include "UTIL/Operators.h"
#include "UTIL/CellIDDecoder.h"
#include "UTIL/ILDConf.h"

#include "GEAR.h"

#include "MarlinTrk/IMarlinTrack.h"
#include "MarlinTrk/IMarlinTrkSystem.h"

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"


/** Helper structs that should go to LCIo to make extraction of layer (and subdetector etc easier )
 */
namespace lcio{
  struct ILDDecoder : public CellIDDecoder<TrackerHit>{
    ILDDecoder() :  lcio::CellIDDecoder<TrackerHit>( ILDCellID0::encoder_string ) {} 
  } ;
  static ILDDecoder ILD_cellID ;

  struct ILDTrackTypeBit{
    static const int SEGMENT ;
    static const int COMPOSITE ;
  } ;
} 

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
 
  struct GHit : lcrtrel::LCExtension<GHit, Hit > {} ;
  
  //------------------------------------------------------------------------------------------

  struct MarTrk : lcrtrel::LCExtension<MarTrk, MarlinTrk::IMarlinTrack> {} ;

  //------------------------------------------------------------------------------------------
  
  struct DChi2 : lcrtrel::LCFloatExtension<DChi2> {} ; 

  //----------------------------------------------------------------

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

      streamlog_out( DEBUG ) << "  add hit << "  <<   *first << " to layer " <<  (*first)->first->layer  << std::endl ;

      hLV[ (*first)->first->layer ].push_back( *first )  ;
      ++first ;
    }
  }

  //------------------------------------------------------------------------------------------

  /** Predicate class for 'distance' of NN clustering. */
  class HitDistance{
  public:

    HitDistance(float dCut, float caCut = -1.0 ) : _dCutSquared( dCut*dCut ) , _caCut( caCut ) {} 

    /** Merge condition: true if distance  is less than dCut */ 
    inline bool operator()( Hit* h0, Hit* h1){
    
      if( std::abs( h0->Index0 - h1->Index0 ) > 1 ) return false ;
    
      if( h0->first->layer == h1->first->layer )
	return false ;

      if(  _caCut > 0.  && std::abs( h0->first->layer - h1->first->layer ) == 1 ){

	gear::Vector3D& p0 =  h0->first->pos   ;
	gear::Vector3D& p1 =  h1->first->pos   ;
	
	double cosAlpha = p0.dot( p1 ) / p0.r() / p1.r()  ;

	// merge hits that seem to come from stiff track from the IP
	//fixme: make parameter
	if( cosAlpha > _caCut ) return true ;
      }

      return ( h0->first->pos - h1->first->pos).r2()  < _dCutSquared ;
    }
  protected:
    HitDistance() ;
    float _dCutSquared, _caCut  ;
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
  


  //-----------------------------------------------
  
  struct LCIOTrackConverter{
    
    bool UsePropagate ;
    LCIOTrackConverter() : UsePropagate(false ) {} 

    lcio::Track* operator() (CluTrack* c) ;

  } ;

  //------------------------------------------------------------------------------------------
  /** Predicate class for identifying small clusters. */
  struct ClusterSize { 
    ClusterSize( unsigned n) : _n(n) {} 
    bool operator()(const CluTrack* cl) const {   return cl->size() < _n ;  } 
    unsigned _n ;
  };
  
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
    double _maxChi2Increment ; 
    
    IMarlinTrkFitter(MarlinTrk::IMarlinTrkSystem* ts, double maxChi2Increment=DBL_MAX ) : 
      _ts( ts ) , 
      _maxChi2Increment(maxChi2Increment) {}
    

    MarlinTrk::IMarlinTrack* operator() (CluTrack* clu) ;
  };

  //-------------------------------------------------------------------------------------

  /** Try to add hits from hLV (hit lists per layer) to the cluster. The cluster needs to have a fitted KalTrack associated to it.
   *  Hits are added if the resulting delta Chi2 is less than dChiMax - a maxStep is the maximum number of steps (layers) w/o 
   *  successfully merging a hit.
   */
  int addHitsAndFilter( CluTrack* clu, HitListVector& hLV , double dChiMax, double chi2Cut, unsigned maxStep, ZIndex& zIndex,  bool backward=false, 
			MarlinTrk::IMarlinTrkSystem* trkSys=0) ; 
  //------------------------------------------------------------------------------------------
  
  /** Try to add a hit from the given HitList in layer of subdetector to the track.
   *  A hit is added if the resulting delta Chi2 is less than dChiMax.
   */
  bool addHitAndFilter( int detectorID, int layer, CluTrack* clu, HitListVector& hLV , double dChiMax, double chi2Cut) ; 
  
  //------------------------------------------------------------------------------------------
  /** Split up clusters that have a hit multiplicity of 2,3,4,...,N in at least layersWithMultiplicity. 
   */
  void split_multiplicity( Clusterer::cluster_list& cluList, int layersWithMultiplicity , int N=5) ;

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


  /** Split the cluster into N clusters.
   */
  void create_n_clusters( Clusterer::cluster_type& clu, Clusterer::cluster_list& cluVec ,  unsigned n ) ;



  //=======================================================================================

  struct TrackInfoStruct{  
    TrackInfoStruct() : zMin(0.), zAvg(0.), zMax(0.), startsInner(false), isCentral(false), isForward(false), isCurler(false) {}
    float zMin ;
    float zAvg ;
    float zMax ;
    bool startsInner ;
    bool isCentral   ;
    bool isForward   ;
    bool isCurler    ;
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

  /** helper class for merging split track segments */
  
  class TrackSegmentMerger{
    
  public:
    /** C'tor takes merge distance */
    TrackSegmentMerger(float chi2Max,  MarlinTrk::IMarlinTrkSystem* trksystem, float b) : _chi2Max( chi2Max ) , _trksystem( trksystem), _b(b) {}
    
    float _chi2Max ;
    MarlinTrk::IMarlinTrkSystem* _trksystem ;
    float _b ;
    
    /** Merge condition: ... */
    inline bool operator()( nnclu::Element<lcio::Track>* h0, nnclu::Element<lcio::Track>* h1){
      
      lcio::Track* trk0 = h0->first ;
      lcio::Track* trk1 = h1->first ;


      // protect against merging multiple segments (and thus complete tracks) 
      if(  h0->second || h1->second ) 
	return false ;

      // const TrackInfoStruct* ti0 =  trk0->ext<TrackInfo>() ;
      // const TrackInfoStruct* ti1 =  trk1->ext<TrackInfo>() ;

      // const lcio::TrackState* tsF0 = trk0->getTrackState( lcio::TrackState::AtFirstHit  ) ;
      // const lcio::TrackState* tsL0 = trk0->getTrackState( lcio::TrackState::AtLastHit  ) ;

      // const lcio::TrackState* tsF1 = trk1->getTrackState( lcio::TrackState::AtFirstHit  ) ;
      // const lcio::TrackState* tsL1 = trk1->getTrackState( lcio::TrackState::AtLastHit  ) ;

      // --- get three hits per track : first last midddle

      unsigned nhit0 = trk0->getTrackerHits().size() ;
      unsigned nhit1 = trk1->getTrackerHits().size() ;

      lcio::TrackerHit* thf0 = trk0->getTrackerHits()[ 0 ] ;
      lcio::TrackerHit* thf1 = trk1->getTrackerHits()[ 0 ] ;

      lcio::TrackerHit* thl0 = trk0->getTrackerHits()[ nhit0 - 1 ] ;
      lcio::TrackerHit* thl1 = trk1->getTrackerHits()[ nhit1 - 1 ] ;

      // lcio::TrackerHit* thm1 = trk1->getTrackerHits()[ nhit1 / 2 ] ;
      // lcio::TrackerHit* thm0 = trk0->getTrackerHits()[ nhit0 / 2 ] ;

      int lthf0 = ILD_cellID(  thf0 )[ ILDCellID0::layer ] ;
      int lthf1 = ILD_cellID(  thf1 )[ ILDCellID0::layer ] ;

      int lthl0 = ILD_cellID(  thl0 )[ ILDCellID0::layer ] ;
      int lthl1 = ILD_cellID(  thl1 )[ ILDCellID0::layer ] ;
      
      //      if( lthf0 <= lthl1 && lthf1 <= lthl0 )   return false ; 

      // allow the track segements to overlap slightly  - FIXME: make a parameter ...
      const int overlapRows = 4 ;
      if( lthf0 + overlapRows <= lthl1 && lthf1  + overlapRows  <= lthl0 )   return false ; 

      // now we take the larger segment and see if we can add the three hits from the other segment...

      lcio::Track* trk = ( nhit0 > nhit1 ? trk0 :  trk1 ) ; 
      lcio::Track* oth = ( nhit0 > nhit1 ? trk1 :  trk0 ) ;

      bool  outward = ( nhit0 > nhit1  ?  lthl0 <= lthf1 + overlapRows :  lthl1 <= lthf0 + overlapRows ) ;
      
      unsigned n = oth->getTrackerHits().size() ;
      
      lcio::TrackerHit* th0 =  ( outward ? oth->getTrackerHits()[ 0 ]     :  oth->getTrackerHits()[ n - 1 ] ) ;
      lcio::TrackerHit* th1 =              oth->getTrackerHits()[ n / 2 ] ;
      lcio::TrackerHit* th2 =  ( outward ? oth->getTrackerHits()[ n -1 ] :  oth->getTrackerHits()[ 0 ]     );
      

      // track state at last hit migyt be rubish....
      //      const lcio::TrackState* ts = ( outward ? trk->getTrackState( lcio::TrackState::AtLastHit  ) : trk->getTrackState( lcio::TrackState::AtFirstHit  ) ) ;
      const lcio::TrackState* ts = trk->getTrackState( lcio::TrackState::AtFirstHit  )  ;
      


      streamlog_out( DEBUG3 ) << " *******  TrackSegmentMerger : will extrapolate track " << ( outward ? " outwards\t" : " inwards\t" ) 
			      <<  lcio::lcshort( trk  ) << "     vs:  [" <<   std::hex << oth->id() << std::dec << "]"  << std::endl ;  
      
      // if( trk->id() == 0x0004534d &&  oth->id() == 0x000454a6 ){
      // 	streamlog_out( DEBUG3 )  << " &&&&&&&&&&&&&& Track 1 : \n" << *trk 
      // 				 << " &&&&&&&&&&&&&& Track 2 : \n" << *oth 
      // 				 <<  std::endl ;
      // }

      std::auto_ptr<MarlinTrk::IMarlinTrack> mTrk( _trksystem->createTrack()  ) ;
      
      int nHit = trk->getTrackerHits().size() ;
      
      if( nHit == 0 || ts ==0 )
	return false ;
      
      // float initial_chi2 = trk->getChi2() ;
      // float initial_ndf  = trk->getNdf() ;
      
      streamlog_out( DEBUG3  )  << "               -- extrapolate TrackState : " << lcshort( ts )    << std::endl ;
      
      //need to add a dummy hit to the track
      mTrk->addHit(  trk->getTrackerHits()[0] ) ;  // is this the right hit ??????????
      
      mTrk->initialise( *ts ,  _b ,  MarlinTrk::IMarlinTrack::backward ) ;
      
      double deltaChi ;
      int addHit = 0 ;

      //-----   now try to add the three hits : ----------------
      addHit = mTrk->addAndFit(  th0 , deltaChi, _chi2Max ) ;
      
      streamlog_out( DEBUG3 ) << "    ****  adding first hit : " <<  gear::Vector3D( th0->getPosition() )  
			      << "         added : " << MarlinTrk::errorCode( addHit )
			      << "         deltaChi2: " << deltaChi 
			      << std::endl ;
      
      if( addHit !=  MarlinTrk::IMarlinTrack::success ) return false ;

      //---------------------
      addHit = mTrk->addAndFit(  th1 , deltaChi, _chi2Max ) ;
      
      streamlog_out( DEBUG3 ) << "    ****  adding second hit : " <<  gear::Vector3D( th1->getPosition() )  
			      << "         added : " << MarlinTrk::errorCode( addHit )
			      << "         deltaChi2: " << deltaChi 
			      << std::endl ;
      
      if( addHit !=  MarlinTrk::IMarlinTrack::success ) return false ;

      //--------------------
      addHit = mTrk->addAndFit(  th2 , deltaChi, _chi2Max ) ;
      
      streamlog_out( DEBUG3 ) << "    ****  adding third hit : " <<  gear::Vector3D( th2->getPosition() )  
			      << "         added : " << MarlinTrk::errorCode( addHit )
			      << "         deltaChi2: " << deltaChi 
			      << std::endl ;
      
      if( addHit !=  MarlinTrk::IMarlinTrack::success ) return false ;



      ////////////
      return true ;
      /////////// 
    }

  };
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


      streamlog_out( DEBUG2 ) << "TrackCircleDistance::operator() : " <<  trk0->id() << " <-> "  << trk1->id() 
			      << "  (  ti0->zAvg > ti1->zAvg ) = " << (  ti0->zAvg > ti1->zAvg )
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
      
      if( trk0->getTanLambda() * trk1->getTanLambda()   < 0. ) 
	return false ; // require the same sign


      double dtl = 2. * std::abs( tl0 - tl1 ) / ( tl0 + tl1 ) ;

      streamlog_out( DEBUG2 ) << "TrackCircleDistance::operator() : " <<  trk0->id() << " <-> "  << trk1->id()
			     << " dtl : " << dtl << "   std::abs( tl0 + tl1 ) = " <<  std::abs( tl0 + tl1 ) 
			     << " (  dtl > 2.  * _dCut  &&  std::abs( tl0 + tl1 ) > 1.e-2  ) = " << (  dtl > 2.  * _dCut  &&  std::abs( tl0 + tl1 ) > 1.e-2  )
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

      // double z0 = trk0->getZ0() ;
      // double z1 = trk1->getZ0() ;

      double rIP = 20. ; 

      // streamlog_out( DEBUG2 ) << "TrackCircleDistance::operator() : " <<  trk0->id() << " <-> "  << trk1->id() 
      //  			      << " (  std::abs( d0 ) < rIP &&  std::abs( z0 ) < rIP  && std::abs( d1 ) < rIP &&  std::abs( z1 ) < rIP     ) " 
      //  			      << (  std::abs( d0 ) < rIP &&  std::abs( z0 ) < rIP  && std::abs( d1 ) < rIP &&  std::abs( z1 ) < rIP     )
      //  			      << std::endl ;
      // // // don't merge tracks that come from an area of 20 mm around the IP
      // if(  std::abs( d0 ) < rIP &&  std::abs( z0 ) < rIP  &&
      //  	   std::abs( d1 ) < rIP &&  std::abs( z1 ) < rIP     )
      // 	return false ;

      
      const TrackState* ts0 =  trk0->getTrackState( TrackState::AtFirstHit ) ;
      const TrackState* ts1 =  trk1->getTrackState( TrackState::AtFirstHit ) ;

      double z0 = ( ts0 ? ts0->getReferencePoint()[2]  : 0. ) ;
      double z1 = ( ts1 ? ts1->getReferencePoint()[2]  : 0. ) ;
     

      streamlog_out( DEBUG2 ) << "TrackCircleDistance::operator() : " <<  trk0->id() << " <-> "  << trk1->id() 
			      << " (  std::abs( z0 ) < rIP  &&  std::abs( z1 ) < rIP     ) " 
			      << ( std::abs( z0 ) < rIP  &&  std::abs( z1 ) < rIP     )
			      << std::endl ;
      
      if(  std::abs( z0 ) < rIP  && std::abs( z1 ) < rIP     )
      	return false ;



      double p0 = trk0->getPhi() ;
      double p1 = trk1->getPhi() ;

      double x0 = ( r0 - d0 ) * sin( p0 ) ;
      double x1 = ( r1 - d1 ) * sin( p1 ) ;

      double y0 = ( d0 - r0 ) * cos( p0 ) ;
      double y1 = ( d1 - r1 ) * cos( p1 ) ;
    
      double dr = 2. * std::abs( r0abs - r1abs )  / (r0abs + r1abs )  ;

      double distMS = sqrt ( ( x0 - x1 ) * ( x0 - x1 ) + ( y0 - y1 ) * ( y0 - y1 )  ) ;
    
      streamlog_out( DEBUG2 ) << "TrackCircleDistance:: operator() : " <<  trk0->id() << " <-> "  << trk1->id() 
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
  
  struct TrackZSort {  // sort tracks wtr to abs(z_average )  
    inline bool operator()( lcio::Track* l, lcio::Track* r) {      
      return (  std::abs( l->ext<TrackInfo>()->zAvg )   <   std::abs( r ->ext<TrackInfo>()->zAvg )  ) ; 
    }
  };
  
  
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
