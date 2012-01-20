#include "ClupatraProcessor.h"

#include <time.h>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

#include <math.h>
#include <cmath>

//---- MarlinUtil 
#include "MarlinCED.h"

//---- LCIO ---
#include "IMPL/LCCollectionVec.h"
#include "IMPL/TrackImpl.h"
#include "IMPL/TrackerHitImpl.h"
#include "EVENT/SimTrackerHit.h"
#include "IMPL/LCFlagImpl.h"
#include "UTIL/Operators.h"
#include "UTIL/LCTOOLS.h"

#include "UTIL/CellIDDecoder.h"

#include "LCIterator.h"


//---- GEAR ----
#include "marlin/Global.h"
#include "gear/GEAR.h"
#include "gear/TPCParameters.h"
#include "gear/ZPlanarParameters.h"
#include "gear/ZPlanarLayerLayout.h"
#include "gear/PadRowLayout2D.h"
#include "gear/BField.h"


#include "MarlinTrk/Factory.h"
#include "MarlinTrk/IMarlinTrack.h"
#include "MarlinTrk/IMarlinTrkSystem.h"


using namespace lcio ;
using namespace marlin ;


#include <float.h>

#include "clupatra_new.h"
using namespace clupatra_new ;



/** helper method to create a track collections and add it to the event */
inline LCCollectionVec* newTrkCol(const std::string& name, LCEvent * evt ){

  LCCollectionVec* col = new LCCollectionVec( LCIO::TRACK ) ;  

  LCFlagImpl hitFlag(0) ;
  hitFlag.setBit( LCIO::TRBIT_HITS ) ;
  col->setFlag( hitFlag.getFlag()  ) ;

  evt->addCollection( col , name ) ;

  return col ;
}
//----------------------------------------------------------------
/** helper method to get the collection from the event */
inline LCCollection* getCollection(  const std::string& name, LCEvent * evt ){

  LCCollection* col = 0 ;
  try{   col = evt->getCollection( name )  ; 
  } catch( lcio::DataNotAvailableException& e) { 
    streamlog_out( DEBUG4 ) <<  " input collection not in event : " << name << "  !!!  " << std::endl ;  
  } 
  return col ;
}

//----------------------------------------------------------------

template <class In, class Pred> In find_smallest(In first, In last, Pred p, double& d){
  
  In res =  last ;

  double min = DBL_MAX ;
 
  while( first != last ){

    double val = p( *first) ;

    if(  val < min ){

      res = first ;
      min = val ;
    }
 
    ++first ;
  }

  d = min ;

  streamlog_out( DEBUG ) << " found smallest " << *res << " min : " << d << std::endl ;

  return res ;
}
//----------------------------------------------------------------
struct Distance3D2{
  gear::Vector3D _pos ;
  Distance3D2( const gear::Vector3D& pos) : _pos( pos ) {}
  template <class T>
  double operator()( const T* t) { 
    gear::Vector3D p( t->getPosition() ) ;
    return ( p - _pos ).r2() ; 

  }
};
//----------------------------------------------------------------
struct MeanAbsZOfTrack{
  double operator()( const Track* t){
    double z = 0 ;
    int hitCount = 0 ;
    const TrackerHitVec& hV = t->getTrackerHits() ;
    for(unsigned i=0, N = hV.size() ; i<N ; ++i){
      z += hV[i]->getPosition()[2]  ;
      ++hitCount ;
    }
    return std::abs(  z )  / hitCount ;
  }
};

//----------------------------------------------------------------


ClupatraProcessor aClupatraProcessor ;


ClupatraProcessor::ClupatraProcessor() : Processor("ClupatraProcessor") {
  
  // modify processor description
  _description = "ClupatraProcessor : nearest neighbour clustering seeded pattern recognition" ;
  
  
  
  registerInputCollection( LCIO::TRACKERHIT,
			   "TPCHitCollection" , 
			   "Name of the tpc hit input collections"  ,
			   _colName ,
			   "AllTPCTrackerHits"  ) ;
  
  
  registerOutputCollection( LCIO::TRACK,
			    "OutputCollection" , 
			    "Name of the output collection"  ,
			    _outColName ,
			    std::string("ClupatraTracks" ) ) ;
  
  
  registerProcessorParameter( "DistanceCut" , 
			      "Cut for distance between hits in mm"  ,
			      _distCut ,
			      (float) 40.0 ) ;
  
  
  registerProcessorParameter( "MinimumClusterSize" , 
			      "minimum number of hits per cluster"  ,
			      _minCluSize ,
			      (int) 6) ;
  
  
  registerProcessorParameter( "DuplicatePadRowFraction" , 
			      "allowed fraction of hits in same pad row per track"  ,
			      _duplicatePadRowFraction,
			      (float) 0.1 ) ;
  
  // registerProcessorParameter( "RCut" , 
  // 			      "Cut for r_min in mm"  ,
  // 			      _rCut ,
  // 			      (float) 0.0 ) ;
  
  registerProcessorParameter( "MaxDeltaChi2" , 
 			      "the maximum delta Chi2 for which a hit is added to a track segement"  ,
 			      _dChi2Max ,
 			      (float) 35. ) ;

  registerProcessorParameter( "Chi2Cut" , 
 			      "the maximum chi2-distance for which a hit is considered for merging "  ,
 			      _chi2Cut ,
 			      (float) 100. ) ;
  
  registerProcessorParameter( "MaxStepWithoutHit" , 
 			      "the maximum number of layers without finding a hit before hit search search is stopped "  ,
 			      _maxStep ,
 			      (int) 3 ) ;

  registerProcessorParameter( "PadRowRange" , 
			      "number of pad rows used in initial seed clustering"  ,
			      _padRowRange ,
			      (int) 12) ;
 
  registerProcessorParameter( "NumberOfZBins" , 
			      "number of bins in z over total length of TPC - hits from different z bins are nver merged"  ,
			      _nZBins,
			      (int) 150 ) ;


  registerProcessorParameter( "MinLayerFractionWithMultiplicity" , 
			      "minimum fraction of layers that have a given multiplicity, when forcing a cluster into sub clusters"  ,
			      _minLayerFractionWithMultiplicity,
			      (float) 0.5 ) ;

  registerProcessorParameter( "MinLayerNumberWithMultiplicity" , 
			      "minimum number of layers that have a given multiplicity, when forcing a cluster into sub clusters"  ,
			      _minLayerNumberWithMultiplicity,
			      (int) 3 ) ;

  registerProcessorParameter("MultipleScatteringOn",
			     "Use MultipleScattering in Fit",
			     _MSOn,
			     bool(true));
  
  registerProcessorParameter("EnergyLossOn",
			     "Use Energy Loss in Fit",
			     _ElossOn,
			     bool(true));
  
  registerProcessorParameter("SmoothOn",
			     "Smooth All Mesurement Sites in Fit",
			     _SmoothOn,
			     bool(false));

  registerProcessorParameter("pickUpSiHits",
			     "try to pick up hits from Si-trackers",
			     _pickUpSiHits,
			     bool(false));
  

  registerOptionalParameter( "SITHitCollection" , 
			     "name of the SIT hit collections - used to extend TPC tracks if (pickUpSiHits==true)"  ,
			     _sitColName ,
			     std::string("SITTrackerHits")  ) ;

  registerOptionalParameter( "VXDHitCollection" , 
			     "name of the VXD hit collections - used to extend TPC tracks if (pickUpSiHits==true)"  ,
			     _vxdColName ,
			     std::string("VTXTrackerHits")  ) ;
  
  registerProcessorParameter("CreateDebugCollections",
			     "optionally create some debug collection with intermediate track segments and used and unused hits",
			     _createDebugCollections,
			     bool(false));
  
}


void ClupatraProcessor::init() { 

  // usually a good idea to
  printParameters() ;
  
  _trksystem =  MarlinTrk::Factory::createMarlinTrkSystem( "KalTest" , marlin::Global::GEAR , "" ) ;
  
  _trksystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useQMS,        _MSOn ) ;
  _trksystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::usedEdx,       _ElossOn) ;
  _trksystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing,  _SmoothOn) ;
  _trksystem->init() ;  
  
  _nRun = 0 ;
  _nEvt = 0 ;
  
  // CEDPickingHandler::getInstance().registerFunction( LCIO::TRACKERHIT , &printTrackerHit ) ; 
  // CEDPickingHandler::getInstance().registerFunction( LCIO::TRACK , &printTrackShort ) ; 
  // CEDPickingHandler::getInstance().registerFunction( LCIO::SIMTRACKERHIT , &printSimTrackerHit ) ; 
  
}

void ClupatraProcessor::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 


void ClupatraProcessor::processEvent( LCEvent * evt ) { 
  
  //  clock_t start =  clock() ; 
  Timer timer ;
  unsigned t_init       = timer.registerTimer(" initialization      " ) ;
  unsigned t_seedtracks = timer.registerTimer(" extend seed tracks  " ) ;
  unsigned t_recluster  = timer.registerTimer(" recluster leftovers " ) ;
  unsigned t_split      = timer.registerTimer(" split clusters      " ) ;
  unsigned t_finalfit   = timer.registerTimer(" final refit         " ) ;
  unsigned t_merge      = timer.registerTimer(" merge segments      " ) ;
  unsigned t_pickup     = timer.registerTimer(" pick up Si hits     " ) ;
  
  timer.start() ;
  
  // the clupa wrapper hits that hold pointers to LCIO hits plus some additional parameters
  // create them in a vector for convenient memeory mgmt 
  std::vector<ClupaHit> clupaHits ;
  
  // on top of the clupahits we need the tiny wrappers for clustering - they are created on the heap
  // and we put them in a vector of pointers that takes ownership (i.e. deletes them at the end)
  HitVec nncluHits ;        
  nncluHits.setOwner( true ) ;


  // this is the final list of cluster tracks
  Clusterer::cluster_list cluList ;    
  cluList.setOwner() ;
  

  HitDistance dist( _distCut ) ;
  
  LCIOTrackConverter converter ;
  
  const gear::TPCParameters& gearTPC = Global::GEAR->getTPCParameters() ;
  const gear::PadRowLayout2D& padLayout = gearTPC.getPadLayout() ;

  unsigned maxTPCLayers =  padLayout.getNRows() ;
  
  double driftLength = gearTPC.getMaxDriftLength() ;
  ZIndex zIndex( -driftLength , driftLength , _nZBins  ) ; 
  

  LCCollection* col = 0 ;

  try{   col =  dynamic_cast<LCCollection*> (evt->getCollection( _colName )  ); 
    
  } catch( lcio::DataNotAvailableException& e) { 
    
    streamlog_out( WARNING ) <<  " input collection not in event : " << _colName << "   - nothing to do  !!! " << std::endl ;  
    
    return ;
  } 
      
  //===============================================================================================
  //   create clupa and clustering hits for every lcio hit
  //===============================================================================================
  
  CellIDDecoder<TrackerHit> idDec( col ) ;
  
  int nHit = col->getNumberOfElements() ;
  
  clupaHits.resize( nHit ) ;       // creates clupa hits (w/ default c'tor)
  nncluHits.reserve( nHit ) ;


  for(int i=0 ; i < nHit ; ++i ) {
    
    //------
    
    TrackerHit* th = (TrackerHit*) col->getElementAt(i) ;
    
    ClupaHit* ch  = & clupaHits[i] ; 
    
    Hit* gh =  new Hit( ch ) ;
    
    nncluHits.push_back( gh ) ;
    
    //-------
    // th->ext<HitInfo>() = ch ;  // assign the clupa hit to the LCIO hit for memory mgmt
    
    ch->lcioHit = th ; 
    
    ch->pos = gear::Vector3D(  th->getPosition() ) ;
    
    //  int padIndex = padLayout.getNearestPad( ch->pos.rho() , ch->pos.phi() ) ;
    //    ch->layer = padLayout.getRowNumber( padIndex ) ;
    ch->layer = idDec( th )[ ILDCellID0::layer ] ;

    ch->zIndex = zIndex( th ) ;
    
    //ch->phiIndex = ....
    
  } 

  //--------------------------------------------------------------------------------------------------------- 
  
  std::sort( nncluHits.begin(), nncluHits.end() , ZSort() ) ;
  
  //--------------------------------------------------------------------------------------------------------- 
  
  HitListVector hitsInLayer( maxTPCLayers ) ;
  addToHitListVector(  nncluHits.begin(), nncluHits.end() , hitsInLayer  ) ;
  
  //---------------------------------------------------------------------------------------------------------

  //===============================================================================================
  //   create output collections  ( some optional )
  //===============================================================================================

  const bool writeSeedCluster      = _createDebugCollections ;
  const bool writeCluTrackSegments = _createDebugCollections ;
  const bool writeLeftoverClusters = _createDebugCollections ;
  
  LCCollectionVec* seedCol =  ( writeSeedCluster      ?  newTrkCol( "ClupatraSeedCluster"      , evt )  :   0   )  ; 
  LCCollectionVec* cluCol  =  ( writeCluTrackSegments ?  newTrkCol( "ClupatraInitialTrackSegments" , evt )  :   0   )  ; 
  LCCollectionVec* locCol  =  ( writeCluTrackSegments ?  newTrkCol( "ClupatraLeftoverClusters" , evt )  :   0   )  ; 
  

  LCCollectionVec* tsCol  =  newTrkCol( "ClupatraTrackSegments" , evt ) ;
  
  LCCollectionVec* outCol =  newTrkCol( _outColName  , evt )  ; 

  //---------------------------------------------------------------------------------------------------------
  
  timer.time(t_init ) ; 

  //===============================================================================================
  // first main step of clupatra:
  //   * cluster in pad row range - starting from the outside - to find clean cluster segments
  //   * extend the track segments with best matching hits, based on extrapolation to next layer(s)
  //   * add the hits and apply a Kalman filter step ( track segement is always best estimate )
  //   * repeat in backward direction ( after smoothing back, to get a reasonable track 
  //     state for extrapolating backwards )
  //===============================================================================================

  Clusterer nncl ;
  
  int outerRow = maxTPCLayers - 1 ;

  nnclu::PtrVector<MarlinTrk::IMarlinTrack> seedTrks ;
  seedTrks.setOwner() ; // memory mgmt - will delete MarlinTrks at the end
  
  IMarlinTrkFitter fitter( _trksystem ) ;

  
  while( outerRow > _padRowRange ) {
    
    HitVec hits ;
    hits.reserve( nHit ) ;
    
    // add all hits in pad row range to hits
    for(int iRow = outerRow ; iRow > ( outerRow - _padRowRange) ; --iRow ) {
      
      std::copy( hitsInLayer[ iRow ].begin() , hitsInLayer[ iRow ].end() , std::back_inserter( hits )  ) ;
    }
    
    //-----  cluster in given pad row range  -----------------------------
    Clusterer::cluster_list sclu ;    
    sclu.setOwner() ;  
    
    nncl.cluster_sorted( hits.begin(), hits.end() , std::back_inserter( sclu ), dist , _minCluSize ) ;
    
    if( writeSeedCluster ) {
      std::transform( sclu.begin(), sclu.end(), std::back_inserter( *seedCol ) , converter ) ;
    }
    
    // remove clusters whith too many duplicate hits per pad row
    Clusterer::cluster_list bclu ;    // bad clusters  
    bclu.setOwner() ;      
    
    split_list( sclu, std::back_inserter(bclu),  DuplicatePadRows( maxTPCLayers, _duplicatePadRowFraction  ) ) ;
    
    // free hits from bad clusters 
    std::for_each( bclu.begin(), bclu.end(), std::mem_fun( &CluTrack::freeElements ) ) ;
    
    
    // ---- now we also need to remove the hits from good cluster seeds from the hitsInLayers:
    for( Clusterer::cluster_list::iterator sci=sclu.begin(), end= sclu.end() ; sci!=end; ++sci ){
      for( Clusterer::cluster_type::iterator ci=(*sci)->begin(), end= (*sci)->end() ; ci!=end;++ci ){
	
	// this is not cheap ...
	hitsInLayer[ (*ci)->first->layer ].remove( *ci )  ; 
      }
    }
    
    
    
    
    std::transform( sclu.begin(), sclu.end(), std::back_inserter( seedTrks) , fitter ) ;
    
    
    for( Clusterer::cluster_list::iterator icv = sclu.begin(), end =sclu.end()  ; icv != end ; ++ icv ) {
      

     addHitsAndFilter( *icv , hitsInLayer , _dChi2Max, _chi2Cut , _maxStep , zIndex ) ; 
      
      static const bool backward = true ;
      addHitsAndFilter( *icv , hitsInLayer , _dChi2Max, _chi2Cut , _maxStep , zIndex, backward ) ; 
    } 

    // merge the good clusters to final list
    cluList.merge( sclu ) ;

    outerRow -= _padRowRange ;
    
  } //while outerRow > padRowRange 
  
 
  //---------------------------------------------------------------------------------------------------------

  //  ---- store track segments from the first main step  -----  
  if( writeCluTrackSegments )
    std::transform( cluList.begin(), cluList.end(), std::back_inserter( *cluCol ) , converter ) ;

  //---------------------------------------------------------------------------------------------------------
  
  
  timer.time( t_seedtracks ) ;
  
  timer.time( t_recluster ) ;
  
  //===============================================================================================
  //  do a global reclustering in all leftover hits
  //===============================================================================================
  
  outerRow = maxTPCLayers - 1 ;
  
  int padRangeRecluster = 50 ; // FIXME: make parameter 
  
  while( outerRow > 0 ) {
    

    Clusterer::cluster_list loclu ; // leftover clusters
    loclu.setOwner() ;
    
    HitVec hits ;


    int  minRow = ( ( outerRow - padRangeRecluster ) > -1 ?  ( outerRow - padRangeRecluster ) : -1 ) ;

    streamlog_out( DEBUG4 ) << "   recluster in the range : " << outerRow << " - " <<  minRow << std::endl ;

    // add all hits in pad row range to hits
    for(int iRow = outerRow ; iRow > minRow ; --iRow ) {
      
      std::copy( hitsInLayer[ iRow ].begin() , hitsInLayer[ iRow ].end() , std::back_inserter( hits )  ) ;
    }
    
    
    
    //  CluTrack debug ;
    
    // int startRowRecluster = 100 ; // FIXME: make parameter 
    // for(unsigned iRow = startRowRecluster ; iRow <  maxTPCLayers ; ++iRow ) {
    
    //   std::copy( hitsInLayer[ iRow ].begin() , hitsInLayer[ iRow ].end() , std::back_inserter( hits )  ) ;
    // }
    
    HitDistance distSmall( _distCut ) ; //   / 2. ) ;
    nncl.cluster( hits.begin(), hits.end() , std::back_inserter( loclu ),  distSmall , _minCluSize ) ;
    
    
    if( writeLeftoverClusters )
    std::transform( loclu.begin(), loclu.end(), std::back_inserter( *locCol ) , converter ) ;
    
    
    
    // timer.time( t_recluster ) ;
    
    //===============================================================================================
    //  now we split the clusters based on their hit multiplicities
    //===============================================================================================
    
    
    _dChi2Max = 5. * _dChi2Max ; //FIXME !!!!!!!!!

    for( Clusterer::cluster_list::iterator it= loclu.begin(), end= loclu.end() ; it != end ; ++it ){
      
      CluTrack* clu = *it ;
      
      streamlog_out(  DEBUG4 ) << " **** left over cluster with size : " << clu->size() << std::endl ;
      
      std::vector<int> mult(8) ; 
      // get hit multiplicities up to 6 ( 7 means 7 or higher ) 
      
      getHitMultiplicities( clu , mult ) ;
      
      streamlog_out(  DEBUG4 ) << " **** left over cluster with hit multiplicities: \n" ;
      for( unsigned i=0,n=mult.size() ; i<n ; ++i) {
	streamlog_out(  DEBUG4 ) << "     m["<<i<<"] = " <<  mult[i] << "\n"  ;
      }
      
      
      if( float( mult[3]) / mult[0]  >= _minLayerFractionWithMultiplicity &&  mult[3] >  _minLayerNumberWithMultiplicity ) {
	
	Clusterer::cluster_list reclu ; // reclustered leftover clusters
	reclu.setOwner() ;
	
	create_three_clusters( *clu , reclu ) ;
	
	std::transform( reclu.begin(), reclu.end(), std::back_inserter( seedTrks) , fitter ) ;
	
	for( Clusterer::cluster_list::iterator ir= reclu.begin(), end= reclu.end() ; ir != end ; ++ir ){
	  
	  streamlog_out( DEBUG4 ) << " extending triplet clustre  of length " << (*ir)->size() << std::endl ;
	  
	  addHitsAndFilter( *ir , hitsInLayer , _dChi2Max, _chi2Cut , _maxStep , zIndex) ; 
	  static const bool backward = true ;
	  addHitsAndFilter( *ir , hitsInLayer , _dChi2Max, _chi2Cut , _maxStep , zIndex, backward ) ; 
	}
	
	cluList.merge( reclu ) ;
      } 
      
      else if( float( mult[2]) / mult[0]  >= _minLayerFractionWithMultiplicity &&  mult[2] >  _minLayerNumberWithMultiplicity ) {
	
	Clusterer::cluster_list reclu ; // reclustered leftover clusters
	reclu.setOwner() ;
	
	create_two_clusters( *clu , reclu ) ;
	
	std::transform( reclu.begin(), reclu.end(), std::back_inserter( seedTrks) , fitter ) ;
	
	for( Clusterer::cluster_list::iterator ir= reclu.begin(), end= reclu.end() ; ir != end ; ++ir ){
	  
	  streamlog_out( DEBUG4 ) << " extending doublet clustre  of length " << (*ir)->size() << std::endl ;
	  
	  addHitsAndFilter( *ir , hitsInLayer , _dChi2Max, _chi2Cut , _maxStep , zIndex) ; 
	  static const bool backward = true ;
	  addHitsAndFilter( *ir , hitsInLayer , _dChi2Max, _chi2Cut , _maxStep , zIndex, backward ) ; 
	} 
	
	cluList.merge( reclu ) ;
	
      }
      else if( float( mult[1]) / mult[0]  >= _minLayerFractionWithMultiplicity &&  mult[1] >  _minLayerNumberWithMultiplicity ) {    
	
	
	seedTrks.push_back( fitter( *it )  );
	
	cluList.push_back( *it ) ;
	
	it = loclu.erase( it ) ;
	--it ; // erase returns iterator to next element 
	
      } else {
	
	//  discard cluster and free hits
	clu->freeElements() ; 
      }
      
    }

  
    outerRow -=  padRangeRecluster ; 

  }
  
  //===============================================================================================
  //  try again to gobble up hits at the ends ....
  //===============================================================================================
  
  for( Clusterer::cluster_list::iterator icv = cluList.begin() ; icv != cluList.end() ; ++ icv ) {
    addHitsAndFilter( *icv , hitsInLayer , _dChi2Max, _chi2Cut , _maxStep , zIndex) ; 
    static const bool backward = true ;
    addHitsAndFilter( *icv , hitsInLayer , _dChi2Max, _chi2Cut , _maxStep , zIndex, backward ) ; 
  }
  
  
  timer.time( t_split ) ;

  //===============================================================================================
  //  now refit the tracks 
  //===============================================================================================

  streamlog_out( DEBUG2 ) << " ===========    refitting final " << cluList.size() << " track segments  "   << std::endl ;

  nnclu::PtrVector<MarlinTrk::IMarlinTrack> finalTrks ;
  finalTrks.setOwner() ;  // memory mgmt - will delete MarlinTrks at the end
  finalTrks.reserve( cluList.size() ) ;
  
  std::transform( cluList.begin(), cluList.end(), std::back_inserter( finalTrks) , IMarlinTrkFitter(_trksystem)  ) ;
  
  std::for_each( finalTrks.begin() , finalTrks.end() , std::mem_fun_t< int, MarlinTrk::IMarlinTrack >( &MarlinTrk::IMarlinTrack::smooth )  ) ;


  std::transform( cluList.begin(), cluList.end(), std::back_inserter( *tsCol ) , converter ) ;

  
  timer.time( t_finalfit) ;
  
  //===============================================================================================
  //   optionally create collections of used and unused TPC hits 
  //===============================================================================================
  
  if( _createDebugCollections ) {
    LCCollectionVec* usedHits   = new LCCollectionVec( LCIO::TRACKERHIT ) ;   ;
    LCCollectionVec* unUsedHits = new LCCollectionVec( LCIO::TRACKERHIT ) ;   ;
    evt->addCollection( usedHits ,   "ClupatraUsedTPCHits"   ) ;
    evt->addCollection( unUsedHits , "ClupatraUnUsedTPCHits" ) ;
    usedHits->setSubset() ;
    unUsedHits->setSubset() ;
    usedHits->reserve(   nncluHits.size() ) ;
    unUsedHits->reserve( nncluHits.size() ) ;
    
    for( HitVec::iterator it = nncluHits.begin(), end = nncluHits.end(); it!=end;++it ){
      
      if( (*it)->second != 0 ){   usedHits->push_back( (*it)->first->lcioHit ) ;
      } else {                  unUsedHits->push_back( (*it)->first->lcioHit ) ;          
      }
    }
  }

  //===============================================================================================
  //  merge track segements 
  //===============================================================================================
  
  
  static const int merge_track_segments = true ;

  if( merge_track_segments ) {

    
    typedef nnclu::NNClusterer<Track> TrackClusterer ;
    TrackClusterer nntrkclu ;

    TrackClusterer::element_vector trkVec ;
    trkVec.setOwner() ;
    TrackClusterer::cluster_vector trkCluVec ;
    trkCluVec.setOwner() ;

    trkVec.reserve( tsCol->size()  ) ;
    trkCluVec.reserve(  tsCol->size() ) ;

    std::for_each( tsCol->begin() , tsCol->end() , ComputeTrackerInfo()  ) ;
    
    std::transform( tsCol->begin() , tsCol->end() , std::back_inserter( trkVec ) , MakeLCIOElement<Track>() ) ; 

    TrackCircleDistance trkMerge( 0.1 ) ; 

    nntrkclu.cluster( trkVec.begin() , trkVec.end() , std::back_inserter( trkCluVec ), trkMerge , 2  ) ;


    streamlog_out( DEBUG4 ) << " ===== merged tracks - # cluster: " << trkCluVec.size()   
			    << " from " << tsCol->size() << " track segments "    << "  ============================== " << std::endl ;
    
    for(  TrackClusterer::cluster_vector::iterator it= trkCluVec.begin() ; it != trkCluVec.end() ; ++it) {
      
      streamlog_out( DEBUG4 ) <<  lcio::header<Track>() << std::endl ;
      
      TrackClusterer::cluster_type*  trkClu = *it ;

      std::list<Track*> mergedTrk ;

      for( TrackClusterer::cluster_type::iterator itC = trkClu->begin() ; itC != trkClu->end() ; ++ itC ){
	
	streamlog_out( DEBUG4 ) << lcshort(  (*itC)->first ) << std::endl ; 
	
	mergedTrk.push_back( (*itC)->first ) ; 
      }
      
      //      std::sort(  mergedTrk.begin() , mergedTrk.end() , TrackZSort() ) ;
      mergedTrk.sort( TrackZSort() ) ;
      
      // ====== create a new LCIO track for the merged cluster ...
      TrackImpl* trk = new TrackImpl ;
      
      // == and copy all the hits 
      unsigned hitCount = 0 ;
      for( std::list<Track*>::iterator itML = mergedTrk.begin() ; itML != mergedTrk.end() ; ++ itML ){
	
      	const TrackerHitVec& hV = (*itML)->getTrackerHits() ;
	for(unsigned i=0, n=hV.size() ; i<n ; ++i){
	  
      	  trk->addHit( hV[i] ) ;
      	}
	hitCount  += hV.size()  ;
      }
      
      // take track states from first and last track :
      Track* firstTrk = mergedTrk.front() ;
      Track* lastTrk  = mergedTrk.back() ;
      
      const TrackState* ts = 0 ; 
      ts = firstTrk->getTrackState( lcio::TrackState::AtIP  ) ;
      if( ts ) trk->addTrackState( new TrackStateImpl( ts )  ) ;

      ts = firstTrk->getTrackState( lcio::TrackState::AtFirstHit  ) ;
      if( ts ) trk->addTrackState( new TrackStateImpl( ts )  ) ;

      ts = lastTrk->getTrackState( lcio::TrackState::AtLastHit  ) ;
      if( ts ) trk->addTrackState( new TrackStateImpl( ts )  ) ;

      ts = lastTrk->getTrackState( lcio::TrackState::AtCalorimeter  ) ;
      if( ts ) trk->addTrackState( new TrackStateImpl( ts )  ) ;


      trk->ext<MarTrk>() = firstTrk->ext<MarTrk>() ;
      
      int hitsInFit = firstTrk->getTrackerHits().size() ;

      trk->subdetectorHitNumbers().resize( 2 * ILDDetID::ETD ) ;
      trk->subdetectorHitNumbers()[ 2 * ILDDetID::TPC - 1 ] =  hitsInFit ;  
      trk->subdetectorHitNumbers()[ 2 * ILDDetID::TPC - 2 ] =  hitCount ;  
      

      streamlog_out( DEBUG2 ) << "   create new merged track from bestTrack parameters - ptr to MarlinTrk : " << trk->ext<MarTrk>()  
			      << "   with subdetector hit numbers  : " <<  trk->subdetectorHitNumbers()[0 ] << " , " <<  trk->subdetectorHitNumbers()[1] 
			      << std::endl ;
      

      outCol->addElement( trk )  ;
    }

    // for(  TrackClusterer::cluster_vector::iterator it= trkCluVec.begin() ; it != trkCluVec.end() ; ++it) {
      
    //   streamlog_out( DEBUG2 ) <<  lcio::header<Track>() << std::endl ;
      
    //   TrackClusterer::cluster_type*  trkClu = *it ;

    //   std::list<Track*> mergedTrk ;
    //   for( TrackClusterer::cluster_type::iterator itC = trkClu->begin() ; itC != trkClu->end() ; ++ itC ){
	
    // 	streamlog_out( DEBUG2 ) << lcshort(  (*itC)->first ) << std::endl ; 

    // 	mergedTrk.push_back( (*itC)->first ) ; 
    //   }

    //   // create a new LCIO track for the merged cluster ...
    //   TrackImpl* trk = new TrackImpl ;
    //   Track* bestTrk = 0 ;
      
    //   double min = 1.e99 ;
    //   int hitCount = 0 ;
      
    //   for( std::list<Track*>::iterator itML = mergedTrk.begin() ; itML != mergedTrk.end() ; ++ itML ){
	
    // 	double z = 0. ;
    //   	const TrackerHitVec& hV = (*itML)->getTrackerHits() ;
	
    //   	for(unsigned i=0 ; i < hV.size() ; ++i){
	  
    //   	  trk->addHit( hV[i] ) ;

    //   	  ++hitCount ;

    // 	  z+= hV[i]->getPosition()[2] ;
	  
    //   	}

    // 	double val = std::abs( z) / hV.size() ;  // take parameters from track segment that is closests to IP .....

    // 	//double val = (*itML)->getChi2() / (*itML)->getNdf() ; 
	
    // 	streamlog_out( DEBUG ) << "   ---  " <<   &mergedTrk <<   " new  value : " << val << "   -- " << lcshort( *itML ) << std::endl ;

    //   	if( val < min ){
    //   	  bestTrk = (*itML) ;
    //   	  min  = val  ;
    //   	}
    //   }
      
      
    //   if( bestTrk != 0 ){ 
	
    // 	int hitsInFit = bestTrk->getTrackerHits().size() ;
	
    // 	trk->setD0( bestTrk->getD0() ) ;
    // 	trk->setOmega( bestTrk->getOmega() ) ;
    // 	trk->setPhi( bestTrk->getPhi() ) ;
    // 	trk->setZ0( bestTrk->getZ0() ) ;
    // 	trk->setTanLambda( bestTrk->getTanLambda() ) ;
    // 	trk->setCovMatrix( bestTrk->getCovMatrix()  ) ;
    // 	// ...
    // 	trk->ext<MarTrk>() = bestTrk->ext<MarTrk>() ;
	

    // 	trk->subdetectorHitNumbers().resize( 2 * ILDDetID::ETD ) ;
    // 	trk->subdetectorHitNumbers()[ 2 * ILDDetID::TPC - 1 ] =  hitsInFit ;  
    // 	trk->subdetectorHitNumbers()[ 2 * ILDDetID::TPC - 2 ] =  hitCount ;  
	

    // 	streamlog_out( DEBUG2 ) << "   create new merged track from bestTrack parameters - ptr to MarlinTrk : " << trk->ext<MarTrk>()  
    // 				<< "  with subdetector hit numbers  : " <<  trk->subdetectorHitNumbers()[0 ] << " , " <<  trk->subdetectorHitNumbers()[1] 
    // 				<< std::endl ;
	

    // 	outCol->addElement( trk )  ;
    //   }
    //   else{
    // 	streamlog_out( ERROR ) << "   no best track found for merged tracks ... !? " << std::endl ; 
    // 	delete trk ;
    //   }
      
      
    // }

    //---------------------------------------------------------------------------------------------
    // add all tracks that have not been merged :

    for( TrackClusterer::element_vector::iterator it = trkVec.begin(); it != trkVec.end() ;++it){

      if( (*it)->second == 0 ){
	
	TrackImpl* trk = dynamic_cast<TrackImpl*>( (*it)->first ) ;
	
	TrackImpl* t =   new TrackImpl( *trk ) ;
	
	t->ext<TrackInfo>() = 0 ; // set extension to 0 to prevent double free ... 

	t->ext<MarTrk>() = dynamic_cast<TrackImpl*>( (*it)->first )->ext<MarTrk>() ;


	// unsigned nHit = trk->getTrackerHits().size() ;
	// t->subdetectorHitNumbers().resize( 2 * ILDDetID::ETD ) ;
	// t->subdetectorHitNumbers()[ 2 * ILDDetID::TPC - 1 ] =  nHit ;  
	// t->subdetectorHitNumbers()[ 2 * ILDDetID::TPC - 2 ] =  nHit ;  

	streamlog_out( DEBUG2 ) << "   create new track from existing LCIO track  - ptr to MarlinTrk : " << t->ext<MarTrk>()  << std::endl ;
	
	outCol->addElement( t ) ;
      }
    }
    
  }
  timer.time( t_merge ) ;  

  //---------------------------------------------------------------------------------------------------------
  //    pick up hits from Si trackers
  //---------------------------------------------------------------------------------------------------------

  if( _pickUpSiHits ){
    
    pickUpSiTrackerHits( outCol , evt ) ;
    
  }
  //---------------------------------------------------------------------------------------------------------

  timer.time( t_pickup ) ;  

  
  streamlog_out( DEBUG4 )  <<  timer.toString () << std::endl ;

  _nEvt ++ ;

}




/*************************************************************************************************/
void ClupatraProcessor::check( LCEvent * evt ) { 

  /*************************************************************************************************/
  

}

/*************************************************************************************************/
void ClupatraProcessor::pickUpSiTrackerHits( EVENT::LCCollection* trackCol , LCEvent* evt) {
  
  /*************************************************************************************************/
  
  streamlog_out( DEBUG3  ) << " ************ pickUpSiTrackerHits() called - nTracks : " << trackCol->getNumberOfElements() <<std::endl ;
  
  std::map< int , std::list<TrackerHit*> > hLMap ;
  
  UTIL::BitField64 encoder( ILDCellID0::encoder_string ) ; 
  

  if(  parameterSet( "SITHitCollection" ) ) {
    
    LCIterator<TrackerHit> it( evt, _sitColName ) ;
    while( TrackerHit* hit = it.next()  ){
      
      hLMap[ hit->getCellID0() ].push_back(  hit ) ;
    }    
  }
  if(  parameterSet( "VXDHitCollection" ) ) {
    
    LCIterator<TrackerHit> it( evt, _vxdColName ) ;
    while( TrackerHit* hit = it.next()  ){
      
      hLMap[ hit->getCellID0() ].push_back(  hit ) ;
    }    
  }


  streamlog_out( DEBUG3 ) << "  *****  hitMap size : " <<   hLMap.size() << std::endl ;
  
  for( std::map< int , std::list<TrackerHit*> >::iterator it= hLMap.begin(), End = hLMap.end() ; it != End ; ++it ){
    
    encoder.setValue( it->first ) ;
    
    streamlog_out( DEBUG3 ) << "  *****  sensor: " << encoder.valueString()  << " - nHits: " <<  it->second.size()  << std::endl ;
    
  }
  
  int nSITLayers = 0 ;
  try{
    const gear::ZPlanarParameters& gearSIT = Global::GEAR->getSITParameters() ;
    const gear::ZPlanarLayerLayout& gearSITLayout = gearSIT.getZPlanarLayerLayout() ;
    nSITLayers = gearSITLayout.getNLayers() ;
  }catch(gear::UnknownParameterException& ){}

  const gear::ZPlanarParameters& gearVXD = Global::GEAR->getVXDParameters() ;
  const gear::ZPlanarLayerLayout& gearVXDLayout = gearVXD.getZPlanarLayerLayout() ;
  int nVXDLayers = gearVXDLayout.getNLayers() ;

  int nLayers  = nVXDLayers + nSITLayers  ;


  // ============ sort tracks wrt pt (1./omega) ===============
  LCCollectionVec* tv  = dynamic_cast<LCCollectionVec*>(trackCol) ;

  std::sort( tv->begin() , tv->end() ,  PtSort()  ) ;
  
  for(  LCIterator<TrackImpl> it( trackCol ) ;  TrackImpl* trk = it.next()  ; ) 
    {
	

    //--------------------------------------------
    // create a temporary MarlinTrk
    //--------------------------------------------

// #if 0 //===========================================================================================

//       std::auto_ptr<MarlinTrk::IMarlinTrack> mTrk( _trksystem->createTrack()  ) ;

//       double b = Global::GEAR->getBField().at( gear::Vector3D(0.,0.0,0.) ).z()  ;
      
//       const EVENT::TrackState* ts = trk->getTrackState( lcio::TrackState::AtIP ) ; 
      
//       //    const EVENT::TrackState* ts = trk->getClosestTrackState( 0., 0., 0. ) ;
//       //    const IMPL::TrackStateImpl* ts = dynamic_cast<const IMPL::TrackStateImpl*>( trk->getClosestTrackState( 0., 0., 0. ) ) ;
//       //    const IMPL::TrackStateImpl* ts = dynamic_cast<const IMPL::TrackStateImpl*>( trk->getTrackState( EVENT::TrackState::AtOther ) ) ;
//       //  // FIXME:  what do we need here EVENT::TrackState::AtIP  or AtFirstHit ??
      
//       int nHit = trk->getTrackerHits().size() ;
      
//       if( nHit == 0 || ts ==0 )
// 	continue ;
      
//       streamlog_out( DEBUG3  )  << "  -- extrapolate TrackState : " << lcshort( ts )    << std::endl ;
      
//       //need to add a dummy hit to the track
//       mTrk->addHit(  trk->getTrackerHits()[0] ) ; // fixme: make sure we got the right TPC hit here !??
      
      
//       mTrk->initialise( *ts ,  b ,  MarlinTrk::IMarlinTrack::backward ) ;
    
// #else  //===========================================================================================
      // use the MarlinTrk allready stored with the TPC track


      MarlinTrk::IMarlinTrack* mTrk = trk->ext<MarTrk>()  ;
 
      //      const EVENT::TrackState* ts = trk->getTrackState( lcio::TrackState::AtIP ) ; 

      if( mTrk == 0 ){
	
	streamlog_out( DEBUG3  )  << "  ------ null pointer in ext<MarTrk> ! ?? ..... " << std::endl ;
	continue ;
      }

      //#endif //===========================================================================================

    //--------------------------------------------------
    // get intersection points with SIT and VXD layers 
    //-------------------------------------------------

    //    for( int l=nVXDLayers-1 ; l >=0 ; --l) {

    for( int lx=nLayers-1 ; lx >=0 ; --lx) {

      int detID = (  lx >= nVXDLayers  ?  ILDDetID::SIT   :  ILDDetID::VXD  ) ;
      int layer = (  lx >= nVXDLayers  ?  lx - nVXDLayers  :  lx              ) ;

      encoder.reset() ;
      encoder[ ILDCellID0::subdet ] = detID ;
      encoder[ ILDCellID0::layer  ] = layer ;
      int layerID = encoder.lowWord() ;  
      
      gear::Vector3D point ; 
      
      int sensorID = -1 ;

      int intersects = mTrk->intersectionWithLayer( layerID, point, sensorID, MarlinTrk::IMarlinTrack::modeClosest ) ;
      
      encoder.setValue( sensorID )  ;

      streamlog_out( DEBUG3 ) << " *******  pickUpSiTrackerHits - intersection with SIT/VXD layer " << layer 
			      << " intersects:  " << MarlinTrk::errorCode( intersects ) 
			      << " sensorID: " << encoder.valueString() 
			      << std::endl ;
      
      if( intersects == MarlinTrk::IMarlinTrack::success ){
	
	std::list<TrackerHit*>& hL = hLMap[ sensorID ] ;
	
	streamlog_out( DEBUG3 ) << "    **** found candidate hits : " << hL.size()  
				<< "         for point " << point << std::endl ;
	
	double min = 1.e99 ;
	double maxDist = 1. ; //FIXME: make parameter - what is reasonable here ?
	 
	std::list<TrackerHit*>::iterator bestIt = find_smallest( hL.begin(), hL.end() , Distance3D2( point ) , min ) ;
	
	if( bestIt == hL.end() || min  > maxDist ){

	  streamlog_out( DEBUG3 ) << " ######### no close by hit found !! " << std::endl ;
	  continue ; // FIXME: need to limit the number of layers w/o hits !!!!!!
	}

	double deltaChi ;

	streamlog_out( DEBUG3 ) << " will add best matching hit : " << *bestIt << " with distance : " << min << std::endl ;

	int addHit = mTrk->addAndFit( *bestIt , deltaChi, _dChi2Max ) ;
	    
	streamlog_out( DEBUG3 ) << "    ****  best matching hit : " <<  gear::Vector3D( (*bestIt)->getPosition() )  
				<< "         added : " << MarlinTrk::errorCode( addHit )
				<< "   deltaChi2: " << deltaChi 
				<< std::endl ;

	if( addHit ==  MarlinTrk::IMarlinTrack::success ){


	  trk->addHit( *bestIt ) ;
	  hL.erase( bestIt ) ;

	  IMPL::TrackStateImpl tsi ;
	  double chi2N; int ndfN ;

	  mTrk->getTrackState( tsi , chi2N , ndfN ) ; 

	  streamlog_out( DEBUG3  )  << "  -- extrapolate TrackState : " << lcshort( (TrackState*)&tsi )  << "\n" 
				    << " chi2: " << chi2N
				    << " ndfN: " << ndfN    
				    << std::endl ;

	}

      }
    }

    // -------------------------   update the track state ----------------------
    // FIXME: should this be done in processEvent()  ?
    lcio::TrackStateImpl* tsi =  new lcio::TrackStateImpl ;
    double chi2 ;
    int ndf  ;
    const gear::Vector3D ipv( 0.,0.,0. );
    
    
    // get track state at the IP 
    int ret = mTrk->propagate(   ipv, *tsi, chi2, ndf ) ;
    //    int ret = mTrk->extrapolate( ipv, *tsi, chi2, ndf ) ;
    
    if( ret == MarlinTrk::IMarlinTrack::success ){
      
      tsi->setLocation(  lcio::TrackState::AtIP ) ;
      

      // the track state at the IP needs to be the first one
      //  -> we have to copy the whole vector, and then add 
      //     all track states except the old one at the IP ....
      TrackStateVec tsv  = trk->trackStates() ;
      trk->trackStates().clear() ;
      
      trk->addTrackState( tsi ) ;
      
      for( int i=0, N=tsv.size() ; i<N ; ++i ){

	if( tsv[i]->getLocation() == lcio::TrackState::AtIP ) {

	  delete  tsv[i] ;

	}else{

	  trk->addTrackState( tsv[i] ) ;
	}
      } //-----------------------------------------------------------------



      trk->setChi2( chi2 ) ;
      trk->setNdf( ndf ) ;
    } 
  }
}


//====================================================================================================

void ClupatraProcessor::end(){ 
  
  streamlog_out( MESSAGE )  << "ClupatraProcessor::end()  " << name() 
			    << " processed " << _nEvt << " events in " << _nRun << " runs "
			    << std::endl ;
  
}


//====================================================================================================
