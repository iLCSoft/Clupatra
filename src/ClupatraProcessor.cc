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
#include "UTIL/ILDConf.h"

#include "LCIterator.h"


//---- GEAR ----
#include "marlin/Global.h"
#include "gear/GEAR.h"
#include "gear/TPCParameters.h"
#include "gear/PadRowLayout2D.h"
#include "gear/BField.h"


#include "MarlinTrk/Factory.h"
#include "MarlinTrk/IMarlinTrack.h"
#include "MarlinTrk/IMarlinTrkSystem.h"


using namespace lcio ;
using namespace marlin ;
using namespace MarlinTrk ;


#include "clupatra_new.h"
using namespace clupatra_new ;


inline void printSimTrackerHit(const lcio::LCObject* o) ;

inline LCCollectionVec* newTrkCol(const std::string& name, LCEvent * evt ){
  LCCollectionVec* col = new LCCollectionVec( LCIO::TRACK ) ;  
  evt->addCollection( col , name ) ;
  return col ;
}


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
			      (int) 3) ;
  
  
  registerProcessorParameter( "DuplicatePadRowFraction" , 
			      "allowed fraction of hits in same pad row per track"  ,
			      _duplicatePadRowFraction,
			      (float) 0.01 ) ;
  
  
  registerProcessorParameter( "RCut" , 
 			      "Cut for r_min in mm"  ,
 			      _rCut ,
 			      (float) 0.0 ) ;
  

  registerProcessorParameter( "PadRowRange" , 
			      "number of pad rows used in initial seed clustering"  ,
			      _padRowRange ,
			      (int) 12) ;
 
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
  
  
}


void ClupatraProcessor::init() { 

  // usually a good idea to
  printParameters() ;
  
  _trksystem =  MarlinTrk::Factory::createMarlinTrkSystem( "KalTest" , marlin::Global::GEAR , "" ) ;
  
  _trksystem->setOption( IMarlinTrkSystem::CFG::useQMS,        _MSOn ) ;
  _trksystem->setOption( IMarlinTrkSystem::CFG::usedEdx,       _ElossOn) ;
  _trksystem->setOption( IMarlinTrkSystem::CFG::useSmoothing,  _SmoothOn) ;
  _trksystem->init() ;  
  
  _nRun = 0 ;
  _nEvt = 0 ;
  
  // //------ register some debugging print functions for picking in CED :
  // CEDPickingHandler::getInstance().registerFunction( LCIO::TRACKERHIT , &printTrackerHit ) ; 
  // CEDPickingHandler::getInstance().registerFunction( LCIO::TRACK , &printTrackShort ) ; 
  // CEDPickingHandler::getInstance().registerFunction( LCIO::SIMTRACKERHIT , &printSimTrackerHit ) ; 
  
}

void ClupatraProcessor::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 


void ClupatraProcessor::processEvent( LCEvent * evt ) { 
  
  clock_t start =  clock() ; 
  
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
  
  //FIXME: parameter !
  static const int ZBins = 160 ; // with 200 bins we can miss tracks in the very forward region  
  ZIndex zIndex( -2750. , 2750. ,ZBins  ) ; 
  
  HitDistance dist( _distCut ) ;
  // HitDistance dist( 20. ) ;
  
  LCIOTrackConverter converter ;
  
  
  const gear::TPCParameters& gearTPC = Global::GEAR->getTPCParameters() ;
  const gear::PadRowLayout2D& padLayout = gearTPC.getPadLayout() ;

  unsigned maxTPCLayers =  padLayout.getNRows() ;
  

  LCCollectionVec* col = 0 ;

  try{   col =  dynamic_cast<LCCollectionVec*> (evt->getCollection( _colName )  ); 
    
  } catch( lcio::DataNotAvailableException& e) { 
    
    streamlog_out( WARNING ) <<  " input collection not in event : " << _colName << "   - nothing to do  !!! " << std::endl ;  
    
    return ;
  } 
      
  //------ create clupa and clustering hits for every lcio hit -------------------------------------------------

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

  static const bool writeSeedCluster = true ;
  static const bool writeCluTrackSegments = true ;
  
  LCCollectionVec* seedCol =  ( writeSeedCluster      ?  newTrkCol( "SeedCluster"      , evt )  :   0   )  ; 

  LCCollectionVec* cluCol  =  ( writeCluTrackSegments ?  newTrkCol( "CluTrackSegments" , evt )  :   0   )  ; 
  
  
  LCCollectionVec* outCol =  newTrkCol( _outColName  , evt )  ; 
  
  //---------------------------------------------------------------------------------------------------------
  
  // cluster in pad row ranges to find clean cluster segments, statirng from the outside
  // and then moving in to the next more inside pad row range...
  
  Clusterer nncl ;
  
  int outerRow = maxTPCLayers - 1 ;
  
  while( outerRow > _padRowRange ) {
    
    HitVec hits ;
    
    
    // add all hits in pad row range to hits
    for(int iRow = outerRow ; iRow > ( outerRow - _padRowRange) ; --iRow ) {
      
      std::copy( hitsInLayer[ iRow ].begin() , hitsInLayer[ iRow ].end() , std::back_inserter( hits )  ) ;
    }
    
    //-----  cluster in given pad row range  -----------------------------
    Clusterer::cluster_list sclu ;    
    sclu.setOwner() ;  

    nncl.cluster( hits.begin(), hits.end() , std::back_inserter( sclu ), dist , _minCluSize ) ;
    
    if( writeSeedCluster ) {
      std::transform( sclu.begin(), sclu.end(), std::back_inserter( *seedCol ) , converter ) ;
    }

    // remove clusters whith too many duplicate hits per pad row
    Clusterer::cluster_list bclu ;    // bad clusters  
    bclu.setOwner() ;      

    split_list( sclu, std::back_inserter(bclu),  DuplicatePadRows( maxTPCLayers, _duplicatePadRowFraction  ) ) ;

    // free hits from these clusters 
    std::for_each( bclu.begin(), bclu.end(), std::mem_fun( &CluTrack::freeElements ) ) ;
   

    // now fit seed cluster tracks
    nnclu::PtrVector<IMarlinTrack> seedTrks ;
    seedTrks.setOwner() ;

    IMarlinTrkFitter fitter( _trksystem ) ;
    std::transform( sclu.begin(), sclu.end(), std::back_inserter( seedTrks) , fitter ) ;
	

    for( Clusterer::cluster_list::iterator icv = sclu.begin(), end =sclu.end()  ; icv != end ; ++ icv ) {
      
      addHitsAndFilter( *icv , hitsInLayer , 35. , 100.,  3 ) ; 
      
      static const bool backward = true ;
      addHitsAndFilter( *icv , hitsInLayer , 60. , 300.,  3 , backward ) ; 
    } 

    // merge the good clusters to final list
    cluList.merge( sclu ) ;

    outerRow -= _padRowRange ;
    
  } //while outerRow > padRowRange 
  
 
  //---------------------------------------------------------------------------------------------------------

  //  ---- here we have almost complete track segments from the main clupatra algorithm -----  
  if( writeCluTrackSegments )
    std::transform( cluList.begin(), cluList.end(), std::back_inserter( *cluCol ) , converter ) ;

  //---------------------------------------------------------------------------------------------------------
  


  //========  create collections of used and unused TPC hits ===========================================
  
  LCCollectionVec* usedHits   = new LCCollectionVec( LCIO::TRACKERHIT ) ;   ;
  LCCollectionVec* unUsedHits = new LCCollectionVec( LCIO::TRACKERHIT ) ;   ;
  evt->addCollection( usedHits ,   "UsedTPCCluTrackerHits"   ) ;
  evt->addCollection( unUsedHits , "UnUsedTPCCluTrackerHits" ) ;
  usedHits->setSubset() ;
  unUsedHits->setSubset() ;
  usedHits->reserve(   nncluHits.size() ) ;
  unUsedHits->reserve( nncluHits.size() ) ;
  
  for( HitVec::iterator it = nncluHits.begin(), end = nncluHits.end(); it!=end;++it ){
    
    if( (*it)->second != 0 ){   usedHits->push_back( (*it)->first->lcioHit ) ;
    } else {                  unUsedHits->push_back( (*it)->first->lcioHit ) ;          
    }
  }
  

  //---------------------------------------------------------------------------------------------------------

  
  
  //---------------------------------------------------------------------------------------------------------
  // the final track collection:
 
  // now fit the tracks again
  nnclu::PtrVector<IMarlinTrack> finalTrks ;
  finalTrks.setOwner() ;
  finalTrks.reserve( 256 ) ;

  std::transform( cluList.begin(), cluList.end(), std::back_inserter( finalTrks) , IMarlinTrkFitter(_trksystem)  ) ;
  
  //FIXME:   need to merge segments ....

  std::transform( cluList.begin(), cluList.end(), std::back_inserter( *outCol ) , converter ) ;
  
  //---------------------------------------------------------------------------------------------------------









#ifdef IGNORE_THIS_CODE  //===========================================================================================================


  //================================================================================
  
  // create vector with left over hits
  GHitVec leftOverHits ;
  leftOverHits.reserve(  h.size() ) ;
  
  typedef GHitVec::const_iterator GHVI ;
  
  for( GHVI it = h.begin(); it != h.end() ; ++it ){
    
    if ( (*it)->second == 0 ) leftOverHits.push_back( *it ) ;
  }
  
  // add all hits that failed the rcut 
  std::copy( hSmallR.begin() , hSmallR.end() , std::back_inserter( leftOverHits )  ) ;
  
  
  //*********************************************************
  //   run KalTest on track segments (clusters)
  //*********************************************************

  streamlog_out( DEBUG ) <<  "************* fitted segments and KalTest tracks : **********************************" 
			 << std::endl ;


  for( GClusterVec::iterator icv = cluList.begin() ; icv != cluList.end() ; ++ icv ) {
    (*icv)->ext<ClusterInfo>() = new ClusterInfoStruct ;
  }

  std::list< KalTrack* > ktracks ;
  
  //  KalTestFitter<KalTest::OrderIncoming, KalTest::FitForward > fitter( _kalTest ) ;
  KalTestFitter<KalTest::OrderOutgoing, KalTest::FitBackward > fitter( _kalTest ) ;
    
  std::transform( cluList.begin(), cluList.end(), std::back_inserter( ktracks ) , fitter ) ;
  
  // std::for_each( ktracks.begin(), ktracks.end(), std::mem_fun( &KalTrack::findXingPoints ) ) ;
  
  
  if( streamlog_level( DEBUG4 ) ) {
    for( GClusterVec::iterator icv = cluList.begin() ; icv != cluList.end() ; ++ icv ) {
      GCluster* clu  = *icv ;
      KalTrack* trk =  clu->ext<ClusterInfo>()->track ;
      gear::Vector3D xv ;
      int  layer ;
      trk->findNextXingPoint(  xv , layer , 1 ) ;
      
      clu->ext<ClusterInfo>()->nextXPoint = xv ;
      clu->ext<ClusterInfo>()->nextLayer = layer ;
      
      streamlog_out( DEBUG ) <<  "   ----  FINDNEXTXINGPOINT: "  <<  clu
			     <<  " next xing point at layer: "   <<  clu->ext<ClusterInfo>()->nextLayer
			     << " : " <<  clu->ext<ClusterInfo>()->nextXPoint ;
    }
  }

  LCCollectionVec* trksegs = new LCCollectionVec( LCIO::TRACK ) ;
  std::transform( ktracks.begin(), ktracks.end(), std::back_inserter( *trksegs ) , KalTrack2LCIO() ) ;
  evt->addCollection( trksegs , "KalTrackSegments" ) ;

  
  
  //===========  merge track segments based on xing points ==================================================
  
  
  
  //=========== assign left over hits ... ==================================================================
  

  //------------- create vector of left over hits per layer
  GHitListVector hitsInLayer( _kalTest->maxLayerIndex() ) ;

  addToHitListVector(  leftOverHits.begin(), leftOverHits.end() , hitsInLayer ,  _kalTest->indexOfFirstLayer( KalTest::DetID::TPC)  ) ;


  for( GClusterVec::iterator icv = cluList.begin() ; icv != cluList.end() ; ++ icv ) {

    addHitsAndFilter( *icv , hitsInLayer , 35. , 100.,  3 ) ; 

    static const bool backward = true ;
    addHitsAndFilter( *icv , hitsInLayer , 35. , 100.,  3 , backward ) ; 
  }

  //===============================================================================================
  // recluster in the leftover hits
  GClusterVec loclu ; // leftover clusters
 
  static const bool recluster_left_overs = true ;
  if( recluster_left_overs ) {

    HitDistance dist(  _distCut ) ;  // FIXME: make parameter 


    GHitVec oddHits ;
    
    oddHits.clear() ;
    // add all hits in pad row range to oddHits
    for(unsigned iRow = 0 ; iRow <  _kalTest->maxLayerIndex()  ; ++iRow ) {
      
      std::copy( hitsInLayer[ iRow ].begin() , hitsInLayer[ iRow ].end() , std::back_inserter( oddHits )  ) ;
    }
    
    //----- recluster in given pad row range
    loclu.clear() ;
    cluster( oddHits.begin(), oddHits.end() , std::back_inserter( loclu ), &dist , _minCluSize ) ;
    
    LCCollectionVec* leftOverCol = new LCCollectionVec( LCIO::TRACK ) ;
    evt->addCollection( leftOverCol , "LeftOverClusters" ) ;
    std::transform( loclu.begin(), loclu.end(), std::back_inserter( *leftOverCol ) , converter ) ;
    
  }
  //===============================================================================================
  // 
  GClusterVec reclu ; // leftovers reclustered

  static const bool refit_leftover_hits = true ;
  if( refit_leftover_hits ) {

    //    NearestHitDistance nnDist(0.) ;
 
    // compute the hit multiplicities in pad rows

    // use a while loop to remove all clusters at the end (erase the pointer) 
    GClusterVec::iterator icv = loclu.begin() ; 
    while( icv != loclu.end() ) {

      GCluster* clu = *icv ;

      
      GHitListVector hLV( nPadRows )  ;
      addToHitListVector( clu->begin(), clu->end(),  hLV ) ;

      std::vector<int> mult(6) ;

      for( unsigned i=0 ; i < hLV.size() ; ++i ){
	unsigned m =  hLV[i].size() ;  

	if( m >= 4 )  m = 4 ;

	++mult[ m ] ;
	++mult[ 5 ] ;
      }

      float total = (  mult[1] + mult[2] + mult[3] + mult[4] ) ; 
      float clumu1 = mult[1] / total ;  
      float clumu2 = mult[2] / total ;
      float clumu3 = mult[3] / total ;

      streamlog_out(  DEBUG3 ) << " leftover cluster multiplicities: (0,1,2,3,>=4, all) :  [" 
			       <<  mult[0] << " , " <<  mult[1] << " , " <<  mult[2] << " , "    
			       <<  mult[3] << " , " <<  mult[4] << mult[5] << "] "  
			       << " mult_1 = " << clumu1 
			       << " mult_2 = " << clumu2  
			       << std::endl ;    




      // findNearestHits( *clu, nPadRows, 1) ;

      

      if( clumu3 >= 0.5 ) {   //FIXME - make parameter
	

	//---- get hits from cluster into a vector
	GHitVec v ;
	v.reserve( clu->size() ) ;
	clu->takeHits( std::back_inserter( v )  ) ;
	
	
	create_three_clusters( v , reclu, nPadRows ) ;

	GClusterVec::reverse_iterator iC = reclu.rbegin()  ;
	GCluster* clu0 = *iC++ ;
	GCluster* clu1 = *iC++ ;
	GCluster* clu2 = *iC   ;

	clu0->ext<ClusterInfo>() = new ClusterInfoStruct ;
	clu1->ext<ClusterInfo>() = new ClusterInfoStruct ;
	clu2->ext<ClusterInfo>() = new ClusterInfoStruct ;

	KalTestFitter<KalTest::OrderOutgoing, KalTest::FitBackward > fitter( _kalTest ) ;
	
	KalTrack* trk0 = fitter( clu0 ) ;
	KalTrack* trk1 = fitter( clu1 ) ; 
	KalTrack* trk2 = fitter( clu2 ) ; 
	
	// try to extend the clusters with leftover hits (from layers that do not have two hits)
	static const bool backward = true ;
	
	if( trk0->getNHits() > 3 ) {
	  addHitsAndFilter( clu0 , hitsInLayer , 35. , 100.,  3 ) ; 
	  addHitsAndFilter( clu0 , hitsInLayer , 35. , 100.,  3 , backward ) ; 
	}

	if( trk1->getNHits() > 3 ) {
	  addHitsAndFilter( clu1 , hitsInLayer , 35. , 100.,  3 ) ; 
	  addHitsAndFilter( clu1 , hitsInLayer , 35. , 100.,  3 , backward ) ; 
	}

	if( trk2->getNHits() > 3 ) {
	  addHitsAndFilter( clu2 , hitsInLayer , 35. , 100.,  3 ) ; 
	  addHitsAndFilter( clu2 , hitsInLayer , 35. , 100.,  3 , backward ) ; 
	}

	delete trk0 ;
	delete trk1 ;
	delete trk2 ;

	
      } if( clumu2 > 0.5 ) {   //FIXME - make parameter
	

	//---- get hits from cluster into a vector
	GHitVec v ;
	v.reserve( clu->size() ) ;
	clu->takeHits( std::back_inserter( v )  ) ;
	
	
	create_two_clusters( v , reclu, nPadRows ) ;


	GClusterVec::reverse_iterator iC = reclu.rbegin()  ;
	GCluster* clu0 = *iC++ ;
	GCluster* clu1 = *iC   ;

	clu0->ext<ClusterInfo>() = new ClusterInfoStruct ;
	clu1->ext<ClusterInfo>() = new ClusterInfoStruct ;

	KalTestFitter<KalTest::OrderOutgoing, KalTest::FitBackward > fitter( _kalTest ) ;
	
	KalTrack* trk0 = fitter( clu0 ) ;
	KalTrack* trk1 = fitter( clu1 ) ; 
	
	// try to extend the clusters with leftover hits (from layers that do not have two hits)
	static const bool backward = true ;
	
	addHitsAndFilter( clu0 , hitsInLayer , 35. , 100.,  3 ) ; 
	addHitsAndFilter( clu0 , hitsInLayer , 35. , 100.,  3 , backward ) ; 
	
	addHitsAndFilter( clu1 , hitsInLayer , 35. , 100.,  3 ) ; 
	addHitsAndFilter( clu1 , hitsInLayer , 35. , 100.,  3 , backward ) ; 
	
	delete trk0 ;
	delete trk1 ;

	
      } else if( clumu1 > 0.5 ) {
	
	// add to final clusters directly
	reclu.push_back( clu ) ;

      }  else {

	//*********************************************************
	//FIXME: need treatment for leftover clusters with 'undefined' multicplicity .....
	//*********************************************************


	delete clu ;
      }

      loclu.erase( icv++ ) ;  // erase all clusters from the list (GClusterVec)
    }

    // for( GClusterVec::iterator icv = reclu.begin() ; icv != reclu.end() ; ++ icv ) {
    //   (*icv)->ext<ClusterInfo>() = new ClusterInfoStruct ;
    // }
    
    LCCollectionVec* ccol = new LCCollectionVec( LCIO::TRACK ) ;
    std::transform( reclu.begin(), reclu.end(), std::back_inserter( *ccol ) , converter ) ;
    evt->addCollection( ccol , "ReclusteredLeftOvers" ) ;
    
    
    // for now just fit the clusters ...
    cluList.merge( reclu ) ;

  }
  //================================================================================================================ 
  // recluster in the leftover hits
  GClusterVec loclu2 ; // leftover clusters
 

  //===============================================================================================
  //  refit all found tracks 
  //==============================


  std::list< KalTrack* > newKTracks ;

  //KalTestFitter<KalTest::OrderIncoming, KalTest::FitForward, KalTest::PropagateToIP > ipFitter( _kalTest ) ;


  // adding an IP hit does not work for curler segments - as we fit everything leave out the IP hit
  //KalTestFitter<KalTest::OrderOutgoing, KalTest::FitBackward, KalTest::PropagateToIP > ipFitter( _kalTest ) ;
  //FIXME: DEBUG - non ip fitter
  KalTestFitter<KalTest::OrderOutgoing, KalTest::FitBackward > ipFitter( _kalTest ) ;
  
  std::transform( cluList.begin(), cluList.end(), std::back_inserter( newKTracks ) , ipFitter  ) ;


  LCCollectionVec* kaltracks = new LCCollectionVec( LCIO::TRACK ) ;
  LCFlagImpl trkFlag(0) ;
  trkFlag.setBit( LCIO::TRBIT_HITS ) ;
  kaltracks->setFlag( trkFlag.getFlag()  ) ;
  
  std::transform( newKTracks.begin(), newKTracks.end(), std::back_inserter( *kaltracks ) , KalTrack2LCIO() ) ;
  //std::transform( cluList.begin(), cluList.end(), std::back_inserter( *kaltracks ) , converter ) ;

  evt->addCollection( kaltracks , "ClupatraTrackSegments" ) ; 
  
 
  //================================================================================================================ 
  //   merge track segments based on track parameters and errors ...
  //
  static const int merge_track_segments = true ;

  if( merge_track_segments ) {

    GenericHitVec<Track> trkVec ;
    GenericClusterVec<Track> trkCluVec ;
    LCCollectionVec* mergedTracks = new LCCollectionVec( LCIO::TRACK ) ;
    LCFlagImpl trkFlag(0) ;
    trkFlag.setBit( LCIO::TRBIT_HITS ) ;
    mergedTracks->setFlag( trkFlag.getFlag()  ) ;
  
    addToGenericHitVec( trkVec , kaltracks ,  AllwaysTrue()  ) ;

    //    TrackStateDistance trkMerge( 1000. ) ;
    TrackCircleDistance trkMerge( 0.1 ) ; 

    cluster( trkVec.begin() , trkVec.end() , std::back_inserter( trkCluVec ), &trkMerge  , 2 ) ;


    streamlog_out( DEBUG4 ) << " ===== merged tracks - # cluster: " << trkCluVec.size()   << "  ============================== " << std::endl ;
    
    for( GenericClusterVec<Track>::iterator it= trkCluVec.begin() ; it != trkCluVec.end() ; ++it) {
      
      streamlog_out( DEBUG2 ) <<  myheader<Track>() << std::endl ;
      
      GenericCluster<Track>* trkClu = *it ;
      
      std::list<Track*> mergedTrk ;
      for( GenericCluster<Track>::iterator itC = trkClu->begin() ; itC != trkClu->end() ; ++ itC ){
	
	streamlog_out( DEBUG2 ) << lcshort(  (*itC)->first ) << std::endl ; 
	
	mergedTrk.push_back( (*itC)->first ) ; 
      }

      TrackImpl* trk = new TrackImpl ;
      Track* bestTrk = 0 ;
      double chi2Min = 99999999999999999. ;
      int hitCount = 0 ;
      int hitsInFit = 0 ;
      for( std::list<Track*>::iterator itML = mergedTrk.begin() ; itML != mergedTrk.end() ; ++ itML ){
	
	const TrackerHitVec& hV = (*itML)->getTrackerHits() ;
	for(unsigned i=0 ; i < hV.size() ; ++i){
	  trk->addHit( hV[i] ) ;
	  ++hitCount ;
	}

	double chi2ndf = (*itML)->getChi2() / (*itML)->getNdf() ;
	
	if( chi2ndf < chi2Min ){
	  bestTrk = (*itML) ;
	  chi2Min = chi2ndf ;
	}

      }
      if( bestTrk != 0 ){ 

	hitsInFit = trk->getTrackerHits().size() ;

	trk->setD0( bestTrk->getD0() ) ;
	trk->setOmega( bestTrk->getOmega() ) ;
	trk->setPhi( bestTrk->getPhi() ) ;
	trk->setZ0( bestTrk->getZ0() ) ;
	trk->setTanLambda( bestTrk->getTanLambda() ) ;
	trk->setCovMatrix( bestTrk->getCovMatrix()  ) ;
	// ...
	
      }
      else{
	streamlog_out( ERROR ) << "   no best track found for merged tracks ... !? " << std::endl ; 
      }

      trk->subdetectorHitNumbers().push_back( hitCount ) ;  
      trk->subdetectorHitNumbers().push_back( hitsInFit ) ;  

      mergedTracks->addElement( trk )  ;
    }

    // add all tracks that have not been merged :
    for( GenericHitVec<Track>::iterator it = trkVec.begin(); it != trkVec.end() ;++it){

      if( (*it)->second == 0 ){
	
	TrackImpl* t =   new TrackImpl( *dynamic_cast<TrackImpl*>( (*it)->first ) ) ;
	
	t->ext<TrackInfo>() = 0 ; // set extension to 0 to prevent double free ... 
	
	mergedTracks->addElement( t ) ;
      }
    }
    
    evt->addCollection( mergedTracks , _outColName ) ;
  }


  //======================================================================================================

 
  //========== cleanup KalTracks ========
  std::for_each( ktracks.begin() , ktracks.end() , delete_ptr<KalTrack> ) ;

  //FIXME: memory leak - need for debugging....
  std::for_each( newKTracks.begin() , newKTracks.end() , delete_ptr<KalTrack> ) ;
  //=====================================



  //========  create collections of used and unused TPC hits ===========================================

  LCCollectionVec* usedHits = new LCCollectionVec( LCIO::TRACKERHIT ) ;   ;
  LCCollectionVec* unUsedHits = new LCCollectionVec( LCIO::TRACKERHIT ) ;   ;
  usedHits->setSubset() ;
  unUsedHits->setSubset() ;
  usedHits->reserve( h.size() ) ;
  unUsedHits->reserve( h.size() ) ;

  for( GHVI it = h.begin(); it != h.end() ;++it){
    if( (*it)->second != 0 ){
      usedHits->push_back( (*it)->first->lcioHit ) ;
    } else {
      unUsedHits->push_back( (*it)->first->lcioHit ) ;          
    }
  }
  for( GHVI it = hSmallR.begin(); it != hSmallR.end() ;++it){
    if( (*it)->second != 0 ){
      usedHits->push_back( (*it)->first->lcioHit ) ;
    } else {
      unUsedHits->push_back( (*it)->first->lcioHit ) ;          
    }
  }
  evt->addCollection( usedHits ,   "UsedTPCCluTrackerHits" ) ;
  evt->addCollection( unUsedHits , "UnUsedTPCCluTrackerHits" ) ;
  
  //========================================================================================================
#endif // IGNORE_THE_CODE  //===========================================================================================================
  
  _nEvt ++ ;

  clock_t end = clock () ; 
  
  streamlog_out( DEBUG )  << "---  clustering time: " 
 			  <<  double( end - start ) / double(CLOCKS_PER_SEC) << std::endl  ;
  
}


/*************************************************************************************************/
void ClupatraProcessor::check( LCEvent * evt ) { 
  /*************************************************************************************************/
  
  // //  std::string colName( "MergedKalTracks"  ) ;
  // std::string colName(  _outColName ) ;


  // bool checkForDuplicatePadRows =  false ;
  // bool checkForMCTruth          =  true  ;
  // bool checkForSplitTracks      =  true  ; 

  // streamlog_out( MESSAGE ) <<  " check called.... " << std::endl ;

  // const gear::TPCParameters& gearTPC = Global::GEAR->getTPCParameters() ;
  // const gear::PadRowLayout2D& pL = gearTPC.getPadLayout() ;


  // //====================================================================================
  // // check for duplicate padRows 
  // //====================================================================================

  // if( checkForDuplicatePadRows ) {

  //   LCCollectionVec* oddCol = new LCCollectionVec( LCIO::TRACK ) ;
  //   oddCol->setSubset( true ) ;
  //   // try iterator class ...

  //   LCIterator<Track> trIt( evt, colName ) ;
  //   while( Track* tr = trIt.next()  ){

      
  //     // check for duplicate layer numbers
  //     std::vector<int> hitsInLayer( pL.getNRows() ) ; 
  //     const TrackerHitVec& thv = tr->getTrackerHits() ;
  //     typedef TrackerHitVec::const_iterator THI ;
  //     for(THI it = thv.begin() ; it  != thv.end() ; ++it ) {
  // 	TrackerHit* th = *it ;
  // 	++ hitsInLayer.at( th->ext<HitInfo>()->layerID )   ;
  //     } 
  //     unsigned nHit = thv.size() ;
  //     unsigned nDouble = 0 ;
  //     for(unsigned i=0 ; i < hitsInLayer.size() ; ++i ) {
  // 	if( hitsInLayer[i] > 1 ){
  // 	  ++nDouble ;
  // 	  streamlog_out( DEBUG4 ) << " &&&&&&&&&&&&&&&&&&&&&&&&&& duplicate hit in layer : " << i << std::endl ;
  // 	}
  //     }
  //     if( double(nDouble) / nHit > _duplicatePadRowFraction ){
  // 	//if( nDouble  > 0){
  // 	streamlog_out( DEBUG4 ) << " oddTrackCluster found with "<< 100. * double(nDouble) / nHit 
  // 				<< "% of double hits " << std::endl ;
  // 	oddCol->addElement( tr ) ;
  //     }
  //   }
  //   evt->addCollection( oddCol , "OddCluTracks" ) ;
  // }
  // //====================================================================================
  // // check Monte Carlo Truth via SimTrackerHits 
  // //====================================================================================

  // if( checkForMCTruth ) {
 

  //   LCCollectionVec* oddCol = new LCCollectionVec( LCIO::TRACK ) ;
  //   oddCol->setSubset( true ) ;

  //   LCCollectionVec* splitCol = new LCCollectionVec( LCIO::TRACK ) ;
  //   splitCol->setSubset( true ) ;
    
  //   typedef std::map<Track* , unsigned > TRKMAP ; 
    
  //   typedef std::map< MCParticle* , TRKMAP > MCPTRKMAP ;
  //   MCPTRKMAP mcpTrkMap ;
    
  //   typedef std::map< MCParticle* , unsigned > MCPMAP ;
  //   MCPMAP hitMap ;
    
    
  //   // if( streamlog_level( DEBUG4) )
  //   //   LCTOOLS::printTracks( evt->getCollection("KalTestTracks") ) ;


  //   LCIterator<Track> trIt( evt, colName  ) ;  
  //   //    "KalTestTracks" ) ;
  //   //    LCIterator<Track> trIt( evt, _outColName ) ;
  //   //    LCIterator<Track> trIt( evt, "TPCTracks" ) ;

  //   while( Track* tr = trIt.next()  ){
      
  //     MCPMAP mcpMap ;

  //     const TrackerHitVec& thv = tr->getTrackerHits() ;
  //     typedef TrackerHitVec::const_iterator THI ;

  //     // get relation between mcparticles and tracks
  //     for(THI it = thv.begin() ; it  != thv.end() ; ++it ) {

  // 	TrackerHit* th = *it ;
  // 	// FIXME:
  // 	// we know that the digitizer puts the sim hit into the raw hit pointer
  // 	// but of course the proper way is to go through the LCRelation ...
  // 	SimTrackerHit* sh = (SimTrackerHit*) th->getRawHits()[0] ;
  // 	MCParticle* mcp = sh->getMCParticle() ;

	
  // 	hitMap[ mcp ] ++ ;   // count all hits from this mcp
	
  // 	mcpMap[ mcp ]++ ;    // count hits from this mcp for this track
	
  // 	mcpTrkMap[ mcp ][ tr ]++ ;  // map between mcp, tracks and hits
	
  //     } 

  //     // check for tracks with hits from several mcparticles
  //     unsigned nHit = thv.size() ;
  //     unsigned maxHit = 0 ; 
  //     for( MCPMAP::iterator it= mcpMap.begin() ;
  // 	   it != mcpMap.end() ; ++it ){
  // 	if( it->second  > maxHit ){
  // 	  maxHit = it->second ;
  // 	}
  //     }

  //     if( double(maxHit) / nHit < 0.95 ){ // What is acceptable here ???
  // 	//if( nDouble  > 0){
  // 	streamlog_out( MESSAGE ) << " oddTrackCluster found with only "
  // 				 << 100.*double(maxHit)/nHit 
  // 				 << "% of hits  form one MCParticle " << std::endl ;
  // 	oddCol->addElement( tr ) ;
  //     }
  //   }
  //   evt->addCollection( oddCol , "OddMCPTracks" ) ;
    
    
  //   if( checkForSplitTracks ) {
      
  //     streamlog_out( DEBUG ) << " checking for split tracks - mcptrkmap size : " <<  mcpTrkMap.size() << std::endl ;
      
  //     // check for split tracks 
  //     for( MCPTRKMAP::iterator it0 = mcpTrkMap.begin() ; it0 != mcpTrkMap.end() ; ++it0){
	
  // 	streamlog_out( DEBUG ) << " checking for split tracks - map size : " <<  it0->second.size() << std::endl ;
	
	
  // 	if( it0->second.size() > 1 ) {
	  
	  
  // 	  typedef std::list< EVENT::Track* > TL ;
  // 	  TL trkList ;
	  
  // 	  for( TRKMAP::iterator it1 = it0->second.begin() ; it1 != it0->second.end() ; ++it1){
	    
  // 	    double totalHits = hitMap[ it0->first ]  ; // total hits for this track 
	    
  // 	    double thisMCPHits = it1->second ;     //  hits from this mcp
	    
  // 	    double ratio =  thisMCPHits / totalHits  ;
	    
  // 	    streamlog_out( DEBUG ) << " checking for split tracks - ratio : " 
  // 				   << thisMCPHits << " / " << totalHits << " = " << ratio << std::endl ;
	    
  // 	    if( ratio > 0.05 && ratio < 0.95 ){
  // 	      // split track
	      
  // 	      splitCol->addElement( it1->first ) ; 
	      
  // 	      trkList.push_back( it1->first ) ;
  // 	    } 
  // 	  }
  // 	  // chi2 between split track segments :
  // 	  // for( TRKMAP::iterator ist0 = it0->second.begin() ; ist0 != it0->second.end() ; ++ist0){
	    
  // 	  //   KalTrack* sptrk0 = ist0->first->ext<KalTrackLink>() ; 
	    
  // 	  //   TRKMAP::iterator ist0_pp = ist0 ;
  // 	  //   ++ist0_pp ;

  // 	  //   for( TRKMAP::iterator ist1 = ist0_pp ; ist1 != it0->second.end() ; ++ist1){
	  
  // 	  //     KalTrack* sptrk1 = ist1->first->ext<KalTrackLink>() ; 
	      
  // 	  //     double chi2 =  KalTrack::chi2( *sptrk0 ,  *sptrk1 ) ;
	      
  // 	  //     streamlog_out( DEBUG4 ) << " *********************  chi2 between split tracks : "  << chi2 << std::endl 
  // 	  // 			      << myheader< Track >() << std::endl 
  // 	  // 			      << lcshort( ist0->first )  << std::endl 
  // 	  // 			      << lcshort( ist1->first )	 << std::endl ; 
	      
  // 	  //   }
  // 	  // }



  // 	  streamlog_out( DEBUG2 ) << " ------------------------------------------------------ " << std::endl ;
	  
  // 	  for( TL::iterator it0 = trkList.begin() ; it0 != trkList.end() ; ++it0 ){
	    
	    
  // 	    //	    KalTrack* trk0 = (*it0)->ext<KalTrackLink>() ; 
	    
  // 	    HelixClass hel ;
  // 	    hel.Initialize_Canonical( (*it0)->getPhi(),
  // 				      (*it0)->getD0(),
  // 				      (*it0)->getZ0(),
  // 				      (*it0)->getOmega(),
  // 				      (*it0)->getTanLambda(),
  // 				      3.50 ) ;
	    
  // 	    streamlog_out( DEBUG1 ) << hel.getXC() << "\t"
  // 				    << hel.getYC() << "\t"
  // 				    << hel.getRadius() << "\t" 
  // 				    << hel.getTanLambda() << std::endl ; 
	    
	    
  // 	    // streamlog_out( DEBUG1 ) << (*it0)->getPhi() << "\t"
  // 	    // 			  << (*it0)->getD0()  << "\t"
  // 	    // 			  << (*it0)->getOmega()  << "\t"
  // 	    // 			  << (*it0)->getZ0()  << "\t"
  // 	    // 			  << (*it0)->getTanLambda()  << "\t"
  // 	    // 			  << std::endl ;
	    
  // 	    //	    streamlog_out( DEBUG1 ) << " trk0 : " << *trk0 << std::endl ;
	    
  // 	    // TL::iterator its = it0 ;
  // 	    // ++its ;
	    
  // 	    // for( TL::iterator it1 =  its ; it1 != trkList.end() ; ++it1 ){
	      
  // 	    //   KalTrack* trk1 = (*it1)->ext<KalTrackLink>() ; 
	      
  // 	    //   streamlog_out( DEBUG1 ) << "    - trk0 : " << *trk0 << std::endl ;
  // 	    //   streamlog_out( DEBUG1 ) << "    - trk1 : " << *trk1 << std::endl ;
	      
  // 	    //   double chi2 =  KalTrack::chi2( *trk0 ,  *trk1 ) ;
	      
  // 	    //   streamlog_out( DEBUG1 ) << " +++++++++++++++++  chi2 between split tracks : " 
  // 	    // 			      << trk0 << " - " << trk1 << " : " << chi2 << std::endl ; 
	      
	      
  // 	    // }
  // 	  }
	  
  // 	}
  //     }
  //     evt->addCollection( splitCol , "SplitTracks" ) ;
  //   }

  // }
  //====================================================================================

}


void ClupatraProcessor::end(){ 
  
  streamlog_out( MESSAGE )  << "ClupatraProcessor::end()  " << name() 
			    << " processed " << _nEvt << " events in " << _nRun << " runs "
			    << std::endl ;
  
}

//====================================================================================================



inline void printSimTrackerHit(const lcio::LCObject* o){

  lcio::SimTrackerHit* hit = const_cast<lcio::SimTrackerHit*> ( dynamic_cast<const lcio::SimTrackerHit*> (o) ) ; 
  
  if( hit == 0 ) {
    
    streamlog_out( ERROR ) << "  printSimTrackerHit : dynamic_cast<SimTrackerHit*> failed for LCObject : " << o << std::endl ;
    return  ;
  }
  
  streamlog_out( MESSAGE ) << *hit 
			   << " MCParticle id: " << hit->getMCParticle()->id() 
			   << " MCParticle parent id: " 
			   << ( hit->getMCParticle()->getParent(0) ? hit->getMCParticle()->getParent(0)->id() : -1 ) 
			   << std::endl ;
}
