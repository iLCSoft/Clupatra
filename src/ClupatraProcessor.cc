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



/** helper method cto create a track collections and add it to the event */
inline LCCollectionVec* newTrkCol(const std::string& name, LCEvent * evt ){

  LCCollectionVec* col = new LCCollectionVec( LCIO::TRACK ) ;  

  LCFlagImpl hitFlag(0) ;
  hitFlag.setBit( LCIO::TRBIT_HITS ) ;
  col->setFlag( hitFlag.getFlag()  ) ;

  evt->addCollection( col , name ) ;

  return col ;
}
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
  
  //  clock_t start =  clock() ; 
  Timer timer ;
  unsigned t_init       = timer.registerTimer(" initialization      " ) ;
  unsigned t_seedtracks = timer.registerTimer(" extend seed tracks  " ) ;
  unsigned t_recluster  = timer.registerTimer(" recluster leftovers " ) ;
  unsigned t_split      = timer.registerTimer(" split clusters      " ) ;
  unsigned t_finalfit   = timer.registerTimer(" final refit         " ) ;
  unsigned t_merge      = timer.registerTimer(" merge segments      " ) ;
  
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
  
  //FIXME: parameter !
  static const int ZBins = 160 ; // with 200 bins we can miss tracks in the very forward region  
  ZIndex zIndex( -2750. , 2750. ,ZBins  ) ; 
  
  HitDistance dist( _distCut ) ;
  // HitDistance dist( 20. ) ;
  
  LCIOTrackConverter converter ;
  
  
  const gear::TPCParameters& gearTPC = Global::GEAR->getTPCParameters() ;
  const gear::PadRowLayout2D& padLayout = gearTPC.getPadLayout() ;

  unsigned maxTPCLayers =  padLayout.getNRows() ;
  

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
  //   create some additional (optional) output collections
  //===============================================================================================


  static const bool writeSeedCluster = true ;
  static const bool writeCluTrackSegments = true ;
  static const bool writeLeftoverClusters = true ;
  
  LCCollectionVec* seedCol =  ( writeSeedCluster      ?  newTrkCol( "SeedCluster"      , evt )  :   0   )  ; 
  
  LCCollectionVec* cluCol  =  ( writeCluTrackSegments ?  newTrkCol( "CluTrackSegments" , evt )  :   0   )  ; 
  
  LCCollectionVec* locCol  =  ( writeCluTrackSegments ?  newTrkCol( "LeftoverClusters" , evt )  :   0   )  ; 
  
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

  nnclu::PtrVector<IMarlinTrack> seedTrks ;
  seedTrks.setOwner() ;
  
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
      
      addHitsAndFilter( *icv , hitsInLayer , 35. , 100.,  3 ) ; 
      
      static const bool backward = true ;
      addHitsAndFilter( *icv , hitsInLayer , 35. , 100.,  3 , backward ) ; 
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

  //===============================================================================================
  //  do a global reclustering in all leftover hits
  //===============================================================================================

  Clusterer::cluster_list loclu ; // leftover clusters
  loclu.setOwner() ;
  
  HitVec hits ;

  //  CluTrack debug ;
  
  for(unsigned iRow = 0 ; iRow <  maxTPCLayers ; ++iRow ) {
    
    std::copy( hitsInLayer[ iRow ].begin() , hitsInLayer[ iRow ].end() , std::back_inserter( hits )  ) ;

    //    std::copy( hitsInLayer[ iRow ].begin() , hitsInLayer[ iRow ].end() , std::back_inserter( debug )  ) ;    
  }
  
  nncl.cluster( hits.begin(), hits.end() , std::back_inserter( loclu ),  dist , _minCluSize ) ;
  
  //loclu.push_back( &debug ) ;

  if( writeLeftoverClusters )
    std::transform( loclu.begin(), loclu.end(), std::back_inserter( *locCol ) , converter ) ;
  


  timer.time( t_recluster ) ;

  //===============================================================================================
  //  now we split the clusters based on their hit multiplicities
  //===============================================================================================

  
  
  for( Clusterer::cluster_list::iterator it= loclu.begin(), end= loclu.end() ; it != end ; ++it ){
    
    CluTrack* clu = *it ;

    streamlog_out(  DEBUG2 ) << " **** left over cluster with size : " << clu->size() << std::endl ;
    
    std::vector<int> mult(5) ; 
    // get hit multiplicities up to 3 ( 4 means 4 or higher ) 
    
    getHitMultiplicities( clu , mult ) ;
    
    streamlog_out(  DEBUG2 ) << " **** left over cluster with hit multiplicities: " 
			     << "     m[0] = " <<  mult[0] << "\n" 
			     << "     m[1] = " <<  mult[1] << "\n"
			     << "     m[2] = " <<  mult[2] << "\n"
			     << "     m[3] = " <<  mult[3] << "\n"
			     << "     m[4] = " <<  mult[4] << std::endl ;
    
    
    
    if( float( mult[3]) / mult[0]  >= 0.5 &&  mult[3] >  3 ) {   //FIXME - make parameter
      
      Clusterer::cluster_list reclu ; // reclustered leftover clusters
      reclu.setOwner() ;
      
      create_three_clusters( *clu , reclu ) ;
      
      std::transform( reclu.begin(), reclu.end(), std::back_inserter( seedTrks) , fitter ) ;
      
      for( Clusterer::cluster_list::iterator ir= reclu.begin(), end= reclu.end() ; ir != end ; ++ir ){
	
	streamlog_out( DEBUG2 ) << " extending triplet clustre  of length " << (*ir)->size() << std::endl ;
	
	addHitsAndFilter( *ir , hitsInLayer , 35. , 100.,  3 ) ; 
	static const bool backward = true ;
	addHitsAndFilter( *ir , hitsInLayer , 35. , 100.,  3 , backward ) ; 
      } 
      
      cluList.merge( reclu ) ;
    } 
    
    else if( float( mult[2]) / mult[0]  >= 0.5 &&  mult[2] >  3 ) {   //FIXME - make parameter
      
      Clusterer::cluster_list reclu ; // reclustered leftover clusters
      reclu.setOwner() ;
      
      create_two_clusters( *clu , reclu ) ;
      
      std::transform( reclu.begin(), reclu.end(), std::back_inserter( seedTrks) , fitter ) ;
      
      for( Clusterer::cluster_list::iterator ir= reclu.begin(), end= reclu.end() ; ir != end ; ++ir ){
	
	streamlog_out( DEBUG2 ) << " extending doublet clustre  of length " << (*ir)->size() << std::endl ;
	
	addHitsAndFilter( *ir , hitsInLayer , 35. , 100.,  3 ) ; 
	static const bool backward = true ;
	addHitsAndFilter( *ir , hitsInLayer , 35. , 100.,  3 , backward ) ; 
      } 
      
      cluList.merge( reclu ) ;

    }
    else if( float( mult[1]) / mult[0]  >= 0.5 &&  mult[1] >  3 ) {    


      fitter( *it ) ;

      cluList.push_back( *it ) ;

      it = loclu.erase( it ) ;
      --it ; // erase returns iterator to next element 

    } else {
      
      //  discard cluster and free hits
      clu->freeElements() ; 
    }

  }




  //===============================================================================================
  //  try again to gobble up hits at the ends ....
  //===============================================================================================

  for( Clusterer::cluster_list::iterator icv = cluList.begin() ; icv != cluList.end() ; ++ icv ) {
    addHitsAndFilter( *icv , hitsInLayer , 35. , 100.,  3 ) ; 
    static const bool backward = true ;
    addHitsAndFilter( *icv , hitsInLayer , 35. , 100.,  3 , backward ) ; 
  }
  
  
  timer.time( t_split ) ;

  streamlog_out( MESSAGE ) << " ===========    refitting final " << cluList.size() << " track segments  "   << std::endl ;

  //===============================================================================================
  //  now refit the tracks 
  //===============================================================================================

  nnclu::PtrVector<IMarlinTrack> finalTrks ;
  finalTrks.setOwner() ;
  finalTrks.reserve( cluList.size() ) ;
  
  std::transform( cluList.begin(), cluList.end(), std::back_inserter( finalTrks) , IMarlinTrkFitter(_trksystem)  ) ;
  
  std::transform( cluList.begin(), cluList.end(), std::back_inserter( *tsCol ) , converter ) ;

  
  timer.time( t_finalfit) ;
  
  //===============================================================================================
  //   create collections of used and unused TPC hits 
  //===============================================================================================
  
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
  

  //===============================================================================================
  //  merge track segements 
  //===============================================================================================
  
  
  static const int merge_track_segments = true ;

  if( merge_track_segments ) {

    
    typedef nnclu::NNClusterer<Track> TrackClusterer ;
    TrackClusterer nntrkclu ;

    TrackClusterer::element_vector trkVec ;
    TrackClusterer::cluster_vector trkCluVec ;

    trkVec.reserve( tsCol->size()  ) ;
    trkCluVec.reserve(  tsCol->size() ) ;

    std::for_each( tsCol->begin() , tsCol->end() , ComputeTrackerInfo()  ) ;
    
    std::transform( tsCol->begin() , tsCol->end() , std::back_inserter( trkVec ) , MakeLCIOElement<Track>() ) ; 

    TrackCircleDistance trkMerge( 0.1 ) ; 

    nntrkclu.cluster( trkVec.begin() , trkVec.end() , std::back_inserter( trkCluVec ), trkMerge , 2  ) ;


    streamlog_out( DEBUG3 ) << " ===== merged tracks - # cluster: " << trkCluVec.size()   
			    << " from " << tsCol->size() << " track segments "    << "  ============================== " << std::endl ;
    
    for(  TrackClusterer::cluster_vector::iterator it= trkCluVec.begin() ; it != trkCluVec.end() ; ++it) {
      
      streamlog_out( DEBUG2 ) <<  lcio::header<Track>() << std::endl ;
      
      TrackClusterer::cluster_type*  trkClu = *it ;

      std::list<Track*> mergedTrk ;
      for( TrackClusterer::cluster_type::iterator itC = trkClu->begin() ; itC != trkClu->end() ; ++ itC ){
	
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

      outCol->addElement( trk )  ;
    }

    // add all tracks that have not been merged :
    for( TrackClusterer::element_vector::iterator it = trkVec.begin(); it != trkVec.end() ;++it){

      if( (*it)->second == 0 ){
	
	TrackImpl* t =   new TrackImpl( *dynamic_cast<TrackImpl*>( (*it)->first ) ) ;
	
	t->ext<TrackInfo>() = 0 ; // set extension to 0 to prevent double free ... 
	
	outCol->addElement( t ) ;
      }
    }
    
  }

  //---------------------------------------------------------------------------------------------------------


  timer.time( t_merge ) ;  
  
  streamlog_out( DEBUG4 )  <<  timer.toString () << std::endl ;

  _nEvt ++ ;

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
