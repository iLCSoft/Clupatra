#include "ClupatraNew.h"

#include <time.h>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

#include <math.h>
#include <cmath>

//---- MarlinUtil 
//#include "ClusterShapes.h"
#include "MarlinCED.h"

//---- LCIO ---
#include "IMPL/LCCollectionVec.h"
#include "IMPL/TrackImpl.h"
#include "IMPL/TrackerHitImpl.h"
#include "EVENT/SimTrackerHit.h"
#include "IMPL/LCFlagImpl.h"
#include "UTIL/Operators.h"
#include "UTIL/LCTOOLS.h"

#include "LCIterator.h"


//---- GEAR ----
#include "marlin/Global.h"
#include "gear/GEAR.h"
#include "gear/TPCParameters.h"
#include "gear/PadRowLayout2D.h"
#include "gear/BField.h"


#include "KalTest.h"

using namespace lcio ;
using namespace marlin ;

#include "clupatra.h"



using namespace clupatra ;

ClupatraNew aClupatraNew ;


ClupatraNew::ClupatraNew() : Processor("ClupatraNew") {
  
  // modify processor description
  _description = "ClupatraNew : simple nearest neighbour clustering" ;
  
  
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


void ClupatraNew::init() { 

  // usually a good idea to
  printParameters() ;


  _kalTest = new KalTest( *marlin::Global::GEAR ) ;

  _kalTest->setOption( KalTest::CFG::ownsHits , false ) ;
  _kalTest->setOption( KalTest::CFG::useQMS   , true ) ;
  _kalTest->setOption( KalTest::CFG::usedEdx  , true ) ;
  
  _kalTest->init() ;

  _nRun = 0 ;
  _nEvt = 0 ;
}

void ClupatraNew::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 

void ClupatraNew::processEvent( LCEvent * evt ) { 

  clock_t start =  clock() ; 

  GHitVec ghits ;
  ghits.setOwner( true ) ;

  std::vector<ClupaHit*> hits ;
  
  //----
  GHitVec h ;
  GHitVec hSmallR ; 
  //----

  GClusterVec cluList ;
  
  RCut<ClupaHit> rCut( _rCut ) ;
  RCutInverse<ClupaHit> rCutInverse( _rCut ) ;
  
  static const int ZBins = 160 ; // with 200 bins we can miss tracks in the very forward region  
  ZIndex<TrackerHit,ZBins > zIndex( -2750. , 2750. ) ; 
  
  HitDistance dist0( _distCut ) ;
  HitDistance dist( 20. ) ;
  //  HitDistance_2 dist_2( 20. ) ;
  

  LCIOTrack<ClupaHit> converter ;
  
  const gear::TPCParameters& gearTPC = Global::GEAR->getTPCParameters() ;
  const gear::PadRowLayout2D& padLayout = gearTPC.getPadLayout() ;
  unsigned nPadRows = padLayout.getNRows() ;


  //------ create clupa hits for the lcio hits with additional parameters
  //note: there are 3 additional hit objects for created for every  lcio::TrackerHit:
  //   - ClupaHit
  //   - GenericHit<ClupaHit> 
  //   - FitHit ( TVTrackhit ) 
  // the ClupaHit holds additional information and is used as  GenericHit<ClupaHit> in the clustering while
  // it also holds a pointer to the hit used in the fit program (and created by KalTest)


  for( StringVec::iterator it = _colNames.begin() ; it !=  _colNames.end() ; it++ ){  

    LCCollectionVec* col = 0 ;

    try{

      col =  dynamic_cast<LCCollectionVec*> (evt->getCollection( *it )  ); 

    }
    catch( lcio::DataNotAvailableException& e) { 
      continue ;
    } 
      
      int nHit = col->getNumberOfElements() ;
      
      hits.reserve(  hits.size() + nHit )   ;
      ghits.reserve(  hits.size() + nHit )   ;
      
      for(int i=0 ; i < nHit ; ++i ) {
	
	TrackerHit* th = (TrackerHit*) col->getElementAt(i) ;
	
	ClupaHit* ch = new ClupaHit ; 
	hits.push_back( ch ) ;
	
	GHit* gh = new GHit( ch ) ;
	ghits.push_back( gh ) ;
	
	th->ext<HitInfo>() = ch ;  // assign the clupa hit to the LCIO hit for memory mgmt
	
	ch->lcioHit = th ; 
	
	ch->pos = Vector3D(  th->getPosition() ) ;
	
	int padIndex = padLayout.getNearestPad( ch->pos.rho() , ch->pos.phi() ) ;
	
	ch->layerID = padLayout.getRowNumber( padIndex ) ;
	
	// create the hit needed by the fitter 
	ch->fitHit  = _kalTest->createHit( th , ch->layerID  , KalTest::DetID::TPC ) ; 
	
	ch->zIndex = zIndex( th ) ;
	
	// ch->phiIndex = ....
	
      } //-------------------- 
      
    }


    //--------------------------

  std::sort( hits.begin(), hits.end() , ZSort() ) ;


  // create two vectors with inner and outer hits:
  addToGenericHitVec( h, hits.begin() , hits.end() , rCut ,  z_index ) ;
  addToGenericHitVec( hSmallR, hits.begin() , hits.end() , rCutInverse ,  z_index ) ;
  
  
  // cluster the sorted hits  ( if |diff(z_index)|>1 the loop is stopped)
  cluster_sorted( h.begin() , h.end() , std::back_inserter( cluList )  , &dist0 , _minCluSize ) ;
  
  streamlog_out( DEBUG ) << "   ***** clusters: " << cluList.size() << std::endl ; 

  //----- save to lcio collection for debugging
  LCCollectionVec* allClu = new LCCollectionVec( LCIO::TRACK ) ;
  std::transform(cluList.begin(), cluList.end(), std::back_inserter( *allClu ) , converter ) ;
  evt->addCollection( allClu , "AllTrackClusters" ) ;


  // find 'odd' clusters that have duplicate hits in pad rows
  GClusterVec ocs ;

  split_list( cluList, std::back_inserter(ocs),  DuplicatePadRows( nPadRows, _duplicatePadRowFraction  ) ) ;


  LCCollectionVec* oddCol = new LCCollectionVec( LCIO::TRACK ) ;
  std::transform( ocs.begin(), ocs.end(), std::back_inserter( *oddCol ) , converter ) ;
  evt->addCollection( oddCol , "OddClu_1" ) ;


  streamlog_out( DEBUG ) << "   ***** clusters: " << cluList.size()  << "   ****** oddClusters " << ocs.size() 
			 << std::endl ; 



  //-------------------- split up cluster with duplicate rows ---------------------------------------------------------


  typedef GClusterVec::iterator GCVI ;



  LCCollectionVec* oddCol2 = new LCCollectionVec( LCIO::TRACK ) ;
  evt->addCollection( oddCol2 , "OddClu_2" ) ;
  LCCollectionVec* oddCol3 = new LCCollectionVec( LCIO::TRACK ) ;
  evt->addCollection( oddCol3 , "OddClu_3" ) ;
  LCCollectionVec* oddCol2_1 = new LCCollectionVec( LCIO::TRACK ) ;
  evt->addCollection( oddCol2_1 , "OddClu_2_1" ) ;
  LCCollectionVec* oddCol3_1 = new LCCollectionVec( LCIO::TRACK ) ;
  evt->addCollection( oddCol3_1 , "OddClu_3_1" ) ;
  LCCollectionVec* oddCol4 = new LCCollectionVec( LCIO::TRACK ) ;
  evt->addCollection( oddCol4 , "OddClu_4" ) ;

  int _nRowForSplitting = 10 ; //FIXME:  make proc param
  
  for( GCVI it = ocs.begin() ; it != ocs.end() ; ++it ){  // loop over all odd clusters ...
    
    GClusterVec sclu ; // new split clusters
    GClusterVec socs ; // split odd clusters
    
    GHitVec oddHits ;
    oddHits.reserve( (*it)->size() ) ;
    
    (*it)->takeHits( std::back_inserter( oddHits )  ) ;
    delete (*it) ;

    unsigned nOddHits = oddHits.size() ;

    
    GHitListVector hitsInLayer( _kalTest->maxLayerIndex() ) ;
    addToHitListVector(  oddHits.begin(), oddHits.end() , hitsInLayer ,  _kalTest->indexOfFirstLayer( KalTest::DetID::TPC)  ) ;
    
    //---------------------------------------------------------------------------------------------------------------------------------
    const bool use_stub_merger = true ;
    if( use_stub_merger ) {
      
      // try clustering in pad row ranges to find clear cluster segments...
      
      int padRowRange = 12 ;  //FIXME: make parameter
      int outerRow = nPadRows + _kalTest->indexOfFirstLayer( KalTest::DetID::TPC)  ;
      
      while( outerRow > padRowRange ) {
	
	oddHits.clear() ;
   	// add all hits in pad row range to oddHits
	for(int iRow = outerRow ; iRow > ( outerRow - padRowRange) ; --iRow ) {
	  
	  std::copy( hitsInLayer[ iRow ].begin() , hitsInLayer[ iRow ].end() , std::back_inserter( oddHits )  ) ;
	}
	
	//----- recluster in given pad row range
	cluster( oddHits.begin(), oddHits.end() , std::back_inserter( sclu ), &dist , _minCluSize ) ;
	
	split_list( sclu, std::back_inserter(socs),  DuplicatePadRows( nPadRows, _duplicatePadRowFraction  ) ) ;
	
	for( GCVI it = socs.begin() ; it != socs.end() ; ++it ){ (*it)->freeHits() ; delete (*it) ; } socs.clear() ;
	
	for( GClusterVec::iterator icv = sclu.begin() ; icv != sclu.end() ; ++ icv ) {
	  (*icv)->ext<ClusterInfo>() = new ClusterInfoStruct ;
	}
	
	std::list< KalTrack* > ktracks ;
	KalTestFitter<KalTest::OrderOutgoing, KalTest::FitBackward > fitter( _kalTest ) ;
	std::transform( sclu.begin(), sclu.end(), std::back_inserter( ktracks ) , fitter ) ;
	
	for( GClusterVec::iterator icv = sclu.begin() ; icv != sclu.end() ; ++ icv ) {
	  addHitsAndFilter( *icv , hitsInLayer , 35. , 100.,  3 ) ; 
	  static const bool backward = true ;
	  addHitsAndFilter( *icv , hitsInLayer , 35. , 100.,  3 , backward ) ; 
	}
	
	std::for_each( ktracks.begin() , ktracks.end() , delete_ptr<KalTrack> ) ;
	
	// merge the good clusters to final list
	cluList.merge( sclu ) ;

	outerRow -= padRowRange ;
	
      } //while outerRow > padRowRange 

    }

    //---------------------------------------------------------------------------------------------------------------------------------
    ///xxX


    //---------------------------------------------------------------------------------------------------------------------------------

    const bool recluster_in_pad_rows = false ;

    if( recluster_in_pad_rows ) {
      
      //========================== first iteration ================================================
      
      // reset the hits index to row ranges for reclustering
       for(unsigned i=0 ; i< nOddHits ; ++i){
	int layer =  oddHits[i]->first->layerID  ;
	oddHits[i]->Index0 =   2 * int( layer / _nRowForSplitting ) ;
      }
      
      //----- recluster in pad row ranges
      cluster( oddHits.begin(), oddHits.end() , std::back_inserter( sclu ), &dist , _minCluSize ) ;
      
      std::transform( sclu.begin(), sclu.end(), std::back_inserter( *oddCol2 ) , converter ) ;
      
      streamlog_out( DEBUG ) << "   ****** oddClusters fixed" << sclu.size()   << std::endl ; 
      
      //--------- remove pad row range clusters where merge occured 
      split_list( sclu, std::back_inserter(socs), DuplicatePadRows( nPadRows, _duplicatePadRowFraction  ) ) ;
      
      
      std::transform( socs.begin(), socs.end(), std::back_inserter( *oddCol3 ) , converter ) ;
      
      
      for( GCVI it = socs.begin() ; it != socs.end() ; ++it ){
	//      (*it)->takeHits( std::back_inserter( oddHits )  ) ;
	(*it)->freeHits() ;
	delete (*it) ;
      }
      socs.clear() ;
      
      
      //========================== second iteration in shifted pad row ranges ================================================
      
      
      oddHits.clear() ;
      for( GCVI it = sclu.begin() ; it != sclu.end() ; ++it ){
	(*it)->takeHits( std::back_inserter( oddHits )  ) ;
	delete (*it) ;
      }
      sclu.clear() ;
      
      // reset the hits index to row ranges for reclustering
      nOddHits = oddHits.size() ;
      
      streamlog_out( DEBUG ) << "   left over odd hits for second iteration of pad row range clustering " << nOddHits << std::endl ;
      
      for(unsigned i=0 ; i< nOddHits ; ++i){
	int layer =  oddHits[i]->first->layerID  ;
	oddHits[i]->Index0 =  2 * int( 0.5 +  ( (float) layer / (float) _nRowForSplitting ) ) ;
      }
      
      //----- recluster in pad row ranges
      cluster( oddHits.begin(), oddHits.end() , std::back_inserter( sclu ), &dist , _minCluSize ) ;
      
      std::transform( sclu.begin(), sclu.end(), std::back_inserter( *oddCol2_1 ) , converter ) ;
      
      streamlog_out( DEBUG ) << "   ****** oddClusters fixed" << sclu.size() 
			     << std::endl ; 
      
      //--------- remove pad row range clusters where merge occured 
      split_list( sclu, std::back_inserter(socs), DuplicatePadRows( nPadRows, _duplicatePadRowFraction  ) ) ;
      
      
      std::transform( socs.begin(), socs.end(), std::back_inserter( *oddCol3_1 ) , converter ) ;
      
      //----------------end  split up cluster with duplicate rows 
      
      for( GCVI it = socs.begin() ; it != socs.end() ; ++it ){
	(*it)->takeHits( std::back_inserter( oddHits )  ) ;
	delete (*it) ;
      }
      socs.clear() ;
      
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      
      // --- recluster the good clusters w/ all pad rows
      
      oddHits.clear() ;
      for( GCVI it = sclu.begin() ; it != sclu.end() ; ++it ){
	(*it)->takeHits( std::back_inserter( oddHits )  ) ;
	delete (*it) ;
      }
      sclu.clear() ;
      
      //   reset the index for 'good' hits coordinate again...
      nOddHits = oddHits.size() ;
      for(unsigned i=0 ; i< nOddHits ; ++i){
	oddHits[i]->Index0 = oddHits[i]->first->zIndex ; //z_index ( oddHits[i]->first ) ;
      }
      
      cluster( oddHits.begin(), oddHits.end() , std::back_inserter( sclu ), &dist , _minCluSize ) ;
      
      std::transform( sclu.begin(), sclu.end(), std::back_inserter( *oddCol4 ) , converter ) ;
      
      // --- end recluster the good clusters w/ all pad rows

      // merge the good clusters to final list
      cluList.merge( sclu ) ;

    } // recluster_in_pad_rows
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    

  }
  ocs.clear() ;


  LCCollectionVec* cluCol = new LCCollectionVec( LCIO::TRACK ) ;
  std::transform( cluList.begin(), cluList.end(), std::back_inserter( *cluCol ) , converter ) ;
  evt->addCollection( cluCol , "CluTrackSegments" ) ;

  // //DEBUG ..... check if there are really no duplicate pad rows ...
  // ocs.clear() ; 
  // split_list( cluList, std::back_inserter(ocs), DuplicatePadRows( nPadRows, _duplicatePadRowFraction  ) ) ;
  // LCCollectionVec* dupCol = new LCCollectionVec( LCIO::TRACK ) ;
  // std::transform( ocs.begin(), ocs.end(), std::back_inserter( *dupCol ) , converter ) ;
  // evt->addCollection( dupCol , "DuplicatePadRowCluster" ) ;

  // streamlog_out( DEBUG ) << "   DuplicatePadRowCluster.size() : " << dupCol->getNumberOfElements() << std::endl ;


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

  //=======================================================================================================  
  static const bool use_best_track = false ;
  static const bool use_best_hit = false ;
  
  if( use_best_track ) {

    streamlog_out( DEBUG ) << "  ------ assign left over hits - best matching track for every hit ..."  << std::endl ;

    Chi2_RPhi_Z ch2rz( 0.1 , 1. ) ; // fixme - need proper errors ....
    Chi2_RPhi_Z_Hit ch2rzh ;


    HitLayerID  tpcLayerID( _kalTest->indexOfFirstLayer( KalTest::DetID::TPC )  ) ;
    
    for( GHVI ih = leftOverHits.begin() ; ih != leftOverHits.end() ; ++ih ){
      
      GHit* hit = *ih ;
      const gear::Vector3D& hPos =  hit->first->pos ;
      
      double ch2Min = 999999999999999. ;
      KalTrack* bestTrk = 0 ;
      
      for( std::list< KalTrack* >::iterator it = ktracks.begin() ; it != ktracks.end() ; ++it ){
	
	const gear::Vector3D* kPos = (*it)->getXingPointForLayer( tpcLayerID( hit ) ) ;
	
	// double rh  =  hPos.rho() ;
	// double rk  =  kPos->rho() ;
	// if( std::abs( rh - rk ) > 0.1 ) {
	// 	streamlog_out( WARNING ) << " --- different radii for hit and crossing point : " <<  tpcLayerID( hit ) << ": " << rh << " - " << rk 
	// 				 <<  *kPos  << std::endl ;
	// } 
	
	if( kPos != 0 ){
	  
	  //	  double ch2 = ch2rz( hPos , *kPos )  ;
	  double ch2 = ch2rzh( hit->first , *kPos )  ;
	  
	  if( ch2 < ch2Min ){
	    
	    ch2Min = ch2 ;
	    bestTrk = *it ;
	  }
	  
	}
	
	// else {
	// 	streamlog_out( MESSAGE ) << " --- no crossing point found for layer : " <<  tpcLayerID( hit ) << ": " << hPos << std::endl ;
	// }
	
      }
      if( bestTrk ) {
	
	const gear::Vector3D* kPos = bestTrk->getXingPointForLayer( tpcLayerID( hit ) ) ;
	
	// double rh  =  hPos.rho() ;
	// double rk  =  kPos->rho() ;
	// if( std::abs( rh - rk ) > 0.1 ) {
	// 	streamlog_out( WARNING ) << "  different radii for hit and crossing point : " << rh << " - " << rk << std::endl ;
	// } 
	
	//      if( std::abs( hPos.rho() - kPos->rho() ) < 0.5 &&   std::abs( hPos.z() - kPos->z() ) < 5. ) {
	
	if(  (  hPos - *kPos ).r()  < 3. ) {   // check for bad outliers... FIXME: need proper criterion here .....
	  
	  
	  GCluster* clu = bestTrk->getCluster< GCluster >() ;
	  
	  streamlog_out( DEBUG ) << " ---- assigning left over hit : " << hPos << " <-> " << *kPos  
				 <<   " dist: " <<  (  hPos - *kPos ).r()  << std::endl ;
	  
	  clu->addHit( hit ) ;
	}	
	else 
	  streamlog_out( DEBUG ) << " ---- NOT assigning left over hit : " << hPos << " <-> " << *kPos << std::endl ;
      }
      else
	streamlog_out( DEBUG ) << " ---- NO best track found ??? ---- " << std::endl ;
      
    }
    
  }
  //        ==========================================================================================
  //} else { // ================== use best matching hit for every track segment =========================
  //        ==========================================================================================
  if( use_best_hit ) {


    streamlog_out( DEBUG1 ) << "  ------ assign left over hits - best matching hit for every track ..."  << std::endl ;
    
    HitLayerID  tpcLayerID( _kalTest->indexOfFirstLayer( KalTest::DetID::TPC )  ) ;
    

    //------------- create vector of left over hits per layer
    typedef std::list<GHit*> HitList ;
    typedef std::vector< HitList > HitListVector ;
    HitListVector hitsInLayer( _kalTest->maxLayerIndex() ) ;
    
    
    for( GHVI ih = leftOverHits.begin() ; ih != leftOverHits.end() ; ++ih ) {
      
      GHit* hit = *ih ;
      //      std::cout << " ++++++  layerId: " << tpcLayerID( hit ) << " max layer index : " <<  _kalTest->maxLayerIndex() << std::endl  ;
      hitsInLayer[ tpcLayerID( hit ) ].push_back( hit )  ;
    }
    //-----------------------------
    
    std::map< GCluster* , KalTrack* > clu2trkMap ;

    const bool use_segment_hits = false ; //true ;
    
    if( use_segment_hits  ){
      
      // store first and last hit of every segment in map with leftover hits in this layer
      
      for( GClusterVec::iterator icv = cluList.begin() ; icv != cluList.end() ; ++ icv ) {
	
	GHit* h0 = (*icv)->front() ;
	GHit* h1 = (*icv)->back() ;
	
	hitsInLayer[ tpcLayerID( h0 ) ].push_back( h0 )  ;
	hitsInLayer[ tpcLayerID( h1 ) ].push_back( h1 )  ;
      }
      
      // sort the tracks wrt. lenghts (#hits)
      ktracks.sort( KalTrackLengthSort() ) ;

      // store assoaciation between cluster and track 
      for( std::list< KalTrack* >::iterator it = ktracks.begin() ; it != ktracks.end() ; ++it ){
	GCluster* c = (*it)->getCluster< GCluster >() ;
	clu2trkMap[ c ] = *it ;
      }	   
    }
    //-------------------------------
    

    //    Chi2_RPhi_Z ch2rz( 0.1 , 1. ) ; // fixme - need proper errors 
    Chi2_RPhi_Z_Hit  ch2rzh ;
    

    for( std::list< KalTrack* >::iterator it = ktracks.begin() ; it != ktracks.end() ; ++it ){
      
      KalTrack* theTrack = *it ;
      if( theTrack == 0 ) 
	continue ;
      
      
      // ----- define chi2 cut    ~15 for 1 GeV pt 
      double chi2Cut = 100000. / ( std::log(1.) - std::log( std::abs(theTrack->getOmega()) ) ) ;


      streamlog_out( DEBUG3 ) << " ------- searching for leftover hits for track : " << theTrack 
			      << "   chi2 cut : " << chi2Cut  << " -  omega : " << theTrack->getOmega() <<  std::endl ;
      
      int xpLayer = 0 ;
      
      // const PointList& xptList = theTrack->getXingPoints() ;
      // for(PointList::const_iterator itXP = xptList.begin() ; itXP != xptList.end() ; ++itXP , xpLayer++ ) {
      // 	const gear::Vector3D* kPos =  *itXP ;
      
      PointList& xpVec = theTrack->getXingPoints() ;
      for( unsigned ixp=0 ; ixp < xpVec.size() ; ++ixp, xpLayer++  ) {
      	const gear::Vector3D* kPos =  xpVec[ixp]  ;
	
	if( kPos == 0 ) {   // we don't have a xing point
	  continue ;
	}
	
       	double ch2Min = 10e99 ;
	GHit* bestHit = 0 ;
	
	HitList& hLL = hitsInLayer.at( xpLayer ) ;
	
	for( HitList::const_iterator ih = hLL.begin() ; ih != hLL.end() ; ++ih ){
	  
	  GHit* hit = *ih ;
	  
	  //VecFromArray hPos(  hit->first->getPosition() ) ;
	  //double ch2 = ch2rz( hPos.v() , *kPos )  ;
	  double ch2 = ch2rzh( hit->first , *kPos )  ;

	  if( ch2 < ch2Min ){
	    
	    ch2Min = ch2 ;
	    bestHit = hit ;
	  }
	}
	
	if( bestHit != 0 ) {
	  
	  const gear::Vector3D& hPos = bestHit->first->pos  ;
	  
	  //	  if( ch2Min  <  6. ) { // Sum( pdf(ch2,ndf==2) )_0^6 ~ 95% )
	  //	  if( ch2Min  <  20. ) { // Sum( pdf(ch2,ndf==2) )_0^20 ~ 99.x% ?? ) // FIXME: need steering parameter and optimize value
	  
	  
	  bestHit->first->chi2Residual = ch2Min ;


	  if( ch2Min  < chi2Cut ) { 
	    
	    streamlog_out( DEBUG1 ) <<   " ---- assigning left over hit : " << hPos << " <-> " << *kPos
				    <<   " dist: " <<  (  hPos - *kPos ).r()
				    <<   " chi2: " <<  ch2Min 
				    <<   "  hit errors :  rphi=" <<  sqrt( bestHit->first->lcioHit->getCovMatrix()[0] 
									   + bestHit->first->lcioHit->getCovMatrix()[2] ) 
				    <<	 "  z= " <<  sqrt( bestHit->first->lcioHit->getCovMatrix()[5] )
				    << std::endl ;
	    
	    
	    if( bestHit->second != 0 ) { //--------------------------------------------------------------------------------
	      
	      // hit is already part of a track segment 
	      
	      
	      GCluster* c = bestHit->second  ;
	      KalTrack* trk = clu2trkMap[ c ] ;
	      

	      if( trk == theTrack ) {
		streamlog_out( ERROR ) << " =======================  found best matching hit from track itself: " 
				       << *bestHit->first->lcioHit
				       <<     std::endl  
				       <<  "      track has  " << trk->getNHits()  << " hits " << std::endl ;

		for( unsigned ii=0 ; ii < xpVec.size() ; ++ii) {
		  if( xpVec[ii] ) 
		    streamlog_out( ERROR ) << "  xing pt : "  << ii << " - " << *xpVec[ii]  ;
		}
		
		
		for( GCluster::iterator its = c->begin(); its != c->end() ; ++its ){
		  GHit* hit = *its ;
		  const gear::Vector3D& hPos = hit->first->pos  ;
		  streamlog_out( ERROR ) << "  hit  : layer: "  <<   tpcLayerID( hit )   << " - " << hPos  ;
		}
		

	      } else {

		
		streamlog_out( DEBUG3 ) << " +++++++++ found best hit already part of a track segment !!!!!! " 
					<< " trk : " << trk  << " #hits: " << trk->getNHits() 
					<< " cluster " << c << c->size() 
					<< std::endl ;   
		
		
		unsigned goodHits(0), allHits(0) ;
		
		double chi2Max = 10. ; // fixme parameter
		
		for( GCluster::iterator its = c->begin(); its != c->end() ; ++its ){
		  
		  ++allHits ;
		  
		  GHit* hit = *its ;
		  const gear::Vector3D& hPos = hit->first->pos ;
		  
		  const gear::Vector3D& kPos = *theTrack->getXingPointForLayer( tpcLayerID( hit ) ) ;
		  
		  if( &kPos != 0 ) {
		    
		    //double chi2 = ch2rz( hPos , *kPos )  ;
		    double chi2 = ch2rzh( hit->first , kPos )  ;

		    streamlog_out( DEBUG3 ) << " +++++++++ chi2 : " << chi2 << ", " << hPos 
					    << " +++++++++        " << kPos 
					    << " +++++++++  hit id " << std::hex << hit->first->lcioHit->id() << std::dec 
					    << std::endl ;
		    
		    if( chi2 < chi2Max ){
		      
		      ++goodHits ;
		    }
		  }
		}
		
		double goodFraction = double( goodHits ) / double(  allHits ) ;
		
		streamlog_out( DEBUG3 ) << " +++++++++ fraction of matching hits : " << goodFraction 
					<< std::endl ;   
		
		
		// ---------------------  merge the track segements -------------------------
		
		if( goodFraction > 0.5  ) { // fixme: what is reasonable here - make parameter ...
		  
		  
		  for( GCluster::iterator its = c->begin(); its != c->end() ; ++its ){

		    delete  xpVec[  tpcLayerID( *its ) ] ; // erase crossing points for these hit
		    xpVec[  tpcLayerID( *its ) ]  = 0 ;   
		  }
		  GCluster* clu = theTrack->getCluster< GCluster >() ;
		  
		  // merge the cluster into the larger one and delete it - remove the hits from the hitsinlayer vector first
		  
		  GHit* h0 = c->front() ;
		  GHit* h1 = c->back() ;
		  
		  hitsInLayer[ tpcLayerID( h0 ) ].remove( h0 )  ;
		  hitsInLayer[ tpcLayerID( h1 ) ].remove( h1 )  ;
		  
		  clu->mergeClusters( c ) ;
		  
		  cluList.remove( c  ) ;
		  
		  
		  streamlog_out( DEBUG3) << " ************ found matching segment, merged all hits: delete cluster : " << c 
					 << " and track : " << trk << std::endl ;
		  
		  delete c ;
		  
		  ktracks.remove( trk ) ;
		  
		} //-------------------------------------------------------------

	      }

		
	    }  else  {  //--------------------------------------------------------------------------------
	      
	      hLL.remove(  bestHit ) ;
	      
	      GCluster* clu = theTrack->getCluster< GCluster >() ;
	      
	      streamlog_out( DEBUG3) << "    ************ found matching hit, add to  cluster : " << clu  << std::endl ;
	      
	      clu->addHit( bestHit ) ;
	    }


	  }
	} 
          // else {
	  //	  streamlog_out( DEBUG1 ) << "????????????????????? no best Hit found xing pnt  : chi2  " << *xpVec[ixp]  << " : " << ch2Min << std::endl ;
	  //	}

      }
    }
    
  }


  //================================================================================================================ 
  // try to create some z stubs in order to get very forward tracks
  
  const bool use_z_stub_merger = true ;
  if( use_z_stub_merger ) {
    
    GHitVec oddHits ;
    GClusterVec sclu ;
    GClusterVec socs ;

    int zIndexRange = 10 ;  //FIXME: make parameter
    int outerZ = ZBins ; 
    
    while( outerZ > 0 ) {
      
      
      streamlog_out( DEBUG3 ) << "====  use_z_stub_merger : z0 " << outerZ << " z1  " << outerZ-zIndexRange << std::endl ;
      
      oddHits.clear() ;



      // add all hits in z range to oddHits
      for(unsigned iRow = 0 ; iRow < nPadRows ; ++iRow ) {
	
	copy_if(  hitsInLayer[ iRow ].begin() , hitsInLayer[ iRow ].end() , std::back_inserter(oddHits), 
		  ZIndexInRange( outerZ-zIndexRange , outerZ ) ) ;
	
	// ZIndexInRange zRange( outerZ-zIndexRange , outerZ )   ;
	
	// if( hitsInLayer[ iRow ].empty() ) 
	//   continue ;
	
	// streamlog_out( DEBUG ) << "====  use_z_stub_merger : hits in layer  " << iRow << "  : " <<   hitsInLayer[ iRow ].size() << std::endl ;
	
	// for( GHitList::iterator iHL =  hitsInLayer[ iRow ].begin() ; iHL !=  hitsInLayer[ iRow ].end() ; ++iHL ){
	
	//   streamlog_out( DEBUG ) << "====  use_z_stub_merger : cecking hit  with z =  " << (*iHL)->first->zIndex   
	// 			 << "  z0 " << outerZ << " z1  " << outerZ-zIndexRange
	// 			 << " -> " << zRange( *iHL )
	// 			 << std::endl ;
	
	//   if( zRange( *iHL ) ){
	//     oddHits.push_back(  *iHL ) ;
	//   }
	// }
      }
      
      streamlog_out( DEBUG3 ) << "====  use_z_stub_merger : copied hits ... # " << oddHits.size() << std::endl ;
      
      
      //----- recluster in given pad row range
      cluster( oddHits.begin(), oddHits.end() , std::back_inserter( sclu ), &dist0 , _minCluSize ) ;
      
      
      streamlog_out( DEBUG3 ) << "====  use_z_stub_merger : clusters found  : " << sclu.size() << std::endl ;
      
      std::transform( sclu.begin(), sclu.end(), std::back_inserter( *oddCol2 ) , converter ) ;
      
      split_list( sclu, std::back_inserter(socs),  DuplicatePadRows( nPadRows, _duplicatePadRowFraction  ) ) ;
      
      for( GCVI it = socs.begin() ; it != socs.end() ; ++it ){ (*it)->freeHits() ; delete (*it) ; } socs.clear() ;
      
      std::transform( sclu.begin(), sclu.end(), std::back_inserter( *oddCol3 ) , converter ) ;

      for( GClusterVec::iterator icv = sclu.begin() ; icv != sclu.end() ; ++ icv ) {
	(*icv)->ext<ClusterInfo>() = new ClusterInfoStruct ;
      }
      
      std::list< KalTrack* > ktracks ;
      KalTestFitter<KalTest::OrderOutgoing, KalTest::FitBackward > fitter( _kalTest ) ;
      std::transform( sclu.begin(), sclu.end(), std::back_inserter( ktracks ) , fitter ) ;
      
      for( GClusterVec::iterator icv = sclu.begin() ; icv != sclu.end() ; ++ icv ) {
	addHitsAndFilter( *icv , hitsInLayer , 35. , 100.,  3 ) ; 
	static const bool backward = true ;
	addHitsAndFilter( *icv , hitsInLayer , 35. , 100.,  3 , backward ) ; 
      }
      
     std::transform( sclu.begin(), sclu.end(), std::back_inserter( *oddCol4 ) , converter ) ;

     std::for_each( ktracks.begin() , ktracks.end() , delete_ptr<KalTrack> ) ;
      
      // merge the good clusters to final list
      cluList.merge( sclu ) ;
      
      outerZ -= zIndexRange ;
      
    } //while outerRow > padRowRange 
    
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

      //      int mult[5] ; for(unisgned i=0;i<5;++i) mult[i] = 0 ;
      std::vector<int> mult(5) ;

      for( unsigned i=0 ; i < hLV.size() ; ++i ){
	unsigned m =  hLV[i].size() ;  

	if( m > 2 )  m = 3 ;

	++mult[ m ] ;
	++mult[ 4] ;
      }

      float total = (  mult[1] + mult[2] + mult[3] ) ; 
      float clumu1 = mult[1] / total ;  
      float clumu2 = mult[2] / total ;

      streamlog_out(  DEBUG3 ) << " leftover cluster multiplicities: (0,1,2,>=3, all) :  [" 
			       <<  mult[0] << " , " <<  mult[1] << " , " <<  mult[2] << " , "    
			       <<  mult[3] << " , " <<  mult[4] << "] "  
			       << " mult_1 = " << clumu1 
			       << " mult_2 = " << clumu2  
			       << std::endl ;    




      // findNearestHits( *clu, nPadRows, 1) ;

      

      if( clumu2 > 0.8 ) {   //FIXME - make parameter
	

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

	
      } else if( clumu1 > 0.9 ) {
	
	// add to final clusters directly
	reclu.push_back( clu ) ;

      }  else {

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
 
  static const bool recluster_left_overs_again = false ;
  if( recluster_left_overs_again ) {
    
    HitDistance distLarge( 3.0 * _distCut ) ;  // FIXME: make parameter 
    
    
    GHitVec oddHits ;
    oddHits.clear() ;

    // add all hits in pad row range to oddHits
    for(int iRow = 0 ; iRow <  _kalTest->maxLayerIndex()  ; ++iRow ) {
      
      std::copy( hitsInLayer[ iRow ].begin() , hitsInLayer[ iRow ].end() , std::back_inserter( oddHits )  ) ;
    }
    
    streamlog_out( DEBUG3 ) <<  " recluster_left_overs_again : oddhits.size() : " << oddHits.size() << std::endl ;


    //----- recluster in given pad row range
    loclu2.clear() ;
    cluster( oddHits.begin(), oddHits.end() , std::back_inserter( loclu2 ), &distLarge , 3 ) ; //_minCluSize ) ;
    
    streamlog_out( DEBUG3 ) <<  " recluster_left_overs_again : loclu2.size() : " << loclu2.size() << std::endl ;

    LCCollectionVec* leftOverCol = new LCCollectionVec( LCIO::TRACK ) ;
    evt->addCollection( leftOverCol , "LeftOverClustersFinal" ) ;
    std::transform( loclu2.begin(), loclu2.end(), std::back_inserter( *leftOverCol ) , converter ) ;
    
  }

  //===============================================================================================
  //  refit all found tracks 
  //==============================


  std::list< KalTrack* > newKTracks ;

  //KalTestFitter<KalTest::OrderIncoming, KalTest::FitForward, KalTest::PropagateToIP > ipFitter( _kalTest ) ;

  KalTestFitter<KalTest::OrderOutgoing, KalTest::FitBackward, KalTest::PropagateToIP > ipFitter( _kalTest ) ;

  //FIXME: DEBUG - non ip fitter
  //  KalTestFitter<KalTest::OrderOutgoing, KalTest::FitBackward > ipFitter( _kalTest ) ;
  
  std::transform( cluList.begin(), cluList.end(), std::back_inserter( newKTracks ) , ipFitter  ) ;


  LCCollectionVec* kaltracks = new LCCollectionVec( LCIO::TRACK ) ;
  LCFlagImpl trkFlag(0) ;
  trkFlag.setBit( LCIO::TRBIT_HITS ) ;
  kaltracks->setFlag( trkFlag.getFlag()  ) ;
  
  std::transform( newKTracks.begin(), newKTracks.end(), std::back_inserter( *kaltracks ) , KalTrack2LCIO() ) ;
  //std::transform( cluList.begin(), cluList.end(), std::back_inserter( *kaltracks ) , converter ) ;

  evt->addCollection( kaltracks , _outColName ) ;
  
 
  //================================================================================================================ 
  //   merge track segments based on track parameters and errors ...
  //
  static const int merge_track_segments = true ;

  if( merge_track_segments ) {

    GenericHitVec<Track> trkVec ;
    GenericClusterVec<Track> trkCluVec ;
    LCCollectionVec* mergedTracks = new LCCollectionVec( LCIO::TRACK ) ;
  
    addToGenericHitVec( trkVec , kaltracks ,  AllwaysTrue()  ) ;

    //    TrackStateDistance trkMerge( 50. ) ;
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
      for( std::list<Track*>::iterator itML = mergedTrk.begin() ; itML != mergedTrk.end() ; ++ itML ){
	
	const TrackerHitVec& hV = (*itML)->getTrackerHits() ;

	for(unsigned i=0 ; i < hV.size() ; ++i){

	  trk->addHit( hV[i] ) ;
	  double chi2ndf = (*itML)->getChi2() / (*itML)->getNdf() ;

	  if( chi2ndf < chi2Min ){
	    bestTrk = (*itML) ;
	    chi2Min = chi2ndf ;
	  }
	}
      }
      if( bestTrk != 0 ){ 

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
      mergedTracks->addElement( trk )  ;

    }

    // add all tracks that have not been merged :
    for( GenericHitVec<Track>::iterator it = trkVec.begin(); it != trkVec.end() ;++it){

      if( (*it)->second == 0 ){

	mergedTracks->addElement(  new TrackImpl( *dynamic_cast<TrackImpl*>( (*it)->first ) ) ) ;
      }
    }


    evt->addCollection( mergedTracks , "MergedKalTracks" ) ;
  }


  //------ register some debugging print funtctions for picking in CED :

  CEDPickingHandler::getInstance().registerFunction( LCIO::TRACKERHIT , &printTrackerHit ) ; 
  CEDPickingHandler::getInstance().registerFunction( LCIO::TRACK , &printTrackShort ) ; 
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
  
  _nEvt ++ ;

  clock_t end = clock () ; 
  
  streamlog_out( DEBUG )  << "---  clustering time: " 
 			  <<  double( end - start ) / double(CLOCKS_PER_SEC) << std::endl  ;
  
}


/*************************************************************************************************/
void ClupatraNew::check( LCEvent * evt ) { 
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


void ClupatraNew::end(){ 
  
  //   std::cout << "ClupatraNew::end()  " << name() 
  // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
  // 	    << std::endl ;
  
  delete _kalTest ;
}


//====================================================================================================
