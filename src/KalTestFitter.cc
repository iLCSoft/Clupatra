#include "KalTestFitter.h"

#include <time.h>
#include <vector>
#include <map>
#include <algorithm>


#include "NNClusters_clupa.h"
// #include "ClusterShapes.h"

//---- LCIO ---
#include "IMPL/LCCollectionVec.h"
#include "IMPL/TrackImpl.h"
#include "IMPL/TrackerHitImpl.h"
//#include "EVENT/SimTrackerHit.h"
#include "IMPL/LCFlagImpl.h"
#include "UTIL/Operators.h"

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


typedef GenericCluster<TrackerHit> HitCluster ;
typedef GenericHit<TrackerHit> Hit ;

namespace KTF{
  // delete helper
  template<class P>  void delete_ptr(P* p) { delete p;}
  //-------------------------------
  
  //-------------------------------
  
  
  // helper class to assign additional parameters to TrackerHits
  struct HitInfoStruct{
    HitInfoStruct() :layerID(-1), usedInTrack(false) {}
    int layerID ;
    bool usedInTrack ;
  } ;
  struct HitInfo : LCOwnedExtension<HitInfo, HitInfoStruct> {} ;
  
  
  //------------------------------------------------------
  // function to extract position for Kaltest (in cm):
  TVector3 hitPosition( Hit* h)  { 
    // return TVector3( h->first->getPosition()[0],   
    // 		   h->first->getPosition()[1],
    // 		   h->first->getPosition()[2]  ) ; 
    return TVector3( h->first->getPosition()[0] *.1 ,   // convert to cm ...
		     h->first->getPosition()[1] *.1 ,
		     h->first->getPosition()[2] *.1 ) ; 
  }   
  
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
  
  
  //-------------------------------------------------------------------------
  
  struct LCIOTrackerHit{ EVENT::TrackerHit* operator()( Hit* h) { return h->first ; }   } ;
  
  // Helper struct that calls the KalTest fit on the hits in the cluster
  template <bool HitOrder, bool FitOrder, bool PropagateIP=false>
  struct KTFitter{
    
    KalTest* _kt ; 
    
    KTFitter(KalTest* k) : _kt( k ) {}
    
    KalTrack* operator() (HitCluster* clu) {  
      
      static HitLayerID tpcLayerID( _kt->indexOfFirstLayer( KalTest::DetID::TPC )  )  ;
      
      clu->sort( LayerSort<HitOrder>() ) ;
      
      KalTrack* trk = _kt->createKalTrack() ;
      
      trk->setCluster<HitCluster>( clu ) ;
      
      
      if( PropagateIP  && HitOrder == KalTest::OrderOutgoing ) {
	
	trk->addIPHit() ;
      }  
      
      trk->addHits( clu->begin() , clu->end() , hitPosition, tpcLayerID , LCIOTrackerHit() ) ; 
      
      if( PropagateIP  && HitOrder == KalTest::OrderIncoming ) {
	
	trk->addIPHit() ;
      }  
      
      trk->fitTrack( FitOrder  ) ;
      
      return trk;
    }

  };
  
  struct KalTrack2LCIO{
    TrackImpl* operator() (KalTrack* trk) {  
      TrackImpl* lTrk = new TrackImpl ;
      trk->toLCIOTrack( lTrk  ) ;
      return lTrk ;
    }
  };

}  // namespace KTF

using namespace KTF ;

//*************************************************************************************************************************
//
//    KalTest Fitter 
//
//*************************************************************************************************************************



KalTestFitter aKalTestFitter ;


KalTestFitter::KalTestFitter() : Processor("KalTestFitter") {
  
  // modify processor description
  _description = "KalTestFitter : fit lcio::Track collection with KalTest " ;
  
  registerInputCollection( LCIO::TRACK,
			   "InputTrackCollection" , 
			   "Name of the input collections holding Tracks to be fitted"  ,
			   _inColName ,
			   "TPCTracks" ) ;
  
  registerOutputCollection( LCIO::TRACK,
			    "OutputCollection" , 
			    "Name of the output Track collections"  ,
			    _outColName ,
			    std::string("KalTestTracks" ) ) ;
}


void KalTestFitter::init() { 

  // usually a good idea to
  printParameters() ;

  _kalTest = new KalTest( *marlin::Global::GEAR ) ;

  _nRun = 0 ;
  _nEvt = 0 ;
}

void KalTestFitter::processRunHeader( LCRunHeader* run) { 
  _nRun++ ;
} 


//----------------------------------------------------------------------------------------------------------------------------------------------

void KalTestFitter::processEvent( LCEvent * evt ) { 

  clock_t start =  clock() ; 
  
  LCCollection* trkCol=0 ;

  try{  trkCol =  evt->getCollection( _inColName ) ;  } catch( lcio::DataNotAvailableException& e ) {
    
    streamlog_out( ERROR ) << " processEvent : can't find input collection : " << _inColName << std::endl ; 
    return ;
  }
  
  GenericClusterVec<TrackerHit> cluList ;
  

  // // create a vector of generic hits from the collection applying a cut on r_min
  // for( StringVec::iterator it = _colNames.begin() ; it !=  _colNames.end() ; it++ ){  
    
  //   LCCollectionVec* col =  dynamic_cast<LCCollectionVec*> (evt->getCollection( *it )  ); 
    
    
  //   //--- assign the layer number to the TrackerHits
    
  //   int nHit = col->getNumberOfElements() ;
  //   for(int i=0 ; i < nHit ; ++i ) {
      
  //     TrackerHitImpl* th = (TrackerHitImpl*) col->getElementAt(i) ;
  //     gear::Vector3D v( th->getPosition()[0],th->getPosition()[1], 0 ) ; 
  //     int padIndex = padLayout.getNearestPad( v.rho() , v.phi() ) ;
      
  //     th->ext<HitInfo>() = new HitInfoStruct ;

  //     th->ext<HitInfo>()->layerID = padLayout.getRowNumber( padIndex ) ;
      
  //     //       //--- for fixed sized rows this would also work...
  //     //       float rMin = padLayout.getPlaneExtent()[0] ;
  //     //       float rMax = padLayout.getPlaneExtent()[1] ;
  //     //       float nRow  = padLayout.getNRows() ;
  //     //       int lCheck =  ( v.rho() - rMin ) / ((rMax - rMin ) /nRow ) ;

  //     //       streamlog_out( DEBUG ) << " layerID : " << th->ext<HitInfo>()->layerID 
  //     // 			     << " r: " << v.rho() 
  //     // 			     << " lCheck : " << lCheck 
  //     // 			     << " phi : " << v.phi()
  //     // 			     << " rMin : " << rMin 
  //     // 			     << " rMax : " << rMax 
  //     // 			     << std::endl ;

  //   } //-------------------- end assign layernumber ---------
    
  //   //addToGenericHitVec( h , col , rCut , zIndex ) ;
  //   std::list< TrackerHit*> hitList ;
  //   TrackerHitCast cast ;
  //   ZSort zsort ;
  //   std::transform(  col->begin(), col->end(), std::back_inserter( hitList ), cast ) ;
  //   hitList.sort( zsort ) ;
  //   //    std::for_each( hitList.begin() , hitList.end() , printZ ) ;

  //   addToGenericHitVec( h, hitList.begin() , hitList.end() , rCut ,  zIndex ) ;
  // }  
  
  // // cluster the sorted hits  ( if |diff(z_index)|>1 the loop is stopped)
  // cluster_sorted( h.begin() , h.end() , std::back_inserter( cluList )  , &dist0 , _minCluSize ) ;
  // //cluster( h.begin() , h.end() , std::back_inserter( cluList )  , &dist , _minCluSize ) ;
  
  // streamlog_out( DEBUG ) << "   ***** clusters: " << cluList.size() << std::endl ; 

  // LCCollectionVec* allClu = new LCCollectionVec( LCIO::TRACK ) ;
  // std::transform(cluList.begin(), cluList.end(), std::back_inserter( *allClu ) , converter ) ;
  // evt->addCollection( allClu , "AllTrackClusters" ) ;



  // // find 'odd' clusters that have duplicate hits in pad rows
  // GenericClusterVec<TrackerHit> ocs ;


  // //  typedef GenericClusterVec<TrackerHit>::iterator GCVI ;

  // //   for( GCVI it = cluList.begin() ; it != cluList.end() ; ++it ){
  // //     std::cout << " *** cluster :" << (*it)->ID 
  // // 	      << " size() :" << (*it)->size() 
  // // 	      << " at : " << *it << std::endl ;
  // //   }
  // //  GCVI remIt =
  // //     std::remove_copy_if( cluList.begin(), cluList.end(), std::back_inserter( ocs ) ,  DuplicatePadRows() ) ;

  // split_list( cluList, std::back_inserter(ocs),  DuplicatePadRows( nPadRows, _duplicatePadRowFraction  ) ) ;

  // //   for( GCVI it = cluList.begin() ; it != cluList.end() ; ++it ){
  // //     std::cout << " +++ cluster :" << (*it)->ID 
  // // 	      << " size() :" << (*it)->size() 
  // // 	      << " at : " << *it << std::endl ;
  // //   }

  // //   for( GCVI it = cluList.begin() ; it != cluList.end() ; ++it ){
  // //     if( (*it)->size() > 20 ){
  // //       ocs.splice( ocs.begin() , cluList , it  );
  // //     }
  // //   }


  // LCCollectionVec* oddCol = new LCCollectionVec( LCIO::TRACK ) ;
  // std::transform( ocs.begin(), ocs.end(), std::back_inserter( *oddCol ) , converter ) ;
  // evt->addCollection( oddCol , "OddClu_1" ) ;


  // streamlog_out( DEBUG ) << "   ***** clusters: " << cluList.size() 
  // 			 << "   ****** oddClusters " << ocs.size() 
  // 			 << std::endl ; 



  // //   //-------------------- split up cluster with duplicate rows 

  // GenericClusterVec<TrackerHit> sclu ; // new split clusters

  // std::vector< GenericHit<TrackerHit>* > oddHits ;
  // oddHits.reserve( h.size() ) ;

  // typedef GenericClusterVec<TrackerHit>::iterator GCVI ;


  // //========================== first iteration ================================================
  // for( GCVI it = ocs.begin() ; it != ocs.end() ; ++it ){
  //   (*it)->takeHits( std::back_inserter( oddHits )  ) ;
  //   delete (*it) ;
  // }
  // ocs.clear() ;

  // int _nRowForSplitting = 10 ; //FIXME:  make proc param
  // // reset the hits index to row ranges for reclustering
  // unsigned nOddHits = oddHits.size() ;
  // for(unsigned i=0 ; i< nOddHits ; ++i){
  //   int layer =  oddHits[i]->first->ext<HitInfo>()->layerID  ;
  //   oddHits[i]->Index0 =   2 * int( layer / _nRowForSplitting ) ;
  // }

  // //----- recluster in pad row ranges
  // cluster( oddHits.begin(), oddHits.end() , std::back_inserter( sclu ), &dist , _minCluSize ) ;

  // LCCollectionVec* oddCol2 = new LCCollectionVec( LCIO::TRACK ) ;
  // std::transform( sclu.begin(), sclu.end(), std::back_inserter( *oddCol2 ) , converter ) ;
  // evt->addCollection( oddCol2 , "OddClu_2" ) ;


  // streamlog_out( DEBUG ) << "   ****** oddClusters fixed" << sclu.size() 
  // 			 << std::endl ; 
 
  // //--------- remove pad row range clusters where merge occured 
  // split_list( sclu, std::back_inserter(ocs), DuplicatePadRows( nPadRows, _duplicatePadRowFraction  ) ) ;


  // LCCollectionVec* oddCol3 = new LCCollectionVec( LCIO::TRACK ) ;
  // std::transform( ocs.begin(), ocs.end(), std::back_inserter( *oddCol3 ) , converter ) ;
  // evt->addCollection( oddCol3 , "OddClu_3" ) ;


  

  // for( GCVI it = ocs.begin() ; it != ocs.end() ; ++it ){
  //   (*it)->takeHits( std::back_inserter( oddHits )  ) ;
  //   delete (*it) ;
  // }
  // ocs.clear() ;


  // //   //========================== second iteration in shifted pad row ranges ================================================


  // oddHits.clear() ;
  // for( GCVI it = sclu.begin() ; it != sclu.end() ; ++it ){
  //   (*it)->takeHits( std::back_inserter( oddHits )  ) ;
  //   delete (*it) ;
  // }
  // sclu.clear() ;

  // //  int _nRowForSplitting = 10 ; //FIXME:  make proc param
  // // reset the hits index to row ranges for reclustering
  // nOddHits = oddHits.size() ;

  // streamlog_out( DEBUG ) << "   left over odd hits for second iteration of pad row range clustering " << nOddHits << std::endl ;

  // for(unsigned i=0 ; i< nOddHits ; ++i){
  //   int layer =  oddHits[i]->first->ext<HitInfo>()->layerID  ;
  //   oddHits[i]->Index0 =   2 * int( 0.5 +  ( (float) layer / (float) _nRowForSplitting ) ) ;
  // }

  // //----- recluster in pad row ranges
  // cluster( oddHits.begin(), oddHits.end() , std::back_inserter( sclu ), &dist , _minCluSize ) ;

  // LCCollectionVec* oddCol2_1 = new LCCollectionVec( LCIO::TRACK ) ;
  // std::transform( sclu.begin(), sclu.end(), std::back_inserter( *oddCol2_1 ) , converter ) ;
  // evt->addCollection( oddCol2_1 , "OddClu_2_1" ) ;


  // streamlog_out( DEBUG ) << "   ****** oddClusters fixed" << sclu.size() 
  // 			 << std::endl ; 
 
  // //--------- remove pad row range clusters where merge occured 
  // split_list( sclu, std::back_inserter(ocs), DuplicatePadRows( nPadRows, _duplicatePadRowFraction  ) ) ;


  // LCCollectionVec* oddCol3_1 = new LCCollectionVec( LCIO::TRACK ) ;
  // std::transform( ocs.begin(), ocs.end(), std::back_inserter( *oddCol3_1 ) , converter ) ;
  // evt->addCollection( oddCol3_1 , "OddClu_3_1" ) ;

  // //----------------end  split up cluster with duplicate rows 
  
  // for( GCVI it = ocs.begin() ; it != ocs.end() ; ++it ){
  //   (*it)->takeHits( std::back_inserter( oddHits )  ) ;
  //   delete (*it) ;
  // }
  // ocs.clear() ;
  

  // //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  // // --- recluster the good clusters w/ all pad rows

  // oddHits.clear() ;
  // for( GCVI it = sclu.begin() ; it != sclu.end() ; ++it ){
  //   (*it)->takeHits( std::back_inserter( oddHits )  ) ;
  //   delete (*it) ;
  // }
  // sclu.clear() ;

  // //   reset the index for 'good' hits coordinate again...
  // nOddHits = oddHits.size() ;
  // for(unsigned i=0 ; i< nOddHits ; ++i){
  //   oddHits[i]->Index0 = zIndex ( oddHits[i]->first ) ;
  // }

  // cluster( oddHits.begin(), oddHits.end() , std::back_inserter( sclu ), &dist , _minCluSize ) ;

  // LCCollectionVec* oddCol4 = new LCCollectionVec( LCIO::TRACK ) ;
  // std::transform( sclu.begin(), sclu.end(), std::back_inserter( *oddCol4 ) , converter ) ;
  // evt->addCollection( oddCol4 , "OddClu_4" ) ;

  // // --- end recluster the good clusters w/ all pad rows

  // // merge the good clusters to final list
  // cluList.merge( sclu ) ;
  
  // LCCollectionVec* cluCol = new LCCollectionVec( LCIO::TRACK ) ;
  // std::transform( cluList.begin(), cluList.end(), std::back_inserter( *cluCol ) , converter ) ;
  // evt->addCollection( cluCol , "CluTrackSegments" ) ;


  // //================================================================================
   
  // // create vecor with left over hits in each layer
  // std::vector< Hit* > leftOverHits ;
  // leftOverHits.reserve(  h.size() ) ;

  // typedef GenericHitVec<TrackerHit>::const_iterator GHVI ;

  // for( GHVI it = h.begin(); it != h.end() ;++it){

  //   if( (*it)->second == 0 ){

  //     leftOverHits.push_back( *it ) ;
  //   }
  // }

  // //============ now try reclustering with circle segments ======================
  // GenericClusterVec<TrackerHit> mergedClusters ; // new split clusters

  // std::list< ClusterSegment* > segs ;
  
  // // fit circles (helices actually) to the cluster segments...
  // std::transform( cluList.begin(), cluList.end(), std::back_inserter( segs ) , CircleFitter() ) ;


  // //*********************************************************
  // //   try  to run KalTest on circle segments ------------------------------------------------------
  // //*********************************************************
  // //  LCIOTrackFromSegment segToTrack ;

  // streamlog_out( DEBUG ) <<  "************* fitted segments and KalTest tracks : **********************************" 
  // 			 << std::endl ;



  LCCollectionVec* kaltracks = new LCCollectionVec( LCIO::TRACK ) ;

  LCFlagImpl trkFlag(0) ;
  trkFlag.setBit( LCIO::TRBIT_HITS ) ;
  kaltracks->setFlag( trkFlag.getFlag()  ) ;

  evt->addCollection( kaltracks , _outColName ) ;


  std::list< KalTrack* > ktracks ;

  KTFitter<KalTest::OrderOutgoing, KalTest::FitBackward, KalTest::PropagateToIP > ipFitter( _kalTest ) ;
  
  std::transform( cluList.begin(), cluList.end(), std::back_inserter( ktracks ) , ipFitter  ) ;

  std::transform( ktracks.begin(), ktracks.end(), std::back_inserter( *kaltracks ) , KalTrack2LCIO() ) ;


  //========== cleanup KalTracks ========


  std::for_each( ktracks.begin() , ktracks.end() , delete_ptr<KalTrack> ) ;

  //=====================================

  _nEvt ++ ;

  clock_t end = clock () ; 

  streamlog_out( DEBUG )  << "---  clustering time: " 
 			  <<  double( end - start ) / double(CLOCKS_PER_SEC) << std::endl  ;

}


/*************************************************************************************************/
void KalTestFitter::check( LCEvent * evt ) { 
  /*************************************************************************************************/

  //  UTIL::LCTOOLS::dumpEventDetailed( evt ) ;

  // bool checkForDuplicatePadRows =  true ;
  // bool checkForMCTruth =  true ;

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
  //   LCIterator<Track> trIt( evt, _outColName ) ;
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
  // 	if( hitsInLayer[i] > 1 ) 
  // 	  ++nDouble ;
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
    
  //   LCIterator<Track> trIt( evt, "KalTestTracks" ) ;
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

  //     if( double(maxHit) / nHit < 0.99 ){ // What is acceptable here ???
  // 	//if( nDouble  > 0){
  // 	streamlog_out( MESSAGE ) << " oddTrackCluster found with only "
  // 				 << 100.*double(maxHit)/nHit 
  // 				 << "% of hits  form one MCParticle " << std::endl ;
  // 	oddCol->addElement( tr ) ;
  //     }
  //   }
  //   evt->addCollection( oddCol , "OddMCPTracks" ) ;

  //   streamlog_out( DEBUG ) << " checking for split tracks - mcptrkmap size : " <<  mcpTrkMap.size() << std::endl ;

  //   // check for split tracks 
  //   for( MCPTRKMAP::iterator it0 = mcpTrkMap.begin() ; it0 != mcpTrkMap.end() ; ++it0){
      
  //     streamlog_out( DEBUG ) << " checking for split tracks - map size : " <<  it0->second.size() << std::endl ;


  //     if( it0->second.size() > 1 ) {


  // 	typedef std::list< EVENT::Track* > TL ;
  // 	TL trkList ;

  // 	for( TRKMAP::iterator it1 = it0->second.begin() ; it1 != it0->second.end() ; ++it1){
	    
  // 	  double totalHits = hitMap[ it0->first ]  ; // total hits for this track 

  // 	  double thisMCPHits = it1->second ;     //  hits from this mcp

  // 	  double ratio =  thisMCPHits / totalHits  ;

  // 	  streamlog_out( DEBUG ) << " checking for split tracks - ratio : " 
  // 				 << totalHits << " / " << thisMCPHits << " = " << ratio << std::endl ;

  // 	  if( ratio > 0.03 && ratio < 0.95 ){
  // 	    // split track

  // 	    splitCol->addElement( it1->first ) ; 

  // 	    trkList.push_back( it1->first ) ;
  // 	  } 
  // 	}

  // 	streamlog_out( DEBUG2 ) << " ------------------------------------------------------ " << std::endl ;

  // 	for( TL::iterator it0 = trkList.begin() ; it0 != trkList.end() ; ++it0 ){
	  

  // 	  KalTrack* trk0 = (*it0)->ext<KalTrackLink>() ; 
	  
  // 	  HelixClass hel ;
  // 	  hel.Initialize_Canonical( (*it0)->getPhi(),
  // 				    (*it0)->getD0(),
  // 				    (*it0)->getZ0(),
  // 				    (*it0)->getOmega(),
  // 				    (*it0)->getTanLambda(),
  // 				    3.50 ) ;

  // 	  streamlog_out( DEBUG1 ) << hel.getXC() << "\t"
  // 				  << hel.getYC() << "\t"
  // 				  << hel.getRadius() << "\t" 
  // 				  << hel.getTanLambda() << std::endl ; 
	  

  // 	  // streamlog_out( DEBUG1 ) << (*it0)->getPhi() << "\t"
  // 	  // 			  << (*it0)->getD0()  << "\t"
  // 	  // 			  << (*it0)->getOmega()  << "\t"
  // 	  // 			  << (*it0)->getZ0()  << "\t"
  // 	  // 			  << (*it0)->getTanLambda()  << "\t"
  // 	  // 			  << std::endl ;

  // 	  for( TL::iterator it1 =  trkList.begin() ; it1 != trkList.end() ; ++it1 ){
	    
  // 	    KalTrack* trk1 = (*it1)->ext<KalTrackLink>() ; 
	    
  // 	    double chi2 =  KalTrack::chi2( *trk0 ,  *trk1 ) ;

  // 	    // streamlog_out( DEBUG1 ) << "  chi2 between split tracks : " 
  // 	    // 			    << trk0 << " - " << trk1 << " : " << chi2 << std::endl ; 

	    
  // 	  }
  // 	}

  //     }
  //   }
  //   evt->addCollection( splitCol , "SplitTracks" ) ;
    

  // }
  //====================================================================================

}


void KalTestFitter::end(){ 
  
  //   std::cout << "KalTestFitter::end()  " << name() 
  // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
  // 	    << std::endl ;
  
  delete _kalTest ;
}

