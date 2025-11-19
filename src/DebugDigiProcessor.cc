#include "DebugDigiProcessor.h"

//#include <time.h>
// #include <vector>
// #include <map>
// #include <algorithm>
// #include <math.h>
// #include <cmath>

//---- LCIO ---
#include "IMPL/LCCollectionVec.h"
#include "EVENT/Track.h"
#include "EVENT/TrackerHit.h"
#include "UTIL/Operators.h"
// #include "UTIL/LCRelationNavigator.h"
// #include "UTIL/LCTypedVector.h"
#include "UTIL/LCTrackerConf.h"
#include <UTIL/ILDConf.h>
#include "UTIL/CellIDDecoder.h"
// #include "UTIL/LCTOOLS.h"


#include "UTIL/LCIterator.h"

using namespace lcio ;
using namespace marlin ;


//--------------------------------------------------------------------------------------------

template <class T>
int subdet( const T* t){
  
  static CellIDDecoder<T> idDec(  LCTrackerCellID::encoding_string() ) ;
  
  return idDec( t )[ LCTrackerCellID::subdet() ] ;  
}

//--------------------------------------------------------------------------------------------




DebugDigiProcessor aDebugDigiProcessor ;


DebugDigiProcessor::DebugDigiProcessor() : Processor("DebugDigiProcessor") ,
					   _nRun(0), _nEvt(0) {
  
  // modify processor description
  _description = "DebugDigiProcessor: analysis plots for Tracks" ;
  
  
  registerInputCollection( LCIO::TRACK,
			   "TrackCollection" , 
			   "Name of the input collections with  Tracks"  ,
			   _trkColName ,
			   std::string("TrackCollection" ) ) ;
  
  IntVec subDetIDs ;
  subDetIDs.push_back( lcio::ILDDetID::TPC ) ;
  
  registerProcessorParameter( "SubDetectorIDs" , 
			      "The IDs of the subdetectors corresponding to the hit collections in TrackerHitCollectionNames"  ,
			      _subDetIDs ,
			      subDetIDs ) ;
  
  StringVec hitColNames ;
  hitColNames.push_back( std::string( "TPCTrackerHits" ) ) ;
  
  registerProcessorParameter( "TrackerHitCollections" , 
			      "The tracker hit collection names - needs to run parallel to SubDetectorIDs"  ,
			      _hitColNames ,
			      hitColNames ) ;
  



}


void DebugDigiProcessor::init() { 

  // usually a good idea to
  printParameters() ;


 
  if( _subDetIDs.size() != _hitColNames.size() ){

    throw Exception( "DebugDigiProcessor::init() - processor parameters SubDetectorIDs and TrackerHitCollections do not have the same size !" ) ;
  }


  _nRun = 0 ;
  _nEvt = 0 ;
}

void DebugDigiProcessor::processRunHeader( LCRunHeader* ) { 

  _nRun++ ;
} 



void DebugDigiProcessor::modifyEvent( LCEvent * evt ) { 

  for( unsigned i=0,N=_hitColNames.size() ; i<N ; ++i){
    
    LCCollection* col = 0 ;
    try { 
      
      col = evt->getCollection( _hitColNames[ i ]  ) ;
    }
    catch( const DataNotAvailableException& ) {}

    if( col ) {

      evt->removeCollection(  _hitColNames[ i ] ) ;
      
      std::string newName( _hitColNames[ i ] ) ;
      newName += "_OLD" ;

      evt->addCollection( col , newName ) ;

      streamlog_out( DEBUG2 ) << " rename collection " <<  _hitColNames[ i ] << " to  " <<  newName << std::endl ;

    }
  }
  
  LCCollection* col = 0 ;
  try { 
      
    col = evt->getCollection( _trkColName ) ;
  }
  catch( const DataNotAvailableException& ) {}

  if( col ) {
    
    evt->removeCollection(  _trkColName ) ;
      
      std::string newName( _trkColName ) ;
      newName += "_OLD" ;

      evt->addCollection( col , newName ) ;

      streamlog_out( DEBUG2 ) << " rename collection " <<  _trkColName << " to  " <<  newName << std::endl ;

    }
  


}



void DebugDigiProcessor::processEvent( LCEvent * evt ) { 


  //=========================================================================================================================
  
  
  LCCollection* trkCol = 0 ;
  try{
    
    std::string newName( _trkColName ) ;
    newName += "_OLD" ;
    trkCol = evt->getCollection( newName ) ;
  }
  catch( lcio::DataNotAvailableException& e){
    
    streamlog_out( DEBUG2 ) << " *** collection not in event : " << _trkColName << std::endl ;
    return ; // nothing to do in this event (no tracks)
  }

  if( trkCol->getNumberOfElements() < 1 ) {
    
    streamlog_out( DEBUG2) << " not tracks in collection " << _trkColName << " evt : " << evt->getEventNumber()
			   << " run : " << evt->getRunNumber() << std::endl ;
    return ;
  }


  std::map< int , LCCollectionVec* >  colMap ;

  for( unsigned i=0,N=_hitColNames.size() ; i<N ; ++i){

    LCCollectionVec* col = new LCCollectionVec( LCIO::TRACKERHIT )  ;

    col->setSubset( true ) ;

    evt->addCollection(  col , _hitColNames[ i ] ) ;
    
    colMap[ _subDetIDs[i] ] = col ;

    streamlog_out( DEBUG2) << " add collection " <<  _hitColNames[ i ]  << " to event " << std::endl ;
  }


  for( unsigned i=0,N=trkCol->getNumberOfElements() ; i<N ; ++i){
  
    Track* trk = (Track*) ( trkCol->getElementAt(i) ) ;

    const TrackerHitVec& hits = trk->getTrackerHits() ;

    for( unsigned j=0,M=hits.size() ; j<M ; ++j){
      
      int subDetID = subdet( hits[j] )  ;

      LCCollectionVec* col =  colMap[ subDetID ]  ;

      if(  col != 0  )  {

	streamlog_out( DEBUG ) << " add hit " << hits[j] << " to collection for sub detector : " << subDetID << std::endl ;

	col->addElement( hits[j] ) ;
	
      }
    }
  }

  //----------------------------------------------------------------------  
  
  _nEvt ++ ;

}




/*************************************************************************************************/
void DebugDigiProcessor::check( LCEvent *  ) { 
  /*************************************************************************************************/

}
//====================================================================================




void DebugDigiProcessor::end(){ 
  
  streamlog_out( MESSAGE )  << " processed " << _nEvt << " events in " << _nRun << " runs "
			    << std::endl ;
}

//====================================================================================================

