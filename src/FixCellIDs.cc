#include "FixCellIDs.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <Impl/TrackerHitImpl.h>
#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>
#include <UTIL/LCTOOLS.h>

//---- GEAR ----
#include "marlin/Global.h"
#include "gear/GEAR.h"
#include "gear/TPCParameters.h"
#include "gear/PadRowLayout2D.h"
#include "gear/BField.h"

using namespace lcio ;
using namespace marlin ;


FixCellIDs aFixCellIDs ;


FixCellIDs::FixCellIDs() : Processor("FixCellIDs") {
  
  // modify processor description
  _description = "fix CellID0 for old TrackerHits ..." ;
  
  
  // register steering parameters: name, description, class-variable, default value
  registerInputCollection( LCIO::TRACKERHIT,
			   "CollectionName" , 
			   "Name of the TrackerHit collection"  ,
			   _colName ,
			   std::string("AllTPCTrackerHits")
			   );
}


void FixCellIDs::init() { 
  
  streamlog_out(DEBUG) << "   init called  " << std::endl ;
  
  // usually a good idea to
  printParameters() ;
  
  _nRun = 0 ;
  _nEvt = 0 ;
  
}


void FixCellIDs::processRunHeader( LCRunHeader* run) { 
  
  _nRun++ ;
} 



void FixCellIDs::modifyEvent( LCEvent * evt ) { 

  const gear::TPCParameters& gearTPC = Global::GEAR->getTPCParameters() ;
  const gear::PadRowLayout2D& padLayout = gearTPC.getPadLayout() ;


  LCCollection* col = 0 ;
  
  try{   col =  dynamic_cast<LCCollection*> (evt->getCollection( _colName )  ); 
    
  } catch( lcio::DataNotAvailableException& e) { 
    
    streamlog_out( WARNING ) <<  " input collection not in event : " << _colName << "   - nothing to do  !!! " << std::endl ;  
    return ;
  } 
  
  
  int nHit = col->getNumberOfElements()  ;
  
  UTIL::BitField64 encoder( ILDCellID0::encoder_string ) ; 
  
  for(int i=0; i< nHit ; i++){
    
    TrackerHitImpl* h = dynamic_cast<TrackerHitImpl*>( col->getElementAt( i ) ) ;
    

    gear::Vector3D   pos(  h->getPosition() ) ;
    
    int padIndex = padLayout.getNearestPad( pos.rho() , pos.phi() ) ;
    
    int layer = padLayout.getRowNumber( padIndex ) ;


    encoder.reset() ;  // reset to 0
    
    encoder[ILDCellID0::subdet] = ILDDetID::TPC ;
    encoder[ILDCellID0::side]   = 0 ;
    encoder[ILDCellID0::layer]  = layer ;
    encoder[ILDCellID0::module] = 0 ;
    encoder[ILDCellID0::sensor] = 0 ;
    
    int layerID = encoder.lowWord() ;  

    h->setCellID0(  layerID ) ;
    
    col->parameters().setValue( "CellIDEncoding" , ILDCellID0::encoder_string ) ;
  
  } 
  


  //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !
  
  streamlog_out( DEBUG3 ) << "   processing event: " << evt->getEventNumber() 
			  << "   in run:  " << evt->getRunNumber() << std::endl ;


  // if( streamlog_level( DEBUG )  ) {
  //   UTIL::LCTOOLS::printTrackerHits(  col ) ;
  // }
  
  _nEvt ++ ;
}



void FixCellIDs::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void FixCellIDs::end(){ 

  streamlog_out( MESSAGE ) << "FixCellIDs::end()  " << name() 
			   << " processed " << _nEvt << " events in " << _nRun << " runs "
			   << std::endl ;
  
}


