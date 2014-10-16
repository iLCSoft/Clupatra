#include "FixCellIDs.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <IMPL/SimTrackerHitImpl.h>
#include <IMPL/TrackerHitImpl.h>
#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>
#include <UTIL/LCTOOLS.h>

//---- GEAR ----
#include "marlin/Global.h"
#include "gear/GEAR.h"
#include "gear/TPCParameters.h"
#include "gear/ZPlanarParameters.h"
#include "gear/ZPlanarLayerLayout.h"
#include "gear/PadRowLayout2D.h"
#include "gear/BField.h"

using namespace lcio ;
using namespace marlin ;


FixCellIDs aFixCellIDs ;


FixCellIDs::FixCellIDs() : Processor("FixCellIDs") ,
			   _nRun(0), _nEvt(0) {
  
  // modify processor description
  _description = "fix CellID0  TrackerHits  based on the position and the GEAR information ..." ;
  
  
  // register steering parameters: name, description, class-variable, default value
  registerInputCollection( LCIO::TRACKERHIT,
			   "TPCCollectionName" , 
			   "Name of the TPC TrackerHit collection"  ,
			   _tpcColName ,
			   std::string("AllTPCTrackerHits")
			   );
  
  registerInputCollection( LCIO::TRACKERHIT,
			   "VXDCollectionName" , 
			   "Name of the VXD TrackerHit collection"  ,
			   _vxdColName ,
			   std::string("VXDTrackerHits")
			   );
  
  registerInputCollection( LCIO::TRACKERHIT,
			   "SITCollectionName" , 
			   "Name of the SIT TrackerHit collection"  ,
			   _sitColName ,
			   std::string("SITTrackerHits")
			   );
}


bool position( lcio::LCObject* o, gear::Vector3D& pos ){
  
  TrackerHit* th = dynamic_cast<TrackerHit*>( o ) ;
  if( th ) {
    pos = th->getPosition() ;
  }
  SimTrackerHit* sth = dynamic_cast<SimTrackerHit*>( o ) ;
  if( sth ) {
    pos = sth->getPosition() ;
  }

  return ( th != 0 || sth != 0 ) ;
}

void setCellID(  lcio::LCObject* o, const UTIL::BitField64& enc ){

  TrackerHitImpl* th = dynamic_cast<TrackerHitImpl*>( o ) ;
  if( th ) th->setCellID0(   enc.lowWord() ) ;
  
  
  SimTrackerHitImpl* sth = dynamic_cast<SimTrackerHitImpl*>( o ) ;
  if( sth ) sth->setCellID0(   enc.lowWord() ) ;

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



void fixZPlanarHits(LCEvent* evt, const std::string& colName, int detID , const gear::ZPlanarParameters& gearZPLAN ){
  
  LCCollection* zplanCol = 0 ;

  UTIL::BitField64 encoder( ILDCellID0::encoder_string ) ; 
  
  try{   zplanCol =  dynamic_cast<LCCollection*> (evt->getCollection( colName )  ); 
    
  } catch( lcio::DataNotAvailableException& e) { 
    
    streamlog_out( DEBUG5 ) <<  " input collection not in event : " << colName << "   - nothing to do for ZPLAN hits  !!! " << std::endl ;  
  } 
  
  if( zplanCol ) { 
    
    
 
   int nHit = zplanCol->getNumberOfElements()  ;
    
    for(int i=0; i< nHit ; i++){
      
      // TrackerHitImpl* h = dynamic_cast<TrackerHitImpl*>( zplanCol->getElementAt( i ) ) ;
      // if( !h )  continue ;
      // gear::Vector3D   pos(  h->getPosition() ) ;

      lcio::LCObject* h =  zplanCol->getElementAt( i ) ;
      gear::Vector3D   pos ;
      if( ! position( h, pos ) ) continue ;
		 
      //-----------------------------------
      gear::SensorID sensorID ;
      
      bool hitOnLadder =  gearZPLAN.isPointInSensitive( pos ,  &sensorID ) ;

      if( !hitOnLadder ) {
	
	streamlog_out( DEBUG5 ) <<  " hit not in sensitive volume  - distance: " << gearZPLAN.distanceToNearestSensitive( pos ).r() 
				<<  " pos : " << pos << std::endl ;
      }
     //-------------------------------------
      
      encoder.reset() ;  // reset to 0
      
      encoder[ILDCellID0::subdet] = detID  ;
      encoder[ILDCellID0::side]   = 0 ; //sensorID.side   ; 
      encoder[ILDCellID0::layer]  = sensorID.layer  ;
      encoder[ILDCellID0::module] = sensorID.module ;
      encoder[ILDCellID0::sensor] = sensorID.sensor ;
      
      //int layerID = encoder.lowWord() ;  
      // h->setCellID0(  layerID ) ;

      setCellID( h, encoder ) ;
      
    }
    zplanCol->parameters().setValue( "CellIDEncoding" , ILDCellID0::encoder_string ) ;

  }
}




void FixCellIDs::modifyEvent( LCEvent * evt ) { 


   UTIL::BitField64 encoder( ILDCellID0::encoder_string ) ; 


  //=====================================================================================================
  //      TPC 
  //=====================================================================================================
  
  
  LCCollection* tpcCol = 0 ;

  try{   tpcCol =  dynamic_cast<LCCollection*> (evt->getCollection( _tpcColName )  ); 

  } catch( lcio::DataNotAvailableException& e) { 
    
    streamlog_out( DEBUG5 ) <<  " input collection not in event : " << _tpcColName << "   - nothing to do for TPC hits  !!! " << std::endl ;  
  } 
  
  if( tpcCol ) {

    const gear::TPCParameters& gearTPC = Global::GEAR->getTPCParameters() ;
    const gear::PadRowLayout2D& padLayout = gearTPC.getPadLayout() ;
    
    int nHit = tpcCol->getNumberOfElements()  ;
    
    for(int i=0; i< nHit ; i++){
      
      TrackerHitImpl* h = dynamic_cast<TrackerHitImpl*>( tpcCol->getElementAt( i ) ) ;
      
      
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
      
      tpcCol->parameters().setValue( "CellIDEncoding" , ILDCellID0::encoder_string ) ;
      
    } 
  }

  //=====================================================================================================
  //      ZPlanar detectors 
  //=====================================================================================================


  try{
    const gear::ZPlanarParameters& gearVXD = Global::GEAR->getVXDParameters() ;

    fixZPlanarHits( evt, _vxdColName , ILDDetID::VXD , gearVXD ) ;
   
  }catch(gear::UnknownParameterException& ){
    streamlog_out( WARNING ) << " FixCellIDs : no gear parameters found for VXD " << std::endl ;
  }

  try{
    const gear::ZPlanarParameters& gearSIT = Global::GEAR->getSITParameters() ;

    fixZPlanarHits( evt, _sitColName , ILDDetID::SIT , gearSIT ) ;
   
  }catch(gear::UnknownParameterException& ){
    streamlog_out( WARNING ) << " FixCellIDs : no gear parameters found for SIT " << std::endl ;
  }


  
  streamlog_out( DEBUG3 ) << "   processing event: " << evt->getEventNumber() 
			  << "   in run:  " << evt->getRunNumber() << std::endl ;




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


