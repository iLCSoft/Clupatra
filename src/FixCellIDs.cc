#include "FixCellIDs.h"
#include <iostream>

#include <EVENT/LCCollection.h>
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
  _description = "fix CellID0 for old TrackerHits ..." ;
  
  
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
  //      VXD 
  //=====================================================================================================


  LCCollection* vxdCol = 0 ;
  
  try{   vxdCol =  dynamic_cast<LCCollection*> (evt->getCollection( _vxdColName )  ); 
    
  } catch( lcio::DataNotAvailableException& e) { 
    
    streamlog_out( DEBUG5 ) <<  " input collection not in event : " << _vxdColName << "   - nothing to do for VXD hits  !!! " << std::endl ;  
  } 
  
  if( vxdCol ) { 
    
    const gear::ZPlanarParameters& gearVXD = Global::GEAR->getVXDParameters() ;
    
    int nHit = vxdCol->getNumberOfElements()  ;
    
    for(int i=0; i< nHit ; i++){
      
      TrackerHitImpl* h = dynamic_cast<TrackerHitImpl*>( vxdCol->getElementAt( i ) ) ;
      
      if( !h )  continue ;

      gear::Vector3D   pos(  h->getPosition() ) ;
            
      //-----------------------------------
      gear::SensorID sensorID ;
      
      bool hitOnLadder =  gearVXD.isPointInSensitive( pos ,  &sensorID ) ;

      if( !hitOnLadder ) {
	
	streamlog_out( DEBUG5 ) <<  " hit not in sensitive volume  - distance: " << gearVXD.distanceToNearestSensitive( pos ).r() 
				<<  " pos : " << pos << std::endl ;
      }
     //-------------------------------------
      
      encoder.reset() ;  // reset to 0
      
      encoder[ILDCellID0::subdet] = ILDDetID::VXD   ;
      encoder[ILDCellID0::side]   = 0 ; //sensorID.side   ; 
      encoder[ILDCellID0::layer]  = sensorID.layer  ;
      encoder[ILDCellID0::module] = sensorID.module ;
      encoder[ILDCellID0::sensor] = sensorID.sensor ;
      
      int layerID = encoder.lowWord() ;  
      
      h->setCellID0(  layerID ) ;
      
    }
    vxdCol->parameters().setValue( "CellIDEncoding" , ILDCellID0::encoder_string ) ;

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


