#include "FixCellIDs.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <IMPL/SimTrackerHitImpl.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackerHitPlaneImpl.h>

#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>
#include <UTIL/LCTOOLS.h>
#include <UTIL/Operators.h>

//---- GEAR ----
#include "marlin/Global.h"
#include "gear/GEAR.h"
#include "gear/TPCParameters.h"
#include "gear/ZPlanarParameters.h"
#include "gear/ZPlanarLayerLayout.h"
#include "gearimpl/FTDParametersImpl.h"
#include "gear/FTDLayerLayout.h"
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
			   std::string("NONE")
			   );
  
  registerInputCollection( LCIO::TRACKERHIT,
			   "VXDCollectionName" , 
			   "Name of the VXD TrackerHit collection"  ,
			   _vxdColName ,
			   std::string("NONE")
			   );
  
  registerInputCollection( LCIO::TRACKERHIT,
			   "SITCollectionName" , 
			   "Name of the SIT TrackerHit collection"  ,
			   _sitColName ,
			   std::string("NONE")
			   );

  registerInputCollection( LCIO::TRACKERHIT,
			   "SETCollectionName" , 
			   "Name of the SET TrackerHit collection"  ,
			   _setColName ,
			   std::string("NONE")
			   );
  
  registerInputCollection( LCIO::TRACKERHIT,
			   "FTDCollectionName" , 
			   "Name of the FTD TrackerHit collection"  ,
			   _ftdColName ,
			   std::string("NONE")
			   );

  registerProcessorParameter("SetSideForVXD",
			    "if true, the side is decoded as +-1 in the cellID for the VXD (zplanar hits)",
			    _setSideForVXD,
			    bool( false ) ) ;


  registerProcessorParameter("LayerOffset",
			     " add  LayerOffset to layer index in cellID",
			     _layerOffset,
			     int(0) ) ;
  

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

  SimTrackerHitImpl* sth = dynamic_cast<SimTrackerHitImpl*>( o ) ;
  if( sth ){
    sth->setCellID0(   enc.lowWord() ) ;
    return ;
  }
  TrackerHitImpl* th = dynamic_cast<TrackerHitImpl*>( o ) ;
  if( th ){
    th->setCellID0(   enc.lowWord() ) ;
    return ;
  }
  
  TrackerHitPlaneImpl* thp = dynamic_cast<TrackerHitPlaneImpl*>( o ) ;
  if( thp ){
    thp->setCellID0(   enc.lowWord() ) ;
    return ;
  }
  
  throw Exception("setCellID called with type that is neither of: SimTrackerHitImpl, TrackerHitImpl, TrackerHitPlaneImpl " ) ;
}

unsigned cellID0(  lcio::LCObject* o ){
  
  SimTrackerHit* sth = dynamic_cast<SimTrackerHit*>( o ) ;
  if( sth ) return sth->getCellID0() ;

  TrackerHit* th = dynamic_cast<TrackerHit*>( o ) ;
  if( th )  return th->getCellID0() ;

  // TrackerHitPlaneImpl* thp = dynamic_cast<TrackerHitPlaneImpl*>( o ) ;
  // if( thp )  return thp->getCellID0() ;

  return -1 ;
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



void fixZPlanarHits(LCEvent* evt, const std::string& colName, int detID , const gear::ZPlanarParameters& gearZPLAN , double layerOffset, bool setSideZPlanar=false){
  
  LCCollection* zplanCol = 0 ;

  UTIL::BitField64 encoder( ILDCellID0::encoder_string ) ; 
  
  try{   

    zplanCol =  dynamic_cast<LCCollection*> (evt->getCollection( colName )  ); 

    streamlog_out( DEBUG3 ) <<  " process input collection : " << colName  << std::endl ;  

  } catch( lcio::DataNotAvailableException& e) { 
    
    streamlog_out( DEBUG3 ) <<  " input collection not in event : " << colName << "   - nothing to do for ZPLAN hits  !!! " << std::endl ;  
  } 
  
  if( zplanCol ) { 
    
    
 
   int nHit = zplanCol->getNumberOfElements()  ;
    
    for(int i=0; i< nHit ; i++){
      
      // TrackerHitImpl* hit = dynamic_cast<TrackerHitImpl*>( zplanCol->getElementAt( i ) ) ;
      // if( hit ) {
      // 	streamlog_out( DEBUG ) << "  ---- processing hit : " << *hit << std::endl  
      // 			       << " ------------ "  << std::endl ;
	
      // }
      //       gear::Vector3D   pos(  h->getPosition() ) ;

      lcio::LCObject* h =  zplanCol->getElementAt( i ) ;
      gear::Vector3D   pos ;
      if( ! position( h, pos ) ) continue ;
		 
      //-----------------------------------
      gear::SensorID sensorID ;
      
      bool hitOnLadder =  gearZPLAN.isPointInSensitive( pos ,  &sensorID ) ;

      if( !hitOnLadder ) {
	
	streamlog_out( DEBUG3 ) <<  " hit not in sensitive volume  - distance: " << gearZPLAN.distanceToNearestSensitive( pos ).r() 
				<<  " pos : " << pos << std::endl ;
      }
      //-------------------------------------
      
      


      encoder.setValue( cellID0( h ) ) ;
      streamlog_out( DEBUG2 ) <<  " old cellID: " << encoder.valueString() <<  " - " << cellID0( h )  << std::endl ;
      
      //------------------------------------
      encoder.reset() ;  // reset to 0
      
      encoder[ILDCellID0::subdet] = detID  ;
      encoder[ILDCellID0::side]   =  ( setSideZPlanar ? sensorID.side  : 0 )   ;  // old files don't set the side for the VXD hits properly ... 
      encoder[ILDCellID0::layer]  = sensorID.layer  + layerOffset  ; 
      encoder[ILDCellID0::module] = sensorID.module ;
      encoder[ILDCellID0::sensor] = sensorID.sensor ;
      
      //int layerID = encoder.lowWord() ;  
      // h->setCellID0(  layerID ) ;

      setCellID( h, encoder ) ;
      
      streamlog_out( DEBUG2 ) <<  " new cellID: " << encoder.valueString() << std::endl ;
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

    streamlog_out( DEBUG3 ) <<  " process input collection : " << _tpcColName  << std::endl ;  

  } catch( lcio::DataNotAvailableException& e) { 
    
    streamlog_out( DEBUG3 ) <<  " input collection not in event : " << _tpcColName << "   - nothing to do for TPC hits  !!! " << std::endl ;  
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

    fixZPlanarHits( evt, _vxdColName , ILDDetID::VXD , gearVXD , _layerOffset, _setSideForVXD ) ;
   
  }catch(gear::UnknownParameterException& ){
    streamlog_out( WARNING ) << " FixCellIDs : no gear parameters found for VXD " << std::endl ;
  }

  try{
    const gear::ZPlanarParameters& gearSIT = Global::GEAR->getSITParameters() ;

    fixZPlanarHits( evt, _sitColName , ILDDetID::SIT , gearSIT, _layerOffset ) ;
   
  }catch(gear::UnknownParameterException& ){
    streamlog_out( WARNING ) << " FixCellIDs : no gear parameters found for SIT " << std::endl ;
  }

  try{
    const gear::ZPlanarParameters& gearSET = Global::GEAR->getSETParameters() ;

    fixZPlanarHits( evt, _setColName , ILDDetID::SET , gearSET, _layerOffset  ) ;
   
  }catch(gear::UnknownParameterException& ){
    streamlog_out( WARNING ) << " FixCellIDs : no gear parameters found for SET " << std::endl ;
  }


  

  //=====================================================================================================
  //      FTD
  //=====================================================================================================


  LCCollection* ftdCol = 0 ;

  try{   ftdCol =  dynamic_cast<LCCollection*> (evt->getCollection( _ftdColName )  ); 

  } catch( lcio::DataNotAvailableException& e) { 
    
    streamlog_out( DEBUG5 ) <<  " input collection not in event : " << _ftdColName << "   - nothing to do for FTD hits  !!! " << std::endl ;  
  } 
  
  if( ftdCol ) {

    const gear::FTDParametersImpl& gearFTD = dynamic_cast<const gear::FTDParametersImpl&> ( Global::GEAR->getFTDParameters() )  ;
    const gear::FTDLayerLayout& layerLayout = gearFTD.getFTDLayerLayout() ;
    
    int nHit = ftdCol->getNumberOfElements()  ;
    
    for(int i=0; i< nHit ; i++){
      
      lcio::LCObject* h =  ftdCol->getElementAt( i ) ;
      gear::Vector3D   pos ;
      if( ! position( h, pos ) ) continue ;

      // TrackerHitImpl* h = dynamic_cast<TrackerHitImpl*>( ftdCol->getElementAt( i ) ) ;
      // gear::Vector3D   pos(  h->getPosition() ) ;

               
      encoder.setValue( cellID0( h ) ) ;
      streamlog_out( DEBUG5 ) <<  " old cellID: " << encoder.valueString() <<  " - " << cellID0( h )  << std::endl ;
      
      // int module = gearFTD.getPetalIndex( pos ) ;
       
      // int layer = gearFTD.getLayerIndex( pos ) ;

      // encoder.reset() ;  // reset to 0
      
      // encoder[ILDCellID0::subdet] = ILDDetID::FTD ;
      // encoder[ILDCellID0::side]   = pos.z() > 0 ? 1 : -1  ;
      // encoder[ILDCellID0::layer]  = layer ;
      // encoder[ILDCellID0::module] =  module ;
      // encoder[ILDCellID0::sensor] = 0 ; //fixme ???
      
      //      int layerID = encoder.lowWord() ;  
      //      h->setCellID0(  layerID ) ;

      //      setCellID( h, encoder ) ;
      
      ftdCol->parameters().setValue( "CellIDEncoding" , ILDCellID0::encoder_string ) ;
      
    } 
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


