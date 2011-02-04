#include "KalTest.h"

 #include "kaltest/TKalDetCradle.h"

#include "EXVMeasLayer.h"
#include "EXVTXKalDetector.h"
#include "EXTPCKalDetector.h"

#include "gear/GEAR.h"
#include "gear/BField.h"
#include "gear/TPCParameters.h"
#include "gear/PadRowLayout2D.h"

#include <math.h>
#include <cmath>

#include "streamlog/streamlog.h"



KalTest::KalTest( const gear::GearMgr& gearMgr) :  
  _gearMgr( &gearMgr ) , 
  _idOffsets(  DetID::Size , -1 ) {
  
  streamlog_out( DEBUG4 ) << "  KalTest - initializing the detector ..." << std::endl ;
  
  _det = new TKalDetCradle ;
  _det->SetOwner( true ) ; // takes care of deleting subdetector in the end ...
  

  // this could be made a public init() method taking options ....
  init() ;
  
}

KalTest::~KalTest(){
  
  delete _det ;
}


void KalTest::init() {
  
  const gear::TPCParameters& tpcParams = _gearMgr->getTPCParameters();
  
  EXVTXKalDetector* vtxdet = new EXVTXKalDetector ;
  EXTPCKalDetector* tpcdet = new EXTPCKalDetector( tpcParams )  ;
  
  tpcdet->SetBField(   _gearMgr->getBField().at( gear::Vector3D( 0.,0.,0.)  ).z()  ) ; //  Tesla 
 
  // FIXME: DEBUG.....
  _det->Install( *vtxdet ) ;  
  _det->Install( *tpcdet ) ;  
  
  _det->Close() ;          // close the cradle
  _det->Sort() ;           // sort meas. layers from inside to outside
  

  streamlog_out( DEBUG3 ) << " created TPC detector; layerOffset =  " 
			  << indexOfFirstLayer( DetID::TPC ) 
			  << std::endl ; 



  //vtxdet->PowerOff();       // power off vtx not to process hit
  
  // --- possible options :
  // _det->SwitchOffMS();    // switch off multiple scattering
  // _det->SwitchOffDEDX();  // switch off enery loss
  
}

int KalTest::indexOfFirstLayer(int detID) { 
  
  if( _idOffsets.at( detID ) == -1 ){
    
    int nDet = _det->GetEntriesFast() ;
    
    int firstID = detID * DetID::Factor ;
    
    for( int i=0  ; i < nDet ; ++i ){
      
      if( dynamic_cast<EXVMeasLayer*>( _det->At( i ) )->getLayerID() == firstID  ){
	
	_idOffsets[ detID ] = i ;
	break ;
      }
    }

  }
  return _idOffsets[ detID ] ; 
 }


int  KalTest::maxLayerIndex() {

    return _det->GetEntriesFast() ;

}
