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

#include "ConfigFlags.h"



KalTest::KalTest( const gear::GearMgr& gearMgr) :  
  _gearMgr( &gearMgr ) , 
  _idOffsets(  DetID::Size , -1 ) {
  
  streamlog_out( DEBUG4 ) << "  KalTest - initializing the detector ..." << std::endl ;
  
  _det = new TKalDetCradle ;
  _det->SetOwner( true ) ; // takes care of deleting subdetector in the end ...
  

  _cfg = new ConfigFlags(  CFG::size ) ;
  
  _cfg->registerOption( KalTest::CFG::ownsHits    , "KalTest::CFG::ownsHits"  , true ) ;
  _cfg->registerOption( KalTest::CFG::useQMS      , "KalTest::CFG::useQMS"    , false ) ;
  _cfg->registerOption( KalTest::CFG::usedEdx     , "KalTest::CFG::usedEdx"   , false ) ;
  
  streamlog_out( DEBUG4 ) << " ******  KalTest configuration : " << std::endl 
			  <<   *_cfg
 			  <<   std::endl ;


  // this  now a public init() method taking options ....
  // init() 
}


void KalTest::setOption(unsigned CFGOption, bool val) { 

  _cfg->setOption( CFGOption , val ) ;
}   

KalTest::~KalTest(){
  delete _cfg ;
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
  
  // --- set some options :

  if( ! (*_cfg)[ KalTest::CFG::useQMS ]  )   _det->SwitchOffMS();    

  if( ! (*_cfg)[ KalTest::CFG::usedEdx ]  )  _det->SwitchOffDEDX();  
  
}

KalTrack* KalTest::createKalTrack()  { 

  return new KalTrack( _det , _cfg->option( KalTest::CFG::ownsHits )  ) ; 
}



TVTrackHit*  KalTest::createHit(EVENT::TrackerHit* hit, int layer, int detID) {
  
  layer +=  indexOfFirstLayer( detID ) ;
  
  if( layer >= 0 && ( _det->GetEntries() > layer ) ) {
    
    TObject* o =  _det->At( layer ) ;
    
    EXVMeasLayer* ml = dynamic_cast< EXVMeasLayer * >( o ) ;
    
    if (ml != 0 && ml->IsActive() && dynamic_cast<const EXVKalDetector &>( ml->GetParent(kFALSE) ).IsPowerOn() ) {
      
      return ml->createHit(  hit ) ;
    }
  }
  
  return 0 ; 
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
