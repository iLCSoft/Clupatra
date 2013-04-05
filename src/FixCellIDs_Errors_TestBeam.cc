#include "FixCellIDs_Errors_TestBeam.h"
#include <iostream>
#include <cmath>

#include <EVENT/LCCollection.h>
#include <IMPL/TrackerHitImpl.h>
#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>
#include <UTIL/LCTOOLS.h>

//---- GEAR ----
#include "marlin/Global.h"
#include "gear/GEAR.h"
#include "gear/TPCParameters.h"
#include "gear/TPCModule.h"
#include "gear/BField.h"
#include "gearimpl/Vector3D.h"

using namespace lcio ;
using namespace marlin ;


FixCellIDs_Errors_TestBeam aFixCellIDs_Errors_TestBeam ;


FixCellIDs_Errors_TestBeam::FixCellIDs_Errors_TestBeam() : Processor("FixCellIDs_Errors_TestBeam") ,
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
  
}


void FixCellIDs_Errors_TestBeam::init() {
  
  streamlog_out(DEBUG) << "   init called  " << std::endl ;
  
  // usually a good idea to
  printParameters() ;
  
  _nRun = 0 ;
  _nEvt = 0 ;
  
}


void FixCellIDs_Errors_TestBeam::processRunHeader( LCRunHeader* run) {
  
  _nRun++ ;
} 



void FixCellIDs_Errors_TestBeam::modifyEvent( LCEvent * evt ) {


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

    UTIL::BitField64 encoder( lcio::ILDCellID0::encoder_string ) ;

    int nHit = tpcCol->getNumberOfElements()  ;

    for(int i=0; i< nHit ; i++){

      TrackerHitImpl* h = dynamic_cast<TrackerHitImpl*>( tpcCol->getElementAt( i ) ) ;

      const double* pos = h->getPosition();
      // Note: here we take the Cartesian coordinates as specified in the gear file!
      const gear::TPCModule& tpcmodule = gearTPC.getNearestModule(pos[0], pos[1]);


      int mod = tpcmodule.getModuleID();
      int pad = tpcmodule.getNearestPad(pos[0], pos[1]);
      int row = tpcmodule.getRowNumber(pad);

      int row_global =
  	// Modules 0 and 1 get the same row index
  	( mod==0 || mod==1 ) ? row :
		/* For modules 2, 3 and 4 the row get increased
		 * by the number of rows in the modules 0/1
		 */
        ( ( mod==2 || mod==3 || mod==4) ? gearTPC.getModule(0).getNRows() + row :
		/* For modules 5 and 6 the row get increased
		 * by the number of rows in the modules 0/1 and 2/3/4
		 */
  	gearTPC.getModule(0).getNRows() + gearTPC.getModule(2).getNRows() + row );


      encoder.reset() ;  // reset to 0
      encoder[ILDCellID0::subdet] = ILDDetID::TPC ;
      encoder[ILDCellID0::layer] = row_global ;

      int CellID = encoder.lowWord() ;
      h->setCellID0(CellID);

      // Determine the covariance matrix according to the given resolution
      double zdrift = gearTPC.getMaxDriftLength() - (std::abs(pos[2]));

      // the following code for getting the sigma parameters for the uncertainties
      // were copied from GearTPCKalDetector
      //
      // The resolution parameters can vary from module to module.
      // In case they are not there use the jgem values
      double sigmax0, sigmax1, sigmaz0, sigmaz1;
      try{
        sigmax0 = tpcmodule.getDoubleVal("sigmax0");
        sigmax1 = tpcmodule.getDoubleVal("sigmax1");
        sigmaz0 = tpcmodule.getDoubleVal("sigmaz0");
        sigmaz1 = tpcmodule.getDoubleVal("sigmaz1");
      }
      catch( gear::UnknownParameterException & )
      {
        streamlog_out(MESSAGE) << "No resolution parameters found for module "
                   << tpcmodule.getModuleID()<<"."
                   << " Using jgem settings." << std::endl;

        // FIXME: n_eff of argon per distance, adapt to pad length. Or read from GEAR file?
        double neff    = 22.7;

//        // FIXME : I think the sqrt(10) is a leftover from the old cm units....
//        sigmax0 = 38.3e-3;
//        sigmax1 = 101.5e-3 / sqrt(10.) / sqrt(neff);
//        sigmaz0 = 500.e-3 ;
//        sigmaz1 = 154.e-3  / sqrt(10.) / sqrt(neff);

        // Alex: here I actually took the values from Isa's calculations, not from GearTPCKalDetector
        sigmax0 = 0.02;
        sigmax1 = 0.0198;
        sigmaz0 = 0.4; // rough value, not precisely extracted from the plot
        sigmaz1 = 0.0; // assume sigma_z constant
      }

      // Sigma squared along Z
      double zSigmaSquare = sigmaz0 * sigmaz0 + sigmaz1 * sigmaz1 * zdrift;
      // Sigma squared in RPhi plane
      double rphiSigmaSquare = sigmax0 * sigmax0 + sigmax1 * sigmax1 * zdrift;

      // The phi angle of the hit.
      // FIXME: probably the phi should be determined in local coordinate system rather than in the global.
      // It makes differences if the measurement layers are quite displaced from each other.
      // This should however not matter for the total hit uncertainty in the RPhi plane as
      // in the relevant places only the sum of covMatrix[0] and covMatrix[2] is taken,
      // thus the resulting sin^2+cos^2 cancels to 1
      gear::Vector3D v(pos);
      double phi = v.phi();
      
      float covMatrix[6] = {0., 0., 0., 0., 0., 0.};

      covMatrix[0] = pow( cos(phi), 2) * rphiSigmaSquare;
      covMatrix[1] = sin(phi) * cos(phi) * rphiSigmaSquare;
      covMatrix[2] = pow( sin(phi), 2) * rphiSigmaSquare;
      covMatrix[5] = float(zSigmaSquare);

      h->setCovMatrix(covMatrix);

      tpcCol->parameters().setValue( "CellIDEncoding" , ILDCellID0::encoder_string ) ;
    } 
  }

  streamlog_out( DEBUG3 ) << "   processing event: " << evt->getEventNumber() 
			  << "   in run:  " << evt->getRunNumber() << std::endl ;

  _nEvt ++ ;
}



void FixCellIDs_Errors_TestBeam::check( LCEvent * evt ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void FixCellIDs_Errors_TestBeam::end(){

  streamlog_out( MESSAGE ) << "FixCellIDs_Errors_TestBeam::end()  " << name()
			   << " processed " << _nEvt << " events in " << _nRun << " runs "
			   << std::endl ;
  
}


