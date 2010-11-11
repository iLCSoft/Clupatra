#include "EXTPCKalDetector.h"
#include "EXTPCMeasLayer.h"
#include "EXTPCHit.h"

#include "TRandom.h"
#include "TMath.h"
#include "TTUBE.h"
#include "TNode.h"
#include "TVirtualPad.h"

#include "KalTest.h"

#include <sstream>
#include <iomanip>


#include "gear/GEAR.h"
#include "gear/TPCParameters.h"
#include "gear/PadRowLayout2D.h"
#include "gearimpl/Util.h"

#include "streamlog/streamlog.h"

Double_t EXTPCKalDetector::fgVdrift = 5.e-3;
//ClassImp(EXTPCKalDetector)

EXTPCKalDetector::EXTPCKalDetector(Int_t m)
                : EXVKalDetector(m),
		  fNodePtr(0)
{
//    Double_t A, Z, density, radlen;
//    A       = 14.00674 * 0.7 + 15.9994 * 0.3;
//    Z       = 7.3;
//    density = 1.205e-3;
//    radlen  = 3.42e4;
//    TMaterial &air = *new TMaterial("TPCAir", "", A, Z, density, radlen, 0.);

//    A       = 39.948*0.9+(12.011*0.2+1.00794*0.8)*0.1;
//    Z       = 16.4;
//    density = 0.749e-3;
//    radlen  =  1.196e4*2;
//    TMaterial &gas = *new TMaterial("TPCGas", "", A, Z, density, radlen, 0.);

//    A       = 12.0107;
//    Z       =  6.;
//    density = 0.1317;
//    radlen  = 42.7/density;
//    TMaterial &cfrp = *new TMaterial("TPCCFRP", "", A, Z, density, radlen, 0.);

//    static const Double_t inthick   = 2.1075 * 2.;   // thick of inner shell
//    static const Double_t outthick  = 4.1175 * 2.;   // thick of outer shell
//    static const Int_t    nlayers   = 200;           // number of layer
//    static const Double_t lhalf     = 255.;          // half length
//    static const Double_t rmin      = 44.215;        // minimum radius
//    static const Double_t rstep     = 0.76775;       // step length of radius
//    static const Double_t rtub      = 39.5;          // inner r of support tube
//    static const Double_t outerr    = 206.;          // outer radius of TPC
// //   static const Double_t sigmax0   = 55.e-4;
// //   static const Double_t sigmax1   = 166.e-4 / 3 / TMath::Sqrt(28);
// //   static const Double_t sigmaz    = 600.e-4;
//    static const Double_t sigmax0   = .02 ; // 55.e-4;
//    static const Double_t sigmax1   = .02 ; // 166.e-4 / 3 / TMath::Sqrt(28);
//    static const Double_t sigmaz    = .06 ; // 600.e-4;

//    Bool_t active = EXTPCMeasLayer::kActive;
//    Bool_t dummy  = EXTPCMeasLayer::kDummy;
//    Add(new EXTPCMeasLayer(air, cfrp, rtub, lhalf, sigmax0, sigmax1, sigmaz, dummy));
//    Add(new EXTPCMeasLayer(cfrp, gas, rtub+inthick, lhalf, sigmax0, sigmax1, sigmaz, dummy));

//    // create measurement layers of central tracker
//    Double_t r = rmin;
//    for (Int_t layer = 0; layer < nlayers; layer++) {
//       std::stringstream ss;
//       ss << "TPC" << std::setw(3) << std::setfill('0') << layer << std::ends;
//       Add(new EXTPCMeasLayer(gas, gas, r, lhalf, sigmax0, sigmax1, sigmaz, active, ss.str().data()));
//       r += rstep;
//    }
//    Add(new EXTPCMeasLayer(gas, cfrp, outerr-outthick, lhalf, sigmax0, sigmax1, sigmaz, dummy));
//    Add(new EXTPCMeasLayer(cfrp, air, outerr, lhalf, sigmax0, sigmax1, sigmaz, dummy));

//    SetOwner();
}


EXTPCKalDetector::EXTPCKalDetector(const gear::TPCParameters& tpcParams ) : 
  EXVKalDetector(250), // initial size 
  fNodePtr(0) {
  

  //fg: get as many parameters from GEAR as possible ...
  const gear::PadRowLayout2D& pL = tpcParams.getPadLayout() ; 

  Double_t A, Z, density, radlen;
  A       = 14.00674 * 0.7 + 15.9994 * 0.3;
  Z       = 7.3;
  density = 1.205e-3;
  radlen  = 3.42e4;
  TMaterial &air = *new TMaterial("TPCAir", "", A, Z, density, radlen, 0.);
  
  A       = 39.948*0.9+(12.011*0.2+1.00794*0.8)*0.1;
  Z       = 16.4;
  density = 0.749e-3 ;
  radlen  =  1.196e4*2;
  TMaterial &gas = *new TMaterial("TPCGas", "", A, Z, density, radlen, 0.);
  
  A       = 12.0107;
  Z       =  6.;
  density = 0.1317;
  radlen  =  tpcParams.getDoubleVal("TPCWallProperties_RadLen")  ; // 42.7/density;
  TMaterial &cfrp = *new TMaterial("TPCCFRP", "", A, Z, density, radlen, 0.);
  

  // FIXME: need to adjust the dEDx through A and density ....

  
  static const Int_t    nlayers   =  pL.getNRows() ;   // n rows
  static const Double_t lhalf     =  tpcParams.getMaxDriftLength() ;    // 255. // half length

  static const Double_t rstep     =  pL.getRowHeight(0) ;   // 0.76775;       // step length of radius

  // assuming that this is the radius of the first measurment layer ....
  static const Double_t rmin      =  tpcParams.getPlaneExtent()[0]   + rstep/2. ;  // 44.215;        // minimum radius

  streamlog_out( DEBUG4 ) << tpcParams << std::endl ;
  
  static const Double_t rtub      = tpcParams.getDoubleVal("tpcInnerRadius")  ; // 39.5; // inner r of support tube
  static const Double_t outerr    = tpcParams.getDoubleVal("tpcOuterRadius")  ; //206.; // outer radius of TPC

  static const Double_t inthick   =  tpcParams.getDoubleVal("tpcInnerWallThickness")  ; //2.1075 * 2.;   // thick of inner shell
  static const Double_t outthick  =  tpcParams.getDoubleVal("tpcOuterWallThickness")  ; //4.1175 * 2.;   // thick of outer shell
  

  static const Double_t sigmax0   = .02 ; // 55.e-4;
  static const Double_t sigmax1   = .02 ; // 166.e-4 / 3 / TMath::Sqrt(28);
  static const Double_t sigmaz0    = .06 ; // 600.e-4;
  static const Double_t sigmaz1    = .06 ; // 600.e-4;      
  
  //   static const Double_t sigmax0   = 55.e-4;
  //   static const Double_t sigmax1   = 166.e-4 / 3 / TMath::Sqrt(28);
  //   static const Double_t sigmaz    = 600.e-4;
  
  Bool_t active = EXTPCMeasLayer::kActive;
  Bool_t dummy  = EXTPCMeasLayer::kDummy;
  
  //FIXME - test:  add a layer inside the beam pipe 
  Add(new EXTPCMeasLayer(air, air, 1.2 , lhalf, sigmax0, sigmax1, sigmaz0, sigmaz1, active , -1 )) ;  // ,ss.str().data()));


  Add(new EXTPCMeasLayer(air, cfrp, rtub, lhalf, sigmax0, sigmax1, sigmaz0, sigmaz1, dummy, -1 ));
  Add(new EXTPCMeasLayer(cfrp, gas, rtub+inthick, lhalf, sigmax0, sigmax1, sigmaz0, sigmaz1,dummy, -1 ));
  
  // create measurement layers
  Double_t r = rmin;

  static const double gasdEdx  =  tpcParams.getDoubleVal("TPCGasProperties_dEdx")  ; 
  
  streamlog_out( DEBUG ) << " using dEdx for TPC gas : " << gasdEdx << " Gev/mm " << std::endl ;


  for (Int_t layer = 0; layer < nlayers; layer++) {
    
    int layerID = KalTest::DetID::TPC * KalTest::DetID::Factor  + layer ;
    
    EXTPCMeasLayer* tpcL =  new EXTPCMeasLayer(gas, gas, r, lhalf, sigmax0, sigmax1, sigmaz0, sigmaz1, active , layerID ) ;

    tpcL->setdEdx_GeV_mm(  gasdEdx  ) ;

    Add( tpcL ) ;  
    
    //    if( streamlog_level( DEBUG0 ) &&
    if( layer % 10 == 0 ){
      
      streamlog_out( DEBUG0)   << " ***** adding TPC layer : [" << layer + layerID  <<  "] at R = " << r << std::endl ;
    }
    
    r += rstep;
  }
  Add(new EXTPCMeasLayer(gas, cfrp, outerr-outthick, lhalf, sigmax0, sigmax1, sigmaz0, sigmaz1,  dummy, -1 ));
  Add(new EXTPCMeasLayer(cfrp, air, outerr, lhalf, sigmax0, sigmax1, sigmaz0, sigmaz1, dummy, -1 ));
  
  SetOwner();
}





EXTPCKalDetector::~EXTPCKalDetector()
{
}

// -----------------
//  Draw
// -----------------
//    Drawing method for event display
//
void EXTPCKalDetector::Draw(Int_t color, const Char_t *opt)
{
   if (!gPad) return;
   TNode *nodep = GetNodePtr();
   nodep->cd();

   if (!fNodePtr) {
      EXTPCMeasLayer *inp  = static_cast<EXTPCMeasLayer *>(First());
      EXTPCMeasLayer *outp = static_cast<EXTPCMeasLayer *>(Last());
      Double_t rin  = inp->GetR();
      Double_t rout = outp->GetR();
      Double_t hlen = outp->GetZmax();
      const Char_t *name  = "TPC";
      const Char_t *nname = "TPCNode";
      TTUBE *tubep = new TTUBE(name,name,"void",rin,rout,hlen);
      tubep->SetBit(kCanDelete);
      fNodePtr = new TNode(nname,nname,name);
      fNodePtr->SetLineColor(color);
      fNodePtr->SetLineWidth((Width_t)0.01);
   }
   EXVKalDetector::Draw(color,opt);
}
