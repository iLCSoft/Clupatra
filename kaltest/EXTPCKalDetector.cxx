
#include "EXTPCKalDetector.h"
#include "EXTPCMeasLayer.h"
#include "EXTPCHit.h"
#include "TRandom.h"
#include "TMath.h"

#include "TTUBE.h"
#include "TNode.h"
#include "TVirtualPad.h"

#include <sstream>
#include <iomanip>


#include "gear/GEAR.h"
#include "gear/TPCParameters.h"
#include "gear/PadRowLayout2D.h"

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

  static const double cm = 0.1 ; // conversion from mm to cm
  //  fgVdrift = tpcParam.getDriftVelocity() * cm //fg: this is not set in current gear files ! ???  



  Double_t A, Z, density, radlen;
  A       = 14.00674 * 0.7 + 15.9994 * 0.3;
  Z       = 7.3;
  density = 1.205e-3;
  radlen  = 3.42e4;
  TMaterial &air = *new TMaterial("TPCAir", "", A, Z, density, radlen, 0.);
  
  A       = 39.948*0.9+(12.011*0.2+1.00794*0.8)*0.1;
  Z       = 16.4;
  density = 0.749e-3;
  radlen  =  1.196e4*2;
  TMaterial &gas = *new TMaterial("TPCGas", "", A, Z, density, radlen, 0.);
  
  A       = 12.0107;
  Z       =  6.;
  density = 0.1317;
  radlen  =  tpcParams.getDoubleVal("TPCWallProperties_RadLen")  ; // 42.7/density;
  TMaterial &cfrp = *new TMaterial("TPCCFRP", "", A, Z, density, radlen, 0.);
  

  // FIXME: need to adjust the dEDx through A and density ....

  
  static const Int_t    nlayers   =  pL.getNRows() ;   // n rows
  static const Double_t lhalf     =  tpcParams.getMaxDriftLength() * cm  ;    // 255. // half length

  static const Double_t rstep     =  pL.getRowHeight(0) * cm    ;   // 0.76775;       // step length of radius

  // assuming that this is the radius of the first measurment layer ....
  static const Double_t rmin      =  tpcParams.getPlaneExtent()[0] * cm  - rstep/2. ;  // 44.215;        // minimum radius

  streamlog_out( DEBUG4 )  << " ***** TPC rmin  = " << rmin 
			   << " tpcParams.getPlaneExtent()[0] " << tpcParams.getPlaneExtent()[0] 
			   << "  pL.getRowHeight(0) "  <<  pL.getRowHeight(0)
			   << std::endl ;
  
  static const Double_t rtub      = tpcParams.getDoubleVal("tpcInnerRadius") * cm  ; // 39.5; // inner r of support tube
  static const Double_t outerr    = tpcParams.getDoubleVal("tpcOuterRadius") * cm  ; //206.; // outer radius of TPC

  static const Double_t inthick   =  tpcParams.getDoubleVal("tpcInnerWallThickness") * cm  ; //2.1075 * 2.;   // thick of inner shell
  static const Double_t outthick  =  tpcParams.getDoubleVal("tpcOuterWallThickness") * cm  ; //4.1175 * 2.;   // thick of outer shell
  
  //    ... from ILD_00 .... 
//             <parameter name="TPCGasProperties_RadLen" type="double" value="1.155205835e+05" />
//             <parameter name="TPCGasProperties_dEdx" type="double" value="2.669325944e-07" />
//             <parameter name="TPCWallProperties_RadLen" type="double" value="8.896320560e+01" />
//             <parameter name="TPCWallProperties_dEdx" type="double" value="4.329175517e-04" />
//             <parameter name="tpcInnerRadius" type="double" value="3.290000000e+02" />
//             <parameter name="tpcInnerWallThickness" type="double" value="1.160000000e+00" />
//             <parameter name="tpcIonPotential" type="double" value="3.200000000e-08" />
//             <parameter name="tpcOuterRadius" type="double" value="1.808000000e+03" />
//             <parameter name="tpcOuterWallThickness" type="double" value="1.780000000e+00" />


   static const Double_t sigmax0   = .02 ; // 55.e-4;
   static const Double_t sigmax1   = .02 ; // 166.e-4 / 3 / TMath::Sqrt(28);
   static const Double_t sigmaz0    = .06 ; // 600.e-4;
   static const Double_t sigmaz1    = .06 ; // 600.e-4;      

//   static const Double_t sigmax0   = 55.e-4;
//   static const Double_t sigmax1   = 166.e-4 / 3 / TMath::Sqrt(28);
//   static const Double_t sigmaz    = 600.e-4;
  
  Bool_t active = EXTPCMeasLayer::kActive;
  Bool_t dummy  = EXTPCMeasLayer::kDummy;

  //FIXME - test:  add a layer for the beam pipe 
  Add(new EXTPCMeasLayer(air, air, 1.2 , lhalf, sigmax0, sigmax1, sigmaz0, sigmaz1, active )) ;  // ,ss.str().data()));


  Add(new EXTPCMeasLayer(air, cfrp, rtub, lhalf, sigmax0, sigmax1, sigmaz0, sigmaz1, dummy));
  Add(new EXTPCMeasLayer(cfrp, gas, rtub+inthick, lhalf, sigmax0, sigmax1, sigmaz0, sigmaz1,dummy));
  
  // create measurement layers
  Double_t r = rmin;
  for (Int_t layer = 0; layer < nlayers; layer++) {

    //    std::stringstream ss;
      
    Add(new EXTPCMeasLayer(gas, gas, r, lhalf, sigmax0, sigmax1, sigmaz0, sigmaz1, active )) ;  // ,ss.str().data()));
    r += rstep;

    // safe offset of first measurment layer in TPC :
    if( layer == 0 ) 
      _layerOffset = GetLast()  ;
    
    if( streamlog_level( DEBUG0 ) ) { // && layer % 10 == 0 ){
      streamlog_out( DEBUG0)   << " ***** adding TPC layer : [" << layer +_layerOffset <<  "] at R = " << r << std::endl ;
    }
  }
  Add(new EXTPCMeasLayer(gas, cfrp, outerr-outthick, lhalf, sigmax0, sigmax1, sigmaz0, sigmaz1,  dummy));
  Add(new EXTPCMeasLayer(cfrp, air, outerr, lhalf, sigmax0, sigmax1, sigmaz0, sigmaz1, dummy));
  
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
      fNodePtr->SetLineWidth(0.01);
   }
   EXVKalDetector::Draw(color,opt);
}
