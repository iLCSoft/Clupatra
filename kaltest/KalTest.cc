#include "KalTest.h"

#include "TKalTrackState.h"
#include "TKalTrackSite.h"
#include "TVTrackHit.h"
#include "TKalDetCradle.h"

#include "EXVMeasLayer.h"
#include "EXVTXKalDetector.h"
#include "EXTPCKalDetector.h"
#include "EXVTXHit.h"
#include "EXTPCHit.h"
#include "EXHYBTrack.h"

#include "gear/GEAR.h"
#include "gear/BField.h"
#include "gear/TPCParameters.h"
#include "gear/PadRowLayout2D.h"

#include "IMPL/TrackImpl.h"
//#include "lcio.h"

#include <math.h>
#include <cmath>

#include "streamlog/streamlog.h"

//using namespace lcio ;

KalTest::KalTest( const gear::GearMgr& gearMgr) :  _gearMgr( &gearMgr ) {
  
  streamlog_out( DEBUG4 ) << "  KalTest - initializing the detector ..." << std::endl ;
  
  _det = new TKalDetCradle ;
  _det->SetOwner( true ) ; // takes care of deleting subdetector in the end ...
  
  _kalHits = new TObjArray ;
  
  // this could be made a public init() method taking options ....
  init() ;
  
}

KalTest::~KalTest(){
  
  delete _det ;
  delete _kalHits ;
}


void KalTest::init() {
  
  const gear::TPCParameters& tpcParams = _gearMgr->getTPCParameters();
  
  EXVTXKalDetector* vtxdet = new EXVTXKalDetector ;
  EXTPCKalDetector* tpcdet = new EXTPCKalDetector( tpcParams )  ;
  
  tpcdet->SetBField(   _gearMgr->getBField().at( gear::Vector3D( 0.,0.,0.)  ).z() * 10 ) ; // 10 for Tesla -> kGauss

  //_det->Install( *vtxdet ) ;  
  _det->Install( *tpcdet ) ;  
  
  _det->Close() ;          // close the cradle
  _det->Sort() ;           // sort meas. layers from inside to outside
  
  vtxdet->PowerOff();       // power off vtx not to process hit
  
  // --- possible options :
  //_det->SwitchOffMS();    // switch off multiple scattering
  //_det->SwitchOffDEDX();  // switch off enery loss
  
}



void KalTest::addHit( const TVector3& pos, int layer ) {
  
  if( layer >= 0 && ( _det->GetEntries() > layer ) ) {

    TObject* o =  _det->At( layer ) ;

    EXVMeasLayer* ml = dynamic_cast< EXVMeasLayer * >( o ) ;

  if (ml != 0 && ml->IsActive() && dynamic_cast<const EXVKalDetector &>( ml->GetParent(kFALSE) ).IsPowerOn() ) {
    
    ml->ProcessHit( pos, *_kalHits ); // create hit point
    


    if( layer % 10 == 0 ){
      double radius = pos.Perp() ;
      std::cout << " ***** adding TPC hit in layer : [" << layer <<  "] at R = " << radius << std::endl ;
    }


  } else {
    
    streamlog_out( WARNING ) << " hit not added to KalTest at : " 
   			     << pos(0) << "," << pos(1) << "," << pos(2)  
 			     << "  in  layer: " << layer
 			     << " ml : " << ml 
 			     << "  ml->IsActive() " <<  ml->IsActive()
 			     << "dynamic_cast<const EXVKalDetector &>( ml->GetParent(kFALSE) ).IsPowerOn() " 
      //<<   & ( ml->GetParent(kFALSE) )
			     << dynamic_cast<const EXVKalDetector &>( ml->GetParent(kFALSE) ).IsPowerOn()
 			     << std::endl ;
  } 
  
} else {
  
  streamlog_out( ERROR ) << " no measurment layer at : " 
			 << pos(0) << "," << pos(1) << "," << pos(2)  
			 << "  in  layer: " << layer 
			 << " detector has  " << _det->GetEntriesFast() << " layers only "
 			 << " - object at given layer : " <<  _det->At( layer )  << std::endl ;
 }
}


void KalTest::fitTrack(IMPL::TrackImpl* trk) {
  
  const Bool_t gkDir = kIterBackward;
  
  using namespace std ;
  
  // ============================================================
  //  Do Kalman Filter    - copied from EXVKalTest.cxx
  // ============================================================

  if (_kalHits->GetEntries() < 3) {

    streamlog_out( ERROR)  << "<<<<<< KalTest::fitTrack(): Shortage of Hits! nhits = " 
			   << _kalHits->GetEntries() << " >>>>>>>" << endl;
    return ;
  }
      
  Int_t i1, i2, i3; // (i1,i2,i3) = (1st,mid,last) hit to filter
  if (gkDir == kIterBackward) {
    i3 = 0;
    i1 = _kalHits->GetEntries() - 1;
    i2 = i1 / 2;
  } else {
    i1 = 0;
    i3 = _kalHits->GetEntries() - 1;
    i2 = i3 / 2;
  }
      
  // ---------------------------
  //  Create a dummy site: sited
  // ---------------------------
      
  TVTrackHit *ht1p = dynamic_cast<TVTrackHit *>(_kalHits->At(i1));
  TVTrackHit *htdp = 0;
  if (dynamic_cast<EXVTXHit *>(ht1p)) {
    htdp = new EXVTXHit(*dynamic_cast<EXVTXHit *>(ht1p));
    //       } else if (dynamic_cast<EXITHit *>(ht1p)) {
    //          htdp = new EXITHit(*dynamic_cast<EXITHit *>(ht1p));
    //       } else if (dynamic_cast<EXITFBHit *>(ht1p)) {
    //          htdp = new EXITFBHit(*dynamic_cast<EXITFBHit *>(ht1p));
  } else if (dynamic_cast<EXTPCHit *>(ht1p)) {
    htdp = new EXTPCHit(*dynamic_cast<EXTPCHit *>(ht1p));
  }
  TVTrackHit &hitd = *htdp;

  hitd(0,1) = 1.e6;   // give a huge error to d
  hitd(1,1) = 1.e6;   // give a huge error to z

  TKalTrackSite &sited = *new TKalTrackSite(hitd);
  sited.SetHitOwner();// site owns hit
  sited.SetOwner();   // site owns states

  // ---------------------------
  //  Create initial helix
  // ---------------------------

  TVTrackHit &h1 = *dynamic_cast<TVTrackHit *>(_kalHits->At(i1)); // first hit
  TVTrackHit &h2 = *dynamic_cast<TVTrackHit *>(_kalHits->At(i2)); // middle hit
  TVTrackHit &h3 = *dynamic_cast<TVTrackHit *>(_kalHits->At(i3)); // last hit
  TVector3    x1 = h1.GetMeasLayer().HitToXv(h1);
  TVector3    x2 = h2.GetMeasLayer().HitToXv(h2);
  TVector3    x3 = h3.GetMeasLayer().HitToXv(h3);

  double Bfield =   h1.GetBfield() ;

  THelicalTrack helstart(x1, x2, x3, h1.GetBfield(), gkDir); // initial helix 

  // ---------------------------
  //  Set dummy state to sited
  // ---------------------------

  static TKalMatrix svd(kSdim,1);
  svd(0,0) = 0.;                        // dr
  svd(1,0) = helstart.GetPhi0();        // phi0
  svd(2,0) = helstart.GetKappa();       // kappa
  svd(3,0) = 0.;                        // dz
  svd(4,0) = helstart.GetTanLambda();   // tan(lambda)
  if (kSdim == 6) svd(5,0) = 0.;        // t0

  static TKalMatrix C(kSdim,kSdim);
  for (Int_t i=0; i<kSdim; i++) {
    C(i,i) = 1.e4;   // dummy error matrix
  }

  sited.Add(new TKalTrackState(svd,C,sited,TVKalSite::kPredicted));
  sited.Add(new TKalTrackState(svd,C,sited,TVKalSite::kFiltered));

  // ---------------------------
  //  Add sited to the kaltrack
  // ---------------------------

  EXHYBTrack kaltrack;   // a track is a kal system
  kaltrack.SetOwner();   // kaltrack owns sites
  kaltrack.Add(&sited);  // add the dummy site to this track

  // ---------------------------
  //  Prepare hit iterrator
  // ---------------------------

  TIter next(_kalHits, gkDir); // come in to IP, if gkDir = kIterBackward

  // ---------------------------
  //  Start Kalman Filter
  // ---------------------------

  TVTrackHit *hitp = 0;
  while ((hitp = dynamic_cast<TVTrackHit *>(next()))) {
    TKalTrackSite  &site = *new TKalTrackSite(*hitp); // new site

    if (!kaltrack.AddAndFilter(site)) {               // filter it
      cout << " site discarded!" << endl;
      delete &site;                        // delete it if failed
    }
  } // end of Kalman filter

  //fg: set origin as reference point ...
//    TKalTrackSite* origin = new TKalTrackSite ;
//    origin->SetPivot( TVector3( 0., 0., 0. ) ) ;
//    //   kaltrack.AddAndFilter( *origin );

  // ---------------------------
  //  Smooth the track
  // ---------------------------
#define SAVE_RESIDUAL false
#ifndef SAVE_RESIDUAL
  TVKalSite &cursite = kaltrack.GetCurSite();
#else
  Int_t isite = 1;
  kaltrack.SmoothBackTo(isite);
  TVKalSite &cursite = static_cast<TVKalSite &>(*kaltrack[isite]);
#endif


   TVKalState& trkState = cursite.GetCurState() ;
//   trkState.Propagate( *origin ) ;


  // dump fit result for now:
  Int_t    ndf  = kaltrack.GetNDF();
  Double_t chi2 = kaltrack.GetChi2();
  Double_t cl   = TMath::Prob(chi2, ndf);
  Double_t dr   = trkState(0, 0); 
  Double_t fi0  = trkState(1, 0); 
  Double_t cpa  = trkState(2, 0);
  Double_t dz   = trkState(3, 0);
  Double_t tnl  = trkState(4, 0); 
  Double_t cs   = tnl/TMath::Sqrt(1.+tnl*tnl);
  Double_t t0   = trkState(5, 0); 
  
//   svd(0,0) = 0.;                        // dr
//   svd(1,0) = helstart.GetPhi0();        // phi0
//   svd(2,0) = helstart.GetKappa();       // kappa
//   svd(3,0) = 0.;                        // dz
//   svd(4,0) = helstart.GetTanLambda();   // tan(lambda)

//    _FCT = 2.99792458E-4;
//     _radius = _pxy / (_FCT*B);
//     _omega = q/_radius;

  const double alpha =  2.99792458E-5 ; // r in cm !! 
  // 1/ kappa = p_t = alpha * B / omega 
  // omega = alpha * B * kappa

  if( cl > 0.3 ) 
    streamlog_out( DEBUG ) << " GOOD FIT !?  - Bfield "  << Bfield <<  std::endl ;
    


  //============== convert parameters to LCIO convention ====
  double phi = fi0 - M_PI/2. ;  // FIXME : which sign ?
  double omega =  - cpa * alpha *Bfield  ; 
  double d0 = -dr * 10. ;
  double z0 = dz * 10. ;
  double tanLambda =  -tnl ;
  //=========================================================

  streamlog_out( DEBUG ) << " kaltest track parameters: "
			 << " chi2/ndf " << chi2 / ndf  
    //			 << " chi2 " <<  chi2 << std::endl 
			 << " conv. level " << cl  
    			 << "\t D0 " <<  d0 
			 << "\t Phi :" << phi
			 << "\t Omega " <<  omega
			 << "\t Z0 " <<  z0
			 << "\t tan(Lambda) " <<  tanLambda 
			 << std::endl ;

  streamlog_out( DEBUG ) << " ----- PIVOT: " << ((TKalTrackSite&) cursite).GetPivot()(0) 
			 << ", " <<   ((TKalTrackSite&) cursite).GetPivot()(1) 
			 << ", " <<   ((TKalTrackSite&) cursite).GetPivot()(2) << std::endl ;
 


  float pivot[3] ;
  pivot[0] =  ((TKalTrackSite&) cursite).GetPivot()(0) * 10. ;
  pivot[1] =  ((TKalTrackSite&) cursite).GetPivot()(1) * 10. ;
  pivot[2] =  ((TKalTrackSite&) cursite).GetPivot()(2) * 10. ;


  //-------------- convert parameters to origin ---------------------
  // L3 Note 1666 "Helicoidal Tracks" - formulas 33, 44, 52 


  double xPr(0.),yPr(0.) ; // new reference point/pivot
  double dx = xPr - pivot[0]  ;
  double dy = yPr - pivot[1]  ;

  double phiPr = std::atan2(   std::sin( phi ) - ( dx / ( 1./omega - d0 ) ) ,
			       std::cos( phi ) + ( dy / ( 1./omega - d0 ) ) ) ;

  double d0Pr = d0 + dx * std::sin( phi ) - dy * std::cos( phi ) 
    - omega/( 2.*( 1.- omega*d0) )  * ( dx*dx +dy*dy)  ;



  double dPhi = phiPr - phi ;
  double s = dx  * std::cos( phi ) + dy * std::sin( phi ) / ( std::sin( dPhi ) / dPhi  ) ;

  double z0Pr = z0 + s * tanLambda ;

  trk->setD0( d0Pr ) ;  
  trk->setPhi( phiPr  ) ; // fi0  - M_PI/2.  ) ;  
  trk->setOmega( omega  ) ;
  trk->setZ0( dz  ) ;  
  trk->setTanLambda( tanLambda ) ;  


  //  trk->setReferencePoint( pivot ) ;
}
