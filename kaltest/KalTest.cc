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
#include "gear/TPCParameters.h"
#include "gear/PadRowLayout2D.h"


#include "streamlog/streamlog.h"


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
  EXTPCKalDetector* tpcdet = new EXTPCKalDetector ;
  
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

void KalTest::fitTrack(){

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

  // ---------------------------
  //  Smooth the track
  // ---------------------------
#ifndef SAVE_RESIDUAL
  TVKalSite &cursite = kaltrack.GetCurSite();
#else
  Int_t isite = 1;
  kaltrack.SmoothBackTo(isite);
  TVKalSite &cursite = static_cast<TVKalSite &>(*kaltrack[isite]);
#endif


  // dump fit result for now:
  Int_t    ndf  = kaltrack.GetNDF();
  Double_t chi2 = kaltrack.GetChi2();
  Double_t cl   = TMath::Prob(chi2, ndf);
  Double_t fi0  = cursite.GetCurState()(1, 0); 
  Double_t cpa  = cursite.GetCurState()(2, 0); 
  Double_t tnl  = cursite.GetCurState()(4, 0); 
  Double_t cs   = tnl/TMath::Sqrt(1.+tnl*tnl);
  Double_t t0   = cursite.GetCurState()(5, 0); 
  
  streamlog_out( DEBUG ) << " track parameters: " << std::endl 
			 << " ndf " << ndf  << std::endl 
			 << " chi2 " <<  chi2 << std::endl 
			 << " cl " << cl  << std::endl 
			 << " fi0 " <<  fi0 << std::endl 
			 << " cpa " <<  cpa << std::endl 
			 << " tnl " <<  tnl << std::endl 
			 << " cs " <<  cs << std::endl 
			 << " t0 " <<  t0 << std::endl  ;

}
