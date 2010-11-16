//*************************************************************************
//* ======================
//*  EXTPCMeasLayer Class
//* ======================
//*
//* (Description)
//*   User defined measurement layer class
//* (Requires)
//*     EXVMeasLayer
//* (Provides)
//*     class EXTPCMeasLayer
//* (Update Recored)
//*   2010/09     F.Gaede - copied from KalDet
//*   2009/11/23  K.Ikematsu   Derived from KalTest/examples/kaltest/
//*                                         hybrid/tpc/EXTPCMeasLayer.cxx
//*
//* $Id: EXTPCMeasLayer.cxx,v 1.1.1.1 2009/11/24 00:13:59 ikematsu Exp $
//*************************************************************************
//

//#define USE_ZDRIFT 

#include "kaltest/TKalTrack.h"    // from KalTrackLib

#include "EXTPCMeasLayer.h"
#include "EXTPCHit.h"
#include "EXTPCKalDetector.h"
#include "TRandom.h"
#include "TMath.h"
#include <cmath>

#include "EVENT/TrackerHit.h"

#include "streamlog/streamlog.h"

//ClassImp(EXTPCMeasLayer)

EXTPCMeasLayer::EXTPCMeasLayer(TMaterial &min,
                               TMaterial &mout,
                               Double_t   r0,
                               Double_t   lhalf,
                               Double_t   sigmax0,
                               Double_t   sigmax1,
                               Double_t   sigmaz0,
                               Double_t   sigmaz1,
                               Bool_t     type,
			       int        layerID )
  : EXVMeasLayer(min, mout, type, layerID),
    TCylinder(r0, lhalf),
    fPhiMin(-TMath::Pi()/2),
    fPhiMax(+TMath::Pi()/2),
    fSigmaX0(sigmax0),
    fSigmaX1(sigmax1),
    fSigmaZ0(sigmaz0),
    fSigmaZ1(sigmaz1),
    fModule(-1),
    fLayer(-1), 
    _dEdx(0.0) 
{

  //fg: don't think we need phiMin/Max:
  fPhiMin = 0. ;
  fPhiMax = 0. ;

}

EXTPCMeasLayer::EXTPCMeasLayer(TMaterial &min,
                               TMaterial &mout,
                               Double_t   r0,
                               TVector3   xc,
                               Double_t   phimin,
                               Double_t   phimax,
                               Double_t   lhalf,
                               Double_t   sigmax0,
                               Double_t   sigmax1,
                               Double_t   sigmaz0,
                               Double_t   sigmaz1,
                               Bool_t     type,
                               Int_t      module,
                               Int_t      layer)
  : EXVMeasLayer(min, mout, type , -1 ),
    TCylinder(r0, lhalf, xc.X(), xc.Y(), xc.Z()),
    fPhiMin(phimin),
    fPhiMax(phimax),
    fSigmaX0(sigmax0),
    fSigmaX1(sigmax1),
    fSigmaZ0(sigmaz0),
    fSigmaZ1(sigmaz1),
    fModule(module),
    fLayer(layer),
    _dEdx(0.0) 
{
}

EXTPCMeasLayer::~EXTPCMeasLayer()
{
}

TKalMatrix EXTPCMeasLayer::XvToMv(const TVector3 &xv,
                                        Int_t     side) const
{
  TVector3 xxv = xv - GetXc();
#if 1
  Double_t phi = TMath::ATan2(xxv.Y(), xxv.X()) - fPhiMin;
#else
  Double_t phi = TMath::ATan2(xxv.Y(), xxv.X()) + fPhiMin;
#endif

  static Double_t kPi    = TMath::Pi();
  static Double_t kTwoPi = 2 * kPi;
  while (phi < -kPi) phi += kTwoPi;
  while (phi >  kPi) phi -= kTwoPi;

  // Calculate hit coordinate information:
  //   mv(0, 0) = r * phi
  //     (1, 0) = drift distance

  TKalMatrix mv(kMdim, 1);

  mv(0, 0) = GetR() * phi;

#ifdef USE_ZDrift 
  mv(1, 0) = (GetLength() * 0.5) - side * xxv.Z();
#else
  mv(1, 0) = xxv.Z();
#endif
  return mv;
}

TKalMatrix EXTPCMeasLayer::XvToMv(const TVTrackHit &vht,
                                  const TVector3   &xv) const
{
  return XvToMv(xv, dynamic_cast<const EXTPCHit &>(vht).GetSide());
}

TVector3 EXTPCMeasLayer::HitToXv(const TVTrackHit &vht) const
{
  const EXTPCHit &ht = dynamic_cast<const EXTPCHit &>(vht);

#if 1
  Double_t phi = ht(0, 0) / GetR() + fPhiMin;
#else
  Double_t phi = ht(0, 0) / GetR() - fPhiMin;
#endif

#ifdef USE_ZDrift 
  Double_t z   = ht.GetSide() * (GetLength() * 0.5 - ht(1, 0)) + GetXc().Z();
#else
  Double_t z   = ht(1, 0);
#endif
  Double_t x   = GetR() * TMath::Cos(phi) + GetXc().X();
  Double_t y   = GetR() * TMath::Sin(phi) + GetXc().Y();

  return TVector3(x, y, z);
}

void EXTPCMeasLayer::CalcDhDa(const TVTrackHit &vht,
                              const TVector3   &xxv,
                              const TKalMatrix &dxphiada,
                                    TKalMatrix &H) const
{
  const EXTPCHit &ht = dynamic_cast<const EXTPCHit &>(vht);

  // Calculate
  //    H = (@h/@a) = (@phi/@a, @z/@a)^t
  //  where
  //        h(a) = (phi, z)^t: expected meas vector
  //        a = (drho, phi0, kappa, dz, tanl, t0)
  //

  Int_t sdim = H.GetNcols();
  Int_t hdim = TMath::Max(5, sdim - 1);

  TVector3 xxvc = xxv - GetXc();
  Double_t xv   = xxvc.X();
  Double_t yv   = xxvc.Y();
  Double_t xxyy = xv * xv + yv * yv;

  // Set H = (@h/@a) = (@d/@a, @z/@a)^t

  for (Int_t i = 0; i < hdim; i++) {
    H(0, i)  = - (yv / xxyy) * dxphiada(0, i)
               + (xv / xxyy) * dxphiada(1, i);
    H(0, i) *= GetR();

#ifdef USE_ZDRIFT
    H(1, i)  = - ht.GetSide() * dxphiada(2, i); 
#else
    H(1, i)  = dxphiada(2, i);
#endif
  }
  if (sdim == 6) {
    H(0, sdim - 1) = 0.;

// #if 0
//     H(1, sdim - 1) = ht.GetVdrift();
// #else
//     H(1, sdim - 1) = - ht.GetVdrift();
// #endif
  }
}

Double_t EXTPCMeasLayer::GetSigmaX(Double_t zdrift) const
{
  return TMath::Sqrt(fSigmaX0 * fSigmaX0 + fSigmaX1 * fSigmaX1 * zdrift);
}

Double_t EXTPCMeasLayer::GetSigmaZ(Double_t zdrift) const
{
  return TMath::Sqrt(fSigmaZ0 * fSigmaZ0 + fSigmaZ1 * fSigmaZ1 * zdrift);
}

Double_t EXTPCMeasLayer::GetSortingPolicy() const
{
  Double_t value = TCylinder::GetSortingPolicy();
  if (fModule >= 0) {
#if 0
    ///// FIXME: temporary treatment //////////////////
    value = GetR() + GetXc().X(); // LP1 GEAR only!!
    ///////////////////////////////////////////////////
#endif
  }
  return value;
}

void EXTPCMeasLayer::ProcessHit(const TVector3  &xx,
                                      TObjArray &hits, 
				EVENT::TrackerHit* hit)
{

  static const double epsilon = 0.0001 ; // 1 micron 

  Int_t      side = (xx.Z() < 0. ? -1 : 1);
  TKalMatrix h    = XvToMv(xx, side);
  Double_t   rphi = h(0, 0);
  Double_t   d    = h(1, 0);

  //fg: we need zdrift here ! 
  Double_t zDrift   =  GetLength() * 0.5 - std::abs( d )  ;


  // Double_t dx = GetSigmaX( zDrift );
  // Double_t dz = GetSigmaZ( zDrift );

  //   Double_t dx = GetSigmaX(d);
  //   Double_t dz = GetSigmaZ(d);

  // use errors stored in LCIO hit
  double dx = sqrt( hit->getCovMatrix()[0] + hit->getCovMatrix()[2] );
  double dz = sqrt( hit->getCovMatrix()[5] );

  Double_t v = EXTPCKalDetector::GetVdrift();

  //if( streamlog_level( DEBUG ) && std::abs( xx.Perp() - GetR() ) > epsilon ) {
  if( std::abs( xx.Perp() - GetR() ) > epsilon ) {
    
    streamlog_out( ERROR ) << " TPC hit at r = " << xx.Perp() 
			   << " not on measurement layer R = " << GetR() << std::endl ; 
    
  }


  Double_t meas [2];
  Double_t dmeas[2];
  meas [0] = rphi;
  meas [1] = d;
  dmeas[0] = dx;
  dmeas[1] = dz;


  Double_t b = EXTPCKalDetector::GetBfield();
  hits.Add(new EXTPCHit(*this, meas, dmeas, side, v , hit, b  ));

  // TVector3 pv = dynamic_cast<EXTPCHit*>(  hits.Last() )->GetExactX() ; 
  // streamlog_out( DEBUG ) << " EXTPCHit hit[" << hits.Last() << "]at : " << pv[0] << ", "  << pv[1] << ", "  << pv[2] << std::endl ;
}

void  EXTPCMeasLayer::addIPHit(const TVector3   &xx,
			       TObjArray  &hits) {
  
  
  //  static const double epsilon = 0.0001 ; // 1 micron 
  static const double epsilon = 0.001 ; // 1 micron 
  
  Int_t      side = (xx.Z() < 0. ? -1 : 1);
  TKalMatrix h    = XvToMv(xx, side);
  Double_t   rphi = h(0, 0);
  Double_t   d    = h(1, 0);
  
  
  streamlog_out( DEBUG ) << " adding faked TPC hit at r = " << xx.Perp()  << std::endl ;
  
  
  if( std::abs( xx.Perp() - GetR() ) > epsilon ) {
  //if( streamlog_level( DEBUG ) && std::abs( xx.Perp() - GetR() ) > epsilon ) {
    
    streamlog_out( ERROR ) << " TPC hit at r = " << xx.Perp() 
			   << " not on measurement layer R = " << GetR() << std::endl ; 
    
  }

  Double_t meas [2];
  Double_t dmeas[2];
  meas [0] = rphi;
  meas [1] = d;
  dmeas[0] = 100. ; // large errors - should not change track state
  dmeas[1] = 100. ; //    
  
  Double_t v = EXTPCKalDetector::GetVdrift();

  Double_t b = EXTPCKalDetector::GetBfield();

  hits.Add(new EXTPCHit(*this, meas, dmeas, side, v, xx, b));

 }


Double_t EXTPCMeasLayer::GetEnergyLoss( Bool_t    isoutgoing,
					const TVTrack  &hel,
					Double_t  df) const
{
   Double_t cpa    = hel.GetKappa();
   Double_t tnl    = hel.GetTanLambda(); 
   Double_t tnl2   = tnl * tnl;
   Double_t tnl21  = 1. + tnl2;
   Double_t cslinv = TMath::Sqrt(tnl21);
   Double_t mom2   = tnl21 / (cpa * cpa);

   // // -----------------------------------------
   // // Bethe-Bloch eq. (Physical Review D P195.)
   // // -----------------------------------------
   // static const Double_t kK   = 0.307075e-3;     // [GeV*cm^2]
   //static const Double_t kMe  = 0.510998902e-3;  // electron mass [GeV]
   static const Double_t kMpi = 0.13957018;      // pion mass [GeV]

   TKalTrack *ktp  = static_cast<TKalTrack *>(TVKalSystem::GetCurInstancePtr());
   Double_t   mass = ktp ? ktp->GetMass() : kMpi;

   // const TMaterial &mat = GetMaterial(isoutgoing);
   // Double_t dnsty = mat.GetDensity();		// density
   // Double_t A     = mat.GetA();                 // atomic mass
   // Double_t Z     = mat.GetZ();                 // atomic number
   // //Double_t I    = Z * 1.e-8;			// mean excitation energy [GeV]
   // //Double_t I    = (2.4 +Z) * 1.e-8;		// mean excitation energy [GeV]
   // Double_t I    = (9.76 * Z + 58.8 * TMath::Power(Z, -0.19)) * 1.e-9;
   // Double_t hwp  = 28.816 * TMath::Sqrt(dnsty * Z/A) * 1.e-9;
   // Double_t bg2  = mom2 / (mass * mass);
   // Double_t gm2  = 1. + bg2;
   // Double_t meM  = kMe / mass;
   // Double_t x    = log10(TMath::Sqrt(bg2));
   // Double_t C0   = - (2. * log(I/hwp) + 1.);
   // Double_t a    = -C0/27.;
   // Double_t del;
   // if (x >= 3.)            del = 4.606 * x + C0;
   // else if (0.<=x && x<3.) del = 4.606 * x + C0 + a * TMath::Power(3.-x, 3.);
   // else                    del = 0.;
   // Double_t tmax = 2.*kMe*bg2 / (1. + meM*(2.*TMath::Sqrt(gm2) + meM)); 
   // Double_t dedx = kK * Z/A * gm2/bg2 * (0.5*log(2.*kMe*bg2*tmax / (I*I))
   //               - bg2/gm2 - del);

   Double_t path = hel.IsInB()
                 ? TMath::Abs(hel.GetRho()*df)*cslinv
                 : TMath::Abs(df)*cslinv;
   
   //Double_t edep = dedx * dnsty * path;

   Double_t edep = _dEdx * path;

   streamlog_out( DEBUG ) << " dEdx - energy loss " <<  edep << " path lengths: " << path << " dEdx: " << _dEdx << std::endl ;


   Double_t cpaa = TMath::Sqrt(tnl21 / (mom2 + edep
                 * (edep + 2. * TMath::Sqrt(mom2 + mass * mass))));
   Double_t dcpa = TMath::Abs(cpa) - cpaa;

   static const Bool_t kForward  = kTRUE;
   static const Bool_t kBackward = kFALSE;
   Bool_t isfwd = ((cpa > 0 && df < 0) || (cpa <= 0 && df > 0)) ? kForward : kBackward;
   return isfwd ? (cpa > 0 ? dcpa : -dcpa) : (cpa > 0 ? -dcpa : dcpa);
}
