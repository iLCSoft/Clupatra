#include "KalTrack.h"

//#include "LCIOTTree.h" 

#include "kaltest/TKalTrackState.h"
#include "kaltest/TKalTrackSite.h"
#include "kaltest/TVTrackHit.h"
#include "kaltest/TKalDetCradle.h"
#include "kaltest/TKalTrack.h"         // from KalTrackLib

#include "TKalFilterCond.h"

#include "EXVMeasLayer.h"
#include "EXVTXKalDetector.h"
#include "EXTPCKalDetector.h"
#include "EXVTXHit.h"
#include "EXTPCHit.h"
//#include "EXHYBTrack.h"

#include "gear/GEAR.h"
#include "gear/BField.h"
#include "gear/TPCParameters.h"
#include "gear/PadRowLayout2D.h"

#include "IMPL/TrackImpl.h"
#include "IMPL/TrackerHitImpl.h"
#include "UTIL/Operators.h"
#include "lcio.h"

#include "clupatra.h"

#include "TObjArray.h"

#include <math.h>
#include <cmath>

#include "streamlog/streamlog.h"



//---------------------------------------------------------------------------------------------------------------

/** Helper class for defining a filter condition based on the delta chi2 in the AddAndFilter step.
 */
class KalTrackFilter : public TKalFilterCond{

public:
  
  /** C'tor - takes as optional argument the maximum allowed delta chi2 for adding the hit (in IsAccepted() )
   */
  KalTrackFilter(double maxDeltaChi2 = -1.) : _deltaChi2( 0.0 ) , 
					      _maxDeltaChi2( maxDeltaChi2 ) {
  } 
  virtual ~KalTrackFilter() {} 
  
  double deltaChi2() { return _deltaChi2 ; }
  
  virtual Bool_t IsAccepted(const TKalTrackSite &site) {
    
    _deltaChi2 = site.GetDeltaChi2();
    
    streamlog_out( DEBUG0 ) << " KalTrackFilter::IsAccepted called  !  deltaChi2 = "  <<  _deltaChi2  << std::endl;

    return ( _maxDeltaChi2 > -1. ?  _deltaChi2 < _maxDeltaChi2  : true )   ; 
  }

protected:

  double _deltaChi2 ;
  double _maxDeltaChi2 ;

} ;
//---------------------------------------------------------------------------------------------------------------


/** helper function to restrict the range of the azimuthal angle to ]-pi,pi]*/
inline double toBaseRange( double phi){
  while( phi <= -M_PI ){  phi += 2. * M_PI ; }
  while( phi >   M_PI ){  phi -= 2. * M_PI ; }
  return phi ;
}



std::ostream& operator<<(std::ostream& o, const KalTrack& trk) {

  o << " track: \t" <<  trk._trk->GetEntriesFast() << std::endl ;
  // to be done ....
  return o ;
}


/** C'tor - initiale with detector */
KalTrack::KalTrack(TKalDetCradle* det, bool ownHits) : _det( det) , _xingPts( _det->GetEntriesFast() ) {
  _trk = new TKalTrack ;
  _trk->SetOwner();

  _kalHits = new TObjArray ;

  // if( ownHits )
  //   _kalHits->SetOwner() ;

  _kalHits->SetOwner( ownHits ) ;

  streamlog_out( DEBUG ) << " creating Kaltrack  - owning hits = " << ownHits << std::endl ;

  for( unsigned i=0 ; i< _xingPts.size() ; _xingPts[i++] = 0  ) ;
}

KalTrack::~KalTrack(){
  delete _trk ;
  delete _kalHits ;
  for( unsigned i=0 ; i< _xingPts.size() ; delete _xingPts[i++] ) ;
}




void KalTrack::addIPHit(){

  // add a faked Hit for the IP 
  TObject* o =  _det->At( 0 ) ;
  EXTPCMeasLayer* ml = dynamic_cast< EXTPCMeasLayer * >( o ) ;
  
  ml->addIPHit( TVector3( 1.2, 0., 0.) ,   *_kalHits ) ;
}

void KalTrack::addHit( const TVector3& pos, int layer , EVENT::TrackerHit* hit) {
  
  if( layer >= 0 && ( _det->GetEntries() > layer ) ) {

    TObject* o =  _det->At( layer ) ;

    EXVMeasLayer* ml = dynamic_cast< EXVMeasLayer * >( o ) ;

    if (ml != 0 && ml->IsActive() && dynamic_cast<const EXVKalDetector &>( ml->GetParent(kFALSE) ).IsPowerOn() ) {

      ml->ProcessHit( pos, *_kalHits , hit ); // create hit point
    

      if( streamlog_level( DEBUG )  &&  layer % 10 == 0 ){
	double radius = pos.Perp() ;
	streamlog_out( DEBUG )  << " ***** adding a TPC hit in layer : [" << layer <<  "] at R = " 
				<< radius << " id: " <<  hit->id()  <<std::endl ;
      }
    

    } else {
    
      streamlog_out( WARNING ) << " hit not added to KalTrack at : " 
			       << pos(0) << "," << pos(1) << "," << pos(2)  
			       << "  in  layer: " << layer
			       << " ml : " << ml 
			       << "  ml->IsActive() " <<  ml->IsActive()
			       << "detector->IsPowerOn() " 
			       << dynamic_cast<const EXVKalDetector &>( ml->GetParent(kFALSE) ).IsPowerOn()
			       << std::endl ;
    } 
  
  } else {
  
    streamlog_out( ERROR ) << " no measurement layer at : " 
			   << pos(0) << "," << pos(1) << "," << pos(2)  
			   << "  in  layer: " << layer 
			   << " detector has  " << _det->GetEntriesFast() << " layers only "
			   << " - object at given layer : " <<  _det->At( layer )  << std::endl ;
  }
}

unsigned KalTrack::getNHits() const { 
  return _kalHits->GetEntriesFast() ; 
}

void KalTrack::fitTrack( bool fitDirection ) {
  
  const Bool_t gkDir = fitDirection ; 
  
  // ============================================================
  //  Do Kalman Filter    - copied from EXVKalTrack.cxx
  // ============================================================
  
  if (_kalHits->GetEntries() < 3) {
    
    streamlog_out( ERROR) << "<<<<<< KalTrack::fitTrack(): Shortage of Hits! nhits = "  
			  << _kalHits->GetEntries() << " >>>>>>>" << std::endl;
    return ;
  }
  
  Int_t i1, i2, i3; // (i1,i2,i3) = (1st,mid,last) hit to filter
  if (gkDir == kIterBackward) {
    i3 = 0 ; // fg: first index is 0 and not 1 
    i1 = _kalHits->GetEntries() - 1;
    i2 = i1 / 2;
  } else {
    i1 = 0 ; 
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

  //  double Bfield =   h1.GetBfield() ;

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
    //    C(i,i) = 1.e6;   // dummy error matrix

    // NB: if the error is too large the initial helix parameters might be changed extremely by the first three (or so) hits,
    //     such that the fit will not work because the helix curls away and does not hit the next layer !!!
    //     ->  this might need a bit more thought ...
    C(i,i) = 1.e2;   // dummy error matrix
  }

  sited.Add(new TKalTrackState(svd,C,sited,TVKalSite::kPredicted));
  sited.Add(new TKalTrackState(svd,C,sited,TVKalSite::kFiltered));

  // ---------------------------
  //  Add sited to the kaltrack
  // ---------------------------

  //  EXHYBTrack kaltrack;   // a track is a kal system
  TKalTrack& kaltrack = *_trk ;

  kaltrack.SetOwner();   // kaltrack owns sites
  kaltrack.Add(&sited);  // add the dummy site to this track

  streamlog_out( DEBUG0 )  << "Kaltrack::fitTrack :  add dummy site at index : " << htdp->GetMeasLayer().GetIndex() << std::endl ;


  // ---------------------------
  //  Prepare hit iterrator
  // ---------------------------

  TIter next(_kalHits, gkDir); // come in to IP, if gkDir = kIterBackward

  // ---------------------------
  //  Start Kalman Filter
  // ---------------------------

  TVTrackHit *hitp = 0;

  KalTrackFilter filter( 35. )   ; //FIXME: make parameter ?
 
  int counter = -1 ;
  while ( (hitp = dynamic_cast<TVTrackHit *>( next() ) ) ) {
    
    TKalTrackSite  &site = *new TKalTrackSite(*hitp); // new site

    ++counter ;
    if( counter > 4  )   // allow first hits to have large dChi2 ....
      site.SetFilterCond( &filter ) ;
    
    streamlog_out( DEBUG0 )  << "Kaltrack::fitTrack :  add site at index : " << hitp->GetMeasLayer().GetIndex() << std::endl ;


    if ( ! kaltrack.AddAndFilter(site)  ) {               // filter it

      streamlog_out( DEBUG4 )  << "Kaltrack::fitTrack :  site discarded!" 
			       << "  cluster : " << getCluster< clupatra::GCluster >() << " size " << getCluster< clupatra::GCluster >()->size()  
			       << std::endl;
      //FIXME:   should we remove the hit from the cluster here ???????

      delete &site;                        // delete it if failed
    }
  } // end of Kalman filter

#define SMOOTH_BACK 0
#if SMOOTH_BACK
  //  Int_t isite = 1;
  //  kaltrack.SmoothBackTo(isite);
  kaltrack.SmoothAll();
#endif
  

  //  double bla = chi2( *this , *this ) ;

}

bool KalTrack::addAndFilter( TVTrackHit* hit, double maxDeltaChi2){
  
  KalTrackFilter filter( maxDeltaChi2 )   ;

  TKalTrackSite* site = new TKalTrackSite( *hit ); // new site

  site->SetFilterCond( &filter ) ;

  bool success = _trk->AddAndFilter( *site ) ;
  
  if( ! success ) {
    
    streamlog_out( DEBUG4 )  << "Kaltrack::addAndFilter :  site discarded!    - delta Chi2: " <<  filter.deltaChi2()    << std::endl;
    
    delete site;                        // delete it if failed
    
    return false ;
  }
  
  _kalHits->Add( hit ) ;
  
  return true ;
}

double KalTrack::testDeltaChi2( TVTrackHit* hit ){

  KalTrackFilter filter( -0.0 )   ; // filter will always fail

  TKalTrackSite* site = new TKalTrackSite( *hit ); // new site

  site->SetFilterCond( &filter ) ;

  bool success = _trk->AddAndFilter( *site ) ;
  
  if(  success ) {
    
    streamlog_out( ERROR )  << "Kaltrack::testDeltaChi2 :  site accepted    - delta Chi2: " <<  filter.deltaChi2()    << std::endl;
    
  } else {
    
    delete site;                        // delete it if failed
  }    

  return  filter.deltaChi2()  ;
}




double KalTrack::chi2( const KalTrack& t0 , const KalTrack& t1) {

 const TKalTrackState& ts0 = t0.getTrackState() ; 
 const TKalTrackState& ts1 = t1.getTrackState() ; 

 double dPhi ;
 const TMatrixD& c0 =  ts0.GetCovMat() ;
 TMatrixD cov0(5,5) ;  for(int i=0;i<5;++i) for(int j=0;j<5;++j) cov0[i][j] = c0[i][j] ;

 const TMatrixD& c1  = ts1.GetCovMat() ;
 TMatrixD cov1(5,5) ;  for(int i=0;i<5;++i) for(int j=0;j<5;++j) cov1[i][j] = c1[i][j] ;

 // cov0.Print() ;

 THelicalTrack h0 = ts0.GetHelix() ;
 h0.MoveTo(  TVector3( 0., 0., 0. ) , dPhi , 0 , &cov0 ) ;

 //cov0.Print() ;

 THelicalTrack h1 = ts1.GetHelix() ;
 h1.MoveTo(  TVector3( 0., 0., 0. ) , dPhi , 0 , &cov1 ) ;

  double tp0[5] , tp1[5] ;
  tp0[0] =  - h0.GetDrho() ; 
  tp0[1] =    toBaseRange( h0.GetPhi0() + M_PI/2. ) ;
  tp0[2] =    1. /h0.GetRho()  ;              
  tp0[3] =    h0.GetDz()   ;
  tp0[4] =    h0.GetTanLambda()  ;

  tp1[0] =  - h1.GetDrho() ; 
  tp1[1] =    toBaseRange( h1.GetPhi0() + M_PI/2. ) ;
  tp1[2] =    1. /h1.GetRho()  ;              
  tp1[3] =    h1.GetDz()   ;
  tp1[4] =    h1.GetTanLambda()  ;

  //adjust z0 to phi<2Pi :
  double& z0 = tp0[3] ;
  double r0 = std::abs( h0.GetRho() ) ;
  double sz0 =  2.* M_PI * r0 * std::abs( tp0[4] ) ;
  if( z0 > 0.  ){
    //    streamlog_out( DEBUG4 ) << " adjusting z0 : " << z0 << " sz0 : " << sz0 << " r: " << r0 << std::endl ; 
    while( z0 > sz0 )  z0 -= sz0 ;
    //streamlog_out( DEBUG4 ) << " adjusted z0 : " << z0 << " sz0 : " << sz0 << std::endl ; 
  }else{
    while( z0 < -sz0 ) z0 += sz0 ;
  }
  double& z1 = tp1[3] ;
  double r1 = std::abs( h1.GetRho() ) ;
  double sz1 =  2.* M_PI * r1 * std::abs( tp1[4] ) ;
  if( z1 > 0.  ){
    //streamlog_out( DEBUG4 ) << " adjusting z1 : " << z1 << " sz1 : " << sz1 << std::endl ; 
    while( z1 > sz1 )  z1 -= sz1 ;
    // streamlog_out( DEBUG4 ) << " adjusting z1 : " << z1 << " sz1 : " << sz1 << std::endl ; 
  }else{
    while( z1 < -sz1 ) z1 += sz1 ;
  }
    
  // adjust kappa cov parameter for omega
  double al2 = ( 1. / h0.GetRho() ) / h0.GetKappa() ;
  al2 *= al2 ;
  cov0( 2, 2 ) *= al2 ;
  cov1( 2, 2 ) *= al2 ;
  
  // add errors in quadrature
  TKalMatrix cov  = ( cov0 + cov1 ) ;
  
  double chi2( 0.0 ) ;
  
  for(int i=0; i<5 ; ++i ) {
    
    double diff2 = ( std::abs( tp0[i] )- std::abs( tp1[i] )  ) ; 
    
    diff2 *= diff2 ;
    
    streamlog_out( DEBUG )  << " -------    diff : " << i <<  " : " << tp0[i] <<" - " <<  tp1[i] << " = " <<  diff2  
			    << " - cov() " << cov(i,i)  <<  " diff2/cov(i,i) " << diff2/cov(i,i)  << std::endl ;
    
    chi2 +=  ( diff2 / cov( i , i ) )  ;
    
  } 
  
  return chi2 ;
  
  
  // ===== compute 'proper' chi2  ===== DOESN'T WORK :(
  // TKalMatrix cov( ts0.GetCovMat().GetSub(0,4,0,4) ) ; 
  // TKalMatrix covInv( TKalMatrix::kInverted , cov ) ;
  // TKalMatrix s0( ts0.GetSub(0,4,0,0) );
  // TKalMatrix s1( ts1.GetSub(0,4,0,0) );
  // // s0.DebugPrint() ;
  // // s1.DebugPrint() ;
  // s1 -= s0 ;
  // TKalMatrix d0( covInv , TKalMatrix::kMult , s1 ) ;
  // TKalMatrix d1( s1 , TKalMatrix::kTransposeMult , d0 ) ;
  // // covInv.DebugPrint() ;
  // // d0.DebugPrint() ;
  // // d1.DebugPrint() ;
  // return  d1(0,0)  ; 
}
 

TKalTrackState& KalTrack::getTrackState() const {
  
  TVKalSite& cursite = _trk->GetCurSite();
  
  return (TKalTrackState&) cursite.GetCurState() ;
}

double KalTrack::getOmega() { 

  return 1./ getTrackState().GetHelix().GetRho() ; 
} 

 
void KalTrack::toLCIOTrack( IMPL::TrackImpl* trk) {

  TKalTrack& kaltrack = *_trk ;
  
  trk->ext< KalTrackLink >() = this ; // link this KalTrack to LCIO track 

  TKalTrackState& trkState = getTrackState() ;

  TVKalSite &cursite = kaltrack.GetCurSite();

  Int_t    ndf  = kaltrack.GetNDF();
  Double_t chi2 = kaltrack.GetChi2();
  Double_t cl   = TMath::Prob(chi2, ndf);
  
  //   Double_t dr   = trkState(0, 0); 
  //   Double_t fi0  = trkState(1, 0); 
  //   Double_t cpa  = trkState(2, 0);
  //   Double_t dz   = trkState(3, 0);
  //   Double_t tnl  = trkState(4, 0); 
  //   Double_t cs   = tnl/TMath::Sqrt(1.+tnl*tnl);
  //   Double_t t0   = trkState(5, 0); 
  
  //============== convert parameters to LCIO convention ====
  
  //  ---- get parameters at origin 
  
  THelicalTrack helix = trkState.GetHelix() ;

  double dPhi ;

  // need to get the 5x5 sub matrix of the covariance matrix  
  const TMatrixD& c0 =  trkState.GetCovMat() ;
  TMatrixD covK(5,5) ;  for(int i=0;i<5;++i) for(int j=0;j<5;++j) covK[i][j] = c0[i][j] ;
  
  helix.MoveTo(  TVector3( 0., 0., 0. ) , dPhi , 0 , &covK ) ;

  //-------------------------------------

  // double phi       =    toBaseRange( helix.GetPhi0()  - M_PI/2. ) ;
  // double omega     =   -1. /helix.GetRho()  ;              
  // double d0        =    helix.GetDrho() ; 
  // double z0        =    helix.GetDz()   ;
  // double tanLambda =   -helix.GetTanLambda()  ;

  //  this is for incomming tracks ...
  double phi       =    toBaseRange( helix.GetPhi0() + M_PI/2. ) ;
  double omega     =    1. /helix.GetRho()  ;              
  double d0        =  - helix.GetDrho() ; 
  double z0        =    helix.GetDz()   ;
  double tanLambda =    helix.GetTanLambda()  ;

  //=========================================================

  trk->setD0( d0 ) ;  
  trk->setPhi( phi  ) ; // fi0  - M_PI/2.  ) ;  
  trk->setOmega( omega  ) ;
  trk->setZ0( z0  ) ;  
  trk->setTanLambda( tanLambda ) ;  

  trk->subdetectorHitNumbers().push_back( 1 ) ;  // workaround for bug in lcio::operator<<( Tracks ) - used for picking ....

  trk->setChi2( chi2 ) ;

  //covariance matrix in LCIO - stored as lower triangle matrix where the order of parameters is: 
  // d0, phi, omega, z0, tan(lambda). So we have cov(d0,d0), cov( phi, d0 ), cov( phi, phi), ...

  //   Double_t dr   = trkState(0, 0);   
  //   Double_t fi0  = trkState(1, 0); 
  Double_t cpa  = trkState(2, 0);
  //   Double_t dz   = trkState(3, 0);
  //   Double_t tnl  = trkState(4, 0); 

  double alpha = omega / cpa  ;

  //  const TKalMatrix& covK = trkState.GetCovMat() ; 
  
  if( streamlog_level( DEBUG ) ) {
    streamlog_out( DEBUG ) << " KalTrack::toLCIOTrack : returning covariance matrix :  - alpha : " << alpha << std::endl ;
    //    covK.DebugPrint() ;
  }

  EVENT::FloatVec cov( 15 )  ; 
  cov[ 0] =   covK( 0 , 0 )   ; //   d0,   d0

  cov[ 1] = - covK( 1 , 0 )   ; //   phi0, d0
  cov[ 2] =   covK( 1 , 1 )   ; //   phi0, phi

  cov[ 3] = - covK( 2 , 0 ) * alpha   ; //   omega, d0
  cov[ 4] =   covK( 2 , 1 ) * alpha   ; //   omega, phi
  cov[ 5] =   covK( 2 , 2 ) * alpha * alpha  ; //   omega, omega

  cov[ 6] = - covK( 3 , 0 )   ; //   z0  , d0
  cov[ 7] =   covK( 3 , 1 )   ; //   z0  , phi
  cov[ 8] =   covK( 3 , 2 ) * alpha   ; //   z0  , omega
  cov[ 9] =   covK( 3 , 3 )   ; //   z0  , z0

  cov[10] = - covK( 4 , 0 )   ; //   tanl, d0
  cov[11] =   covK( 4 , 1 )   ; //   tanl, phi
  cov[12] =   covK( 4 , 2 ) * alpha  ; //   tanl, omega
  cov[13] =   covK( 4 , 3 )   ; //   tanl, z0
  cov[14] =   covK( 4 , 4 )   ; //   tanl, tanl

 
  trk->setCovMatrix( cov ) ;

  //fixme: shouldn't we set the origin here ??
  float pivot[3] ;
  pivot[0] =  ((TKalTrackSite&) cursite).GetPivot()(0) ;
  pivot[1] =  ((TKalTrackSite&) cursite).GetPivot()(1) ;
  pivot[2] =  ((TKalTrackSite&) cursite).GetPivot()(2) ;

  trk->setReferencePoint( pivot ) ;

  streamlog_out( DEBUG ) << " kaltest track parameters: "
			 << " chi2/ndf " << chi2 / ndf  
    			 << " chi2 " <<  chi2 << std::endl 
			 << " conv. level " << cl  

    			 << "\t D0 "          <<  d0         <<  "[+/-" << sqrt( cov[0] ) << "] " 
			 << "\t Phi :"        <<  phi        <<  "[+/-" << sqrt( cov[2] ) << "] " 
			 << "\t Omega "       <<  omega      <<  "[+/-" << sqrt( cov[5] ) << "] " 
			 << "\t Z0 "          <<  z0         <<  "[+/-" << sqrt( cov[9] ) << "] " 
			 << "\t tan(Lambda) " <<  tanLambda  <<  "[+/-" << sqrt( cov[14]) << "] " 

			 << "\t pivot : [" << pivot[0] << ", " << pivot[1] << ", "  << pivot[2] 
			 << " - r: " << std::sqrt( pivot[0]*pivot[0]+pivot[1]*pivot[1] ) << "]" 
			 << std::endl ;

  trk->setNdf( ndf ) ;


  streamlog_out( DEBUG ) << lcio::header( *trk ) << std::endl ;
  streamlog_out( DEBUG ) << lcio::lcshort( (lcio::Track*)trk ) << std::endl ;



#define ADD_LCIO_HITS
#ifdef ADD_LCIO_HITS

  //==========================  add lcio hits to the track =========================

  int nHit = kaltrack.GetEntriesFast() ;

  streamlog_out( DEBUG ) <<  " ============================= kaltrack.GetEntriesFast() " << kaltrack.GetEntriesFast()   << std::endl ;

  EVENT::TrackerHit* hLast = 0 ; 

  for(int i=0 ; i < nHit ; ++i ){

    EVENT::TrackerHit* h = const_cast<EVENT::TrackerHit*> 
      (static_cast<const EVENT::TrackerHit*>
       (  static_cast<const EXTPCHit&>
	  (dynamic_cast<TKalTrackSite &> ( *kaltrack[i] ).GetHit() ).GetHitPtr() ) ) ;

    if( h   && h != hLast ) {// protect against duplicate hits (from dummy site in fitTrack() ....)  

      trk->addHit( h ) ;
      hLast = h ;

    } else {
      streamlog_out( DEBUG ) <<  " hit [" << i << "] not added to LCIOTrack !!! " << std::endl ;
    }

  }

  if( nHit > ( trk->getTrackerHits().size()  + 1 ) )
    streamlog_out( DEBUG3 ) <<  "  hits not used in KalTest fit : " <<   nHit - trk->getTrackerHits().size() << " from " << nHit << std::endl ;


  if(  trk->getTrackerHits().size() < 6 ) {

    streamlog_out( ERROR ) <<  " Small number of hits in lcio Track " <<   trk->getTrackerHits().size() << " out of " << nHit  << std::endl ;
    

  }

#endif 

  //#define ADD_XING_HITS
#ifdef ADD_XING_HITS
  //==========================  add xing hits to the track for visualization - debugging .... =========================


  for( PointList::const_iterator it = _xingPts.begin() ; it != _xingPts.end() ; ++it ){
    
    if( ! *it ) continue ;

    IMPL::TrackerHitImpl* h = new IMPL::TrackerHitImpl ;  //memory leak - only use for debugging ....
    
    //    const gear::Vector3D& pv = it->second ;
    const gear::Vector3D& pv = **it ;

    double pos[3] ;
    pos[0] = pv[0] ;
    pos[1] = pv[1] ;
    pos[2] = pv[2] ;
    
    h->setPosition( pos ) ;

    //    streamlog_out( DEBUG ) << " KalTrack adding hit at : " << pos[0] << ", "  << pos[1] << ", "  << pos[2] << std::endl ;
    
    trk->addHit( h ) ;
    
  }
#endif 


  //#define ADD_DEBUG_HITS
#ifdef ADD_DEBUG_HITS
  //==========================  add some hits along the helix to the track for visualization - debugging .... ==============
  
  int nHit = 200 ;
  
  double phiStart = 0 ;
  double phiEnd = M_PI/4 ;
  
  double deltaPhi =  ( phiEnd - phiStart ) / nHit ;
  
  if( omega > 0 ) deltaPhi *= -1 ;
  
  for(int i=0 ; i < nHit ; ++i ){
    
    IMPL::TrackerHitImpl* h = new IMPL::TrackerHitImpl ;
    
    TVector3 pv = helix.CalcXAt( phiStart ) ; 
    double pos[3] ;
    pos[0] = pv[0] ;
    pos[1] = pv[1] ;
    pos[2] = pv[2] ;
    // pos[0] = pv[0] *10. ;
    // pos[1] = pv[1] *10. ;
    // pos[2] = pv[2] *10. ;
    
    phiStart += deltaPhi ;
    
    h->setPosition( pos ) ;
    
    trk->addHit( h ) ;
    
  }

#endif

  // ============================================================================================================
}

void KalTrack::findNextXingPoint(gear::Vector3D& v, int& layer, int step, bool backward) {
  
  layer = -1 ; // default return value
  
  // sanity check - require at least 5 hits
  if( _trk->GetEntries() < 5 ) 
    return ;

  // when smoothing backward  the first three sites have a bad cov. matrix and fit 
  // this might change, if the initial track state is computed from a chi2 fit ...
  const int firstGoodIndex = 4 ;
  if( backward )
    _trk->SmoothBackTo( firstGoodIndex ) ;
  
  TKalTrackSite& site0 =  *dynamic_cast<TKalTrackSite*>( _trk->At( firstGoodIndex ) ) ; 
  TKalTrackSite& siteL =  *dynamic_cast<TKalTrackSite*>( _trk->Last() ) ;
  
  int  idx0 = site0.GetHit().GetMeasLayer().GetIndex()  ; 
  int  idxL = siteL.GetHit().GetMeasLayer().GetIndex()  ; 


  int idx = (backward ? idx0 + firstGoodIndex : idxL  ) ;
  TKalTrackSite* site = (backward ? &site0 : &siteL ) ;

  bool isIncoming  = (backward ? idxL > idx0 :  idxL < idx0  ) ; 
  
  streamlog_out( DEBUG2 ) << " KalTrack::findNextXingPoint : " 
			  << " index at track site[0] : "  <<  idx0  
			  << " index at last site :     "  <<  idxL  << std::endl 
			  << "  site0: pivot : " << gear::Vector3D( site0.GetPivot() )
			  << "  siteL: pivot : " << gear::Vector3D( siteL.GetPivot() )
			  << "  isIncoming: " << isIncoming 
			  <<   std::endl ;
  
  std::auto_ptr<TVTrack> help(& static_cast<TKalTrackState &>( site->GetCurState() ).CreateTrack() ); // tmp track
  
  TVector3 xx ;       // expected hit position vector
  double  fid  = 0. ; // deflection angle from the last hit
  
  
  if  (isIncoming )  idx -= step ;     
  else               idx += step ;
  


  int lastIdx =  _det->GetEntriesFast() ;  // need # of tracking layers
  
  
  if(  idx < 0 || idx >= lastIdx ) return ;
  
  
  if (  dynamic_cast<TVSurface *>( _det->At( idx )  )->CalcXingPointWith( *help , xx, fid)  ){
    
    streamlog_out( DEBUG2  ) << " ---- next crossing at layer " << idx <<  " at: " 
			     << xx[0] << ", "  << xx[1] << ", " << xx[2] 
			     << " r: " << xx.Perp()  
			     << std::endl ;
    
    // return values
    layer = idx  ;

    v = gear::Vector3D( xx[0] , xx[1] , xx[2] )   ;
  }

}


void KalTrack::findXingPoints() {

  //  _xingPts.clear() ;
  
  int nSites = _trk->GetEntriesFast() -1 ;

  int nHits = _kalHits->GetEntriesFast()  ;


  if( nSites < 10 ){
    
    streamlog_out( DEBUG4 ) << "  ======= findXingPoints:   less then nSites used in fit:  No Crossing Points Computed !!! " 
			    << "  nSites: " << nSites
			    << std::endl ;
    
    return ; 
  }
  
  // ----- sanity check - don't compute xing points if some sites have been discarded in the fit
  if( nHits !=  nSites ) {

    streamlog_out( DEBUG4 ) << "  ======= findXingPoints:   hits discarded in fit : No Crossing Points Computed !!! " 
			    << "  nHits: " << nHits << " nSites : " << nSites
			    << std::endl ;
    return ;
  }
  
  // get first and last site used in track fit
  //TKalTrackSite& site0 =  *dynamic_cast<TKalTrackSite*>( _trk->At(1)  ) ; 
  //  TKalTrackSite& site10 =  *dynamic_cast<TKalTrackSite*>( _trk->At( 5 )  ) ; // use the 10-th site to have reasonable fit
  TKalTrackSite& siteL =  *dynamic_cast<TKalTrackSite*>( _trk->Last() ) ;


  //FIXME: smoothing back to site 1 does not work if initial cov matrix is just guessed (not fit) 
  // can only reasonably fit back to site 4 (first three needed for the cov. to 'adjust itself') 

  _trk->SmoothBackTo(1);
  TKalTrackSite& site0 =  *dynamic_cast<TKalTrackSite*>( _trk->At(1)  ) ; 
  
  int  idx0 = site0.GetHit().GetMeasLayer().GetIndex()  ; 
  int  idx1 = siteL.GetHit().GetMeasLayer().GetIndex()  ; 
  
  bool isIncoming  = idx1 < idx0 ; // the last site has the best track parameters 

  TKalTrackSite& siteIN =   ( isIncoming  ? siteL : site0 ) ; 
  TKalTrackSite& siteOUT =  ( isIncoming  ? site0 : siteL ) ; 


  streamlog_out( DEBUG4 ) << " KalTrack::findXingPoints : " 
			  << " index at track site[0] : "  <<  idx0  
			  << " index at last site :     "  <<  idx1  << std::endl 
			  << "  siteIN: pivot : " << gear::Vector3D( siteIN.GetPivot() )
			  << "  siteOUT: pivot : " << gear::Vector3D( siteOUT.GetPivot() )
			  <<   std::endl ;
  

  // search inwards first 
  int            idx  =  ( idx1 > idx0  ? idx0  : idx1  ) ;  
  //  TKalTrackSite& site =  ( idx1 > idx0  ? site0 : siteL ) ; 
  //  TKalTrackSite& site =  ( isIncoming  ? siteL : site10 ) ; 
  
  
  std::auto_ptr<TVTrack> help0(& static_cast<TKalTrackState &>( siteIN.GetCurState()).CreateTrack() ); // tmp track
  
  TVector3 xx ;       // expected hit position vector
  double  fid  = 0. ; // deflection angle from the last hit
  
  
  // if( ! isIncoming ){
  //   help0->MoveTo(  site0.GetPivot() , fid   ) ;
  // }
  

  --idx ;  // next site inwards

  while( idx >= 0 ) {

    

    if (  dynamic_cast<TVSurface *>( _det->At(idx)  )->CalcXingPointWith( *help0 , xx, fid)  ){
      
      streamlog_out( DEBUG4  ) << " ---- crossing layer " << idx <<  " at: " 
			     << xx[0] << ", "  << xx[1] << ", " << xx[2] 
			     << " r: " << xx.Perp()  
			     << std::endl ;

      _xingPts[ idx ] = new  gear::Vector3D( xx[0] , xx[1] , xx[2] )   ;

      --idx ;
      
    } else {

      streamlog_out( DEBUG4 ) << " ---- no crossing point found with layer " << idx << std::endl ;      

      break ;
    }
  }

  // now search outwards
  idx  =  ( idx1 > idx0  ? idx1  : idx0  ) ; 
  //site =  ( idx1 > idx0  ? siteL : site0 ) ; 
  //site =  ( isIncoming  ? site10 : siteL ) ; 


  streamlog_out( DEBUG4 ) << " ---- creating track for site at layer : " << idx 
			 << " state " << & siteOUT.GetCurState()
			 << std::endl ;
  
  std::auto_ptr<TVTrack> help(& static_cast<TKalTrackState &>( siteOUT.GetCurState()).CreateTrack() ); // tmp track
  
  // if(  isIncoming ){
  //   help->MoveTo(  site0.GetPivot() , fid   ) ;
  // }
  
  ++idx ;  // next site 

  int lastIdx =  _det->GetEntriesFast() ;  // need # of tracking layers

  while( idx < lastIdx ) { 


    TVector3 xx ;       // expected hit position vector
    double  fid  = 0. ; // deflection angle from the last hit

    if (  dynamic_cast<TVSurface *>( _det->At(idx)  )->CalcXingPointWith( *help0, xx, fid)  ){
      
      streamlog_out( DEBUG4 ) << " ---- crossing layer " << idx <<  " at: " 
			      << xx[0] << ", "  << xx[1] << ", " << xx[2] 
			      << " r: " << xx.Perp()  
			      << std::endl ;
      
      _xingPts[ idx ]  = new  gear::Vector3D( xx[0] , xx[1], xx[2] )  ;
      // _xingPts[ idx ]  = new  gear::Vector3D( xx[0]*10., xx[1]*10., xx[2]*10 )  ;

      ++idx ;
      
    } else {

      streamlog_out( DEBUG4 ) << " ---- no crossing point found with layer " << idx << std::endl ;      

      break ;
    }
  }


  if( streamlog_level( DEBUG4 ) ){

    TIter next(_kalHits, kIterBackward);

    TVTrackHit *hitp = 0;
    
    while ( (hitp = dynamic_cast<TVTrackHit *>( next() ) ) ) {
      
      int layer =  hitp->GetMeasLayer().GetIndex() ;
      
      if(  _xingPts[layer] != 0 ) {
	
	streamlog_out( ERROR )  << "   !!!!!!!!!!!!!!!!!!!! found xingPoint in layer  " << layer << " where we have a hit !!!!???????? " << std::endl ;
	
      }
    }
  }



}
