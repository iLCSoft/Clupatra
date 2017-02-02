#include "TrackEfficiency.h"


#include "marlin/AIDAProcessor.h"

//---- MarlinUtil 
#include "MarlinCED.h"

//#include <time.h>
#include <vector>
#include <map>
#include <algorithm>
#include <math.h>
#include <cmath>

//---- LCIO ---
#include "IMPL/LCCollectionVec.h"
#include "EVENT/Track.h"
#include "EVENT/SimTrackerHit.h"
#include "UTIL/Operators.h"
#include "UTIL/LCRelationNavigator.h"
#include "UTIL/LCTypedVector.h"
#include "UTIL/LCTypedVector.h"
#include "UTIL/ILDConf.h"
#include "UTIL/LCTOOLS.h"

//---- GEAR ----
// #include "marlin/Global.h"
#include "gear/GEAR.h"
// #include "gear/TPCParameters.h"
// #include "gear/PadRowLayout2D.h"
// #include "gear/BField.h"


//MarlinUtil
#include "HelixClass.h"


//---- ROOT -----
#include "TH1F.h" 

#include "UTIL/LCIterator.h"

using namespace lcio ;
using namespace marlin ;


void printSimTrackerHit(const EVENT::LCObject* o) ;
//----------------------------------------------------------------

TrackEfficiency aTrackEfficiency ;

// helper enum defining histogram index in vector 
namespace TrackEfficiencyHistos {
  enum index{
    hd0,   hphi,     homega, hz0, htanL, hpt,
    hed0,  hephi,    heomega, hez0, hetanL, hept,
    hdd0,  hdphi,    hdomega, hdz0, hdtanL, hdpt, 
    hpd0,  hpphi,    hpomega, hpz0, hptanL, hppt, 
    hd0mcp,hphimcp,  homegamcp, hz0mcp, htanLmcp, hptmcp,
    hdptp2,
    hpt_f,hcosth_f, hpt_t, hcosth_t, hacth_t, hacth_f,
    //-----  keep Size as last :
    Size   
  }; 


  class Histograms{
  public:
    Histograms(std::vector<TH1*>& v) : _h(&v) {}

    void create(int idx, const char* n, int nBin=100, double min=0., double max=0. ){
      create( idx , n , n , nBin, min , max ) ; 
    }

    void create(int idx, const char* n, const char* t,  int nBin=100, double min=0., double max=0. ){

      //    _h->resize( idx+1 ) ;
      _h->at( idx ) = new TH1D( n, t , nBin , min, max ) ;

      streamlog_out( DEBUG ) << " create histo " <<  n << " at index " << idx << std::endl ;
    }

    void create(int idx, const char* n, const char* t,  int nBin , double* bins ){

      _h->at( idx ) = new TH1D( n, t , nBin , bins ) ;

      streamlog_out( DEBUG ) << " create histo " <<  n << " at index " << idx << std::endl ;
    }

    void fill( int idx , double val, double weight=1.0 ){  _h->at( idx )->Fill( val , weight ) ; }

  protected:

    std::vector<TH1*>* _h;
  };



}
using namespace TrackEfficiencyHistos ;


//======================================================================================================
// 
#define APPLY_CUT( LEVEL, Cut, Exp )  if( (Exp) == false ) {  if ( Cut ) \
    streamlog_out( LEVEL ) << "  ***** failed cut:  [ " <<  #Exp \
			   <<  " ] in evt: " << evt->getEventNumber()	\
			   << " run: "  << evt->getRunNumber()   << std::endl ; \
    Cut = false ; }

// #define APPLY_CUT( LEVEL, Cut, Exp )  if( (Exp) == false ) { Cut = false ; \
//            streamlog_out( LEVEL ) << "  ***** failed cut:  [ " <<  #Exp \
// 			   <<  " ] in evt: " << evt->getEventNumber()	\
// 			   << " run: "  << evt->getRunNumber()   << std::endl ; }
//======================================================================================================

//======================================================================================================


TrackEfficiency::TrackEfficiency() : Processor("TrackEfficiency"),
				     _nRun(0), _nEvt(0) {
  
  // modify processor description
  _description = "TrackEfficiency: analysis plots for Tracks" ;
  
  
  
  registerInputCollection( LCIO::MCPARTICLE,
			   "MCParticlesCollection" , 
			   "Name of the input collections with MCParticles"  ,
			   _mcpColName ,
			   std::string("MCParticlesSkimmed" ) ) ;

  registerInputCollection( LCIO::SIMTRACKERHIT,
			   "SimTrackerHitCollection" , 
			   "Name of the input collections with SimTrackerHits"  ,
			   _sthColName ,
			   std::string("TPCCollection" ) ) ;

  registerInputCollection( LCIO::LCRELATION,
			   "TrackMCTruthRelation" , 
			   "Name of the input collections with LCRelation Tracks to MCParticles"  ,
			   _t2mColName ,
			   std::string("TrackMCTruthRelation" ) ) ;
  
  registerInputCollection( LCIO::LCRELATION,
			   "MCTruthTrackRelation" , 
			   "Name of the input collections with LCRelation MCParticles to Tracks"  ,
			   _m2tColName ,
			   std::string("MCTruthTrackRelation" ) ) ;
  

  registerInputCollection( LCIO::TRACK,
			   "TrackCollection" , 
			   "Name of the input collections with  Tracks"  ,
			   _trkColName ,
			   std::string("TrackCollection" ) ) ;



  FloatVec ptRange ;
  ptRange.push_back( 0. ) ;
  ptRange.push_back( 1000. ) ;

  registerProcessorParameter( "PtRange" , 
			      "Min and max value of pt range [GeV]"  ,
			      _ptRange ,
			      ptRange ) ;
  
}


void TrackEfficiency::init() { 

  // usually a good idea to
  printParameters() ;


  CEDPickingHandler::getInstance().registerFunction( LCIO::SIMTRACKERHIT , &printSimTrackerHit ) ; 

  if( _ptRange.size() != 2 ) {

    throw Exception("TrackEfficiency::init()  processor parameter 'ptRange' needs to have two values (min and max pt in GeV)" ) ;
  }


  _nRun = 0 ;
  _nEvt = 0 ;
}

void TrackEfficiency::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 

void TrackEfficiency::processEvent( LCEvent * evt ) { 


  clock_t start =  clock() ; 

  static const double bField = 3.5 ;   // should get this from Gear  ....

  static const double alpha =  2.99792458e-4 * bField ;  ;
  // p_t[GeV] = alpha / omega[mm] ;

  //static const double ptMin = _ptRange[0] ;
  //static const double ptMax = _ptRange[1] ;
  //static double ptMean = ( ptMax + ptMin ) / 2. ;


  //----------------------------------------------------------------------
  Histograms h(_h1) ;

  if( isFirstEvent() ){

    // this creates a directory for this processor ....
    AIDAProcessor::histogramFactory( this ) ;
    
    _h1.resize( TrackEfficiencyHistos::Size ) ;
    
    static const int nBin = 50 ;

    //  define binning for pt and omega
    // double ptBinL  = ptMin - ( ptMax - ptMin ) / 10.  ;
    // double ptBinH  = ptMax + ( ptMax - ptMin ) / 10.  ;
    double ptBinL  = 0.01 ; 
    double ptBinH  = 20. ;
    double omBinL  = alpha / ptBinH ;
    double omBinH  = alpha / ptBinL ;
    

    h.create( hd0,    "hd0"    ) ; 
    h.create( hphi,   "hphi"   ) ; 
    h.create( homega, "homega" ,"omega*mm" , nBin ,omBinL ,  omBinH ) ;
    h.create( hz0,    "hz0" , " z0 /mm  " , nBin , -3. , 3. ) ; 
    h.create( htanL,  "htanL"  ) ; 
    h.create( hpt,    "hpt"  , "pt/GeV" , nBin , ptBinL , ptBinH ) ; 

    h.create( hd0mcp,    "hd0mcp"    ) ; 
    h.create( hphimcp,   "hphimcp"   ) ; 
    h.create( homegamcp, "homegamcp" ,"omega*mm" , nBin ,omBinL ,  omBinH ) ;  //, "omega*mm" , 50, 0.0003 , 0.0004 ) ; 
    h.create( hz0mcp,    "hz0mcp" ) ; //, " z0 /mm  " , 50 , -3. , 3. ) ; 
    h.create( htanLmcp,  "htanLmcp"  ) ; 
    h.create( hptmcp,    "hptmcp" , "pt/GeV" , nBin , ptBinL , ptBinH ) ;
    
    h.create( hed0,    "hed0"    ) ; 
    h.create( hephi,   "hephi"   ) ; 
    h.create( heomega, "heomega" ) ;
    h.create( hez0,    "hez0"    ) ;
    h.create( hetanL,  "hetanL"  ) ; 
    h.create( hept,    "hept"  ) ; 

    h.create( hpd0,    "hpd0"    , "d0 pull   ", 50 , -5. , 5. ) ; 
    h.create( hpphi,   "hpphi"   , "phi pull  ", 50 , -5. , 5. ) ; 
    h.create( hpomega, "hpomega" , "omega pull", 50 , -5. , 5. ) ;
    h.create( hpz0,    "hpz0"    , "z0 pull   ", 50 , -5. , 5. ) ;
    h.create( hptanL,  "hptanL"  , "tanL pull ", 50 , -5. , 5. ) ; 
    h.create( hppt,    "hppt"    , "pt pull   ", 50 , -5. , 5. ) ; 

    h.create( hdd0,    "hdd0"    ) ; 
    h.create( hdphi,   "hdphi"   ) ; 
    h.create( hdomega, "hdomega" ) ;
    h.create( hdz0,    "hdz0"    ) ;
    h.create( hdtanL,  "hdtanL"  ) ; 
    h.create( hdpt,    "hdpt"  ) ; 

    h.create( hdptp2,    "hdptp2"  ) ; 

    h.create( hcosth_t,    "hcosth_t"    , "cos theta - true    ", 21 , -1. , 1. ) ; 
    h.create( hcosth_f,    "hcosth_f"    , "cos theta - found   ", 21 , -1. , 1. ) ; 

    h.create( hacth_t,    "hacth_t"    , "|cos theta| - true    ", 20 , 0. , 1. ) ; 
    h.create( hacth_f,    "hacth_f"    , "|cos theta| - found   ", 20 , 0. , 1. ) ; 

    //  void create(int idx, const char* n, const char* t,  int nBin=100, double min=0., double max=0. ){

    // variable  bin size histogram for pt
    static const int nBins = 13 ;
    double bins[nBins+1] = { 0.1, 0.2, 0.4, 0.6 , 0.8 , 1.0 , 2., 5.0 , 10. , 20. , 50. , 100. , 300. , 500. } ;
//     static const int nBins = 9 ;
//     double bins[nBins+1] = { 0.1, 0.2, 0.5 , 1.0 , 2., 5.0 , 10. , 20. , 50. , 100.} ;

    h.create( hpt_t, "hpt_t", "pt - true ", nBins , bins ) ;
    h.create( hpt_f, "hpt_f", "pt - found ", nBins , bins ) ;
    
  }
  //=========================================================================================================================


  LCCollection* m2tCol = 0 ;
  try{
    m2tCol = evt->getCollection( _m2tColName ) ;
  }
  catch( lcio::DataNotAvailableException& e){
    streamlog_out( DEBUG ) << " *** collection not in event : " << _m2tColName << std::endl ;
    return ; // nothing to do in this event (no tracks)
  }
  LCCollection* t2mCol = 0 ;
  try{
    t2mCol = evt->getCollection( _t2mColName ) ;
  }
  catch( lcio::DataNotAvailableException& e){
    streamlog_out( DEBUG ) << " *** collection not in event : " << _t2mColName << std::endl ;
    return ; // nothing to do in this event (no tracks)
  }


  LCCollection* trkCol = 0 ;
  try{
    
    trkCol = evt->getCollection( _trkColName ) ;
  }
  catch( lcio::DataNotAvailableException& e){

    streamlog_out( DEBUG ) << " *** collection not ine event : " << _trkColName << std::endl ;
    return ; // nothing to do in this event (no tracks)
    
  }
  if( trkCol->getNumberOfElements() < 1 ) {

    streamlog_out( ERROR ) << " not tracks in collection " << _trkColName << " evt : " << evt->getEventNumber()
			   << " run : " << evt->getRunNumber() << std::endl ;
    return ;
  }


  LCRelationNavigator navm2t( m2tCol  ) ;
  LCRelationNavigator navt2m( t2mCol  ) ;


  LCIterator<MCParticle> mcpIt( evt, _mcpColName ) ;
  
  LCCollectionVec* mcpTracks = new LCCollectionVec( LCIO::MCPARTICLE )  ;
  mcpTracks->setSubset( true ) ;
  mcpTracks->reserve(  mcpIt.size() ) ;
  
  std::string name = "MCParticleTracks_" ;
  name += this->name() ;
  
  evt->addCollection( mcpTracks , name  ) ;


  streamlog_out( DEBUG6 ) << " ----- adding collection " << name << " to the event " << std::endl ;

  LCCollectionVec* mcpTrksNotFound = new LCCollectionVec( LCIO::MCPARTICLE )  ;
  mcpTrksNotFound->setSubset( true ) ;
  mcpTrksNotFound->reserve(  mcpIt.size() ) ;

  name = "MCParticleTracksNotFound_" ;
  name += this->name() ;
  evt->addCollection( mcpTrksNotFound ,  name  ) ;

  streamlog_out( DEBUG6 ) << " ----- adding collection " << name << " to the event " << std::endl ;


  LCCollectionVec* splitTracks = new LCCollectionVec( LCIO::TRACK )  ;
  splitTracks->setSubset( true ) ;
  splitTracks->reserve(  trkCol->getNumberOfElements() ) ;
  name = "SplitTracks_" ;
  name += this->name() ;
  evt->addCollection( splitTracks ,  name  ) ;
  
  //------ count sim hits from every MCParticle
  typedef std::map< MCParticle* , int > HITMAP ;
  HITMAP hitMap ;
  try{
    LCIterator<SimTrackerHit> sthIt( evt, _sthColName ) ;
    
    while( SimTrackerHit* sth = sthIt.next()  ){
      
      hitMap[ sth->getMCParticle() ]++ ;
    }
  }catch( lcio::DataNotAvailableException& e){

  }

  //--- check if we hav any tracks in the relation (not always for LEP trk)


  while( MCParticle* mcp = mcpIt.next()  ){
    
 
    bool cut = true ;
    //    APPLY_CUT( DEBUG, cut,  std::abs( mcp->getCharge() )  > 0.5  ) ;
    if(  std::abs( mcp->getCharge() ) < 0.5 ) cut = false ; // silent cut
    if( cut )
      streamlog_out( DEBUG ) <<  lcshort( mcp ) << std::endl ;

    //    APPLY_CUT( DEBUG, cut,  mcp->getGeneratorStatus() == 1   ) ;   // no documentation lines

    gear::Vector3D v( mcp->getVertex()[0], mcp->getVertex()[1], mcp->getVertex()[2] );
    gear::Vector3D e( mcp->getEndpoint()[0], mcp->getEndpoint()[1], mcp->getEndpoint()[2] );
    gear::Vector3D p( mcp->getMomentum()[0], mcp->getMomentum()[1], mcp->getMomentum()[2] );

    APPLY_CUT( DEBUG, cut, v.r() < 100.   ) ;   // start at IP+/-10cm

    // ==== vzeros =========
    // APPLY_CUT( DEBUG, cut, v.rho() > 100.   ) ;   // non prompt track
    // bool isVzero = false ;
    // if(   mcp->getParents().size() > 0  ) {
    //   isVzero =  ( std::abs( mcp->getParents()[0]->getCharge() )  < 0.1   ) ;
    // }
    // APPLY_CUT( DEBUG, cut,  isVzero    ) ;   // start at IP+/-10cm
    // ==== vzeros =========

    APPLY_CUT( DEBUG, cut, e.rho()==0.0  || e.rho() > 40.   ) ; // end at rho > 40 cm

    APPLY_CUT( DEBUG, cut, p.rho() > .1 ) ; //FIXME 1. Gev <->  pt> 100 MeV

    APPLY_CUT( DEBUG, cut, std::abs( cos( p.theta() ) )  < 0.99  ) ; //FIXME 0.9 <=> .99  //  | cos( theta ) | > 0.99

    APPLY_CUT( DEBUG, cut, hitMap[ mcp ]  > 5 ) ; //  require at least 5 tracker hits 

    //....

    if( cut ) 
      mcpTracks->push_back( mcp ) ;
  }
  
 

  name = "MCParticleTracks_" ;
  name += this->name() ;

  LCIterator<MCParticle> trmIt( evt, name ) ;

  //  while( MCParticle* trm = trmIt.next()  ){
  for(int ii=0,nn=mcpTracks->getNumberOfElements() ; ii<nn ; ++ii) {

    MCParticle* trm = (MCParticle*) mcpTracks->getElementAt(ii) ; 
    
    //    std::cout << "    ----- searching particle track : " << lcshort( trm ) << std::endl;

    double pxmcp =  trm->getMomentum()[0]  ;
    double pymcp =  trm->getMomentum()[1]  ;
    //    double pzmcp =  trm->getMomentum()[2]  ;
    
    //    double pmcp  = sqrt( pxmcp*pxmcp + pymcp*pymcp + pzmcp*pzmcp ) ;
    double ptmcp = sqrt( pxmcp*pxmcp + pymcp*pymcp ) ;
    
    // double thmcp  = atan2( ptmcp , pzmcp ) ;

    // double d0mcp = 0. ;  //fixme: only true for prompt tracks
    // double phmcp = atan2( pymcp, pxmcp) ;
    // double ommcp = alpha / ptmcp ; 
    // if(  trm->getCharge() < 0. ) 
    //   ommcp = -ommcp ;
    // double z0mcp = 0 ;  //fixme: only true for prompt tracks
    //    double tLmcp = pzmcp / ptmcp  ; 


    HelixClass helix ;

    float pos[3] ;
    pos[0] = trm->getVertex()[0] ;
    pos[1] = trm->getVertex()[1] ;
    pos[2] = trm->getVertex()[2] ;
    float mom[3] ;
    mom[0] = trm->getMomentum()[0] ;
    mom[1] = trm->getMomentum()[1] ;
    mom[2] = trm->getMomentum()[2] ;

    gear::Vector3D p( trm->getMomentum()[0], trm->getMomentum()[1], trm->getMomentum()[2] );
    double costhmcp  = cos( p.theta() ) ;

    float q = trm->getCharge() ;
    
    helix.Initialize_VP( pos , mom, q,  3.5 ) ;

    double d0mcp = helix.getD0() ;
    double phmcp = helix.getPhi0() ;
    double ommcp = helix.getOmega() ; 
    double z0mcp = helix.getZ0() ;
    double tLmcp = helix.getTanLambda() ; 


    h.fill( hd0mcp,    d0mcp ) ;
    h.fill( hphimcp,   phmcp ) ;
    h.fill( homegamcp, ommcp ) ;
    h.fill( hz0mcp,    z0mcp ) ;
    h.fill( htanLmcp,  tLmcp ) ;
    h.fill( hptmcp,    ptmcp ) ;
    

    h.fill( hpt_t , ptmcp ) ;
    h.fill( hcosth_t , costhmcp ) ;
    h.fill( hacth_t , std::abs( costhmcp) ) ;
    

  
    bool useT2M = true ;

    const EVENT::LCObjectVec& trkV = ( useT2M ?   navt2m.getRelatedFromObjects( trm ) :   navm2t.getRelatedToObjects( trm )  ); 
    
    //    std::cout <<  " ------   trkV.size()  : " << trkV.size()  << std::endl ;

    bool foundTrack = false ;

    if( trkV.size() >  0 ){
      

      const FloatVec& wV =  ( useT2M ?  navt2m.getRelatedFromWeights( trm ) :  navm2t.getRelatedToWeights( trm ) ) ;

      double wMax = 0.0 ;
      int iMax = 0 ;

      streamlog_out( DEBUG2 ) <<  " ------   trkV.size()  : " << trkV.size()  << std::endl ;

      for(unsigned i=0 ; i<wV.size() ; ++i ) {
	 
	streamlog_out( DEBUG2 ) <<  "          nhit  : " << ((Track*)trkV[i])->getTrackerHits().size() <<    " ---- weight : " <<   wV[i]   << std::endl ;

	if( wV[i] > wMax ){
	  wMax =  wV[i] ;
	  iMax = i ;
	}  
      } 

      double minGoodHitFraction = 0.75 ; // 90 ;

      int nSplitSegments = 0 ;
      
      // store split tracks collection
      if( trkV.size() >  1 ){
	
	for(unsigned i=0,N=trkV.size() ; i<N ; ++i ) {
	  
	  //	  if(  wV[i] > minGoodHitFraction ) {
	    
	    ++nSplitSegments ;
	    
	    splitTracks->addElement( trkV[i] ) ;
	    //}
	}
      }


      Track* tr = (Track*) ( trkV.at(iMax)  ) ; 
 
      double d0 = tr->getD0() ;
      double ph = tr->getPhi() ;
      double om = tr->getOmega() ;
      double z0 = tr->getZ0() ;
      double tL = tr->getTanLambda() ;
      double pt = std::abs( alpha / om ) ;

      double ed0 = sqrt( tr->getCovMatrix()[0]  );
      double eph = sqrt( tr->getCovMatrix()[2]  );
      double eom = sqrt( tr->getCovMatrix()[5]  );
      double ez0 = sqrt( tr->getCovMatrix()[9]  );
      double etL = sqrt( tr->getCovMatrix()[14] );
      
      double ept = alpha * eom / ( om * om )  ; // ept = al/om^2 * eom 
    
      
      double dd0 = d0 - d0mcp ;
      double dph = ph - phmcp ; 
      if( std::abs(dph) > M_PI )
	dph = 2.* M_PI - std::abs( dph ) ;


      double dom = om - ommcp ; 
      double dz0 = z0 - z0mcp ; 
      double dtL = tL - tLmcp ; 
      double dpt = pt - ptmcp ; 

      bool cut = true ;
      
      //       APPLY_CUT( DEBUG4, cut,  std::abs( dph )  <  (3.*eph)  ) ;
      //       APPLY_CUT( DEBUG4, cut,  std::abs( dom )  <  (5.*eom)  ) ;              // require 5 sigma on omega only for now
      //       APPLY_CUT( DEBUG4, cut,  std::abs( dtL  ) <  (3.*etL)  ) ;


      // trk->subdetectorHitNumbers()[ 2 * ILDDetID::TPC - 1 ] =  hitsInFit ;  
      // trk->subdetectorHitNumbers()[ 2 * ILDDetID::TPC - 2 ] =  hitCount ;  
      int nTPCHit = tr->getSubdetectorHitNumbers()[ 2 * ILDDetID::TPC - 2 ] ;
      int nMCPTPCHit = hitMap[ trm ] ;
      //      double goodHitFraction = ( wMax * nTPCHit )  / nMCPTPCHit ;

      double goodHitFraction = wMax ;
      
      APPLY_CUT( DEBUG4, cut,  goodHitFraction > minGoodHitFraction  ) ;
      //	APPLY_CUT( DEBUG4, cut,  nSplitSegments < 2 ) ;  // don't count split tracks as found....
      
      //      std::cout   << " ========= track  found ?  -> " << cut << lcshort( trm ) << std::endl ;

     if( ! cut ) {
	streamlog_out(DEBUG3) << "  goodHitFraction : " << wMax<< "*" << nTPCHit << "/" << nMCPTPCHit 
			      << "=" << goodHitFraction << std::endl ;
	
	streamlog_out(DEBUG3) << " phi : " << dph << "  -  " <<  3.*eph << " - " << dph /ph << std::endl ;
	streamlog_out(DEBUG3) << " ome : " << dom << "  -  " <<  3.*eom << " - " << dom /om << std::endl ;
	streamlog_out(DEBUG3) << " tanL: " << dtL << "  -  " <<  3.*etL << " - " << dtL /tL << std::endl ;
	streamlog_out(DEBUG3) << " pt: " << pt <<  " dpt: " << dpt << " ept : " << ept << " theta " <<  180. * atan( 1. / tL ) / M_PI  << std::endl ;
	

      }
      

      if( cut == true ){
	
	foundTrack = true ;

	h.fill( hd0,    d0 ) ;
	h.fill( hphi,   ph ) ;
	h.fill( homega, om ) ;
	h.fill( hz0,    z0 ) ;
	h.fill( htanL,  tL ) ;
	h.fill( hpt,    pt ) ;
      
	h.fill( hed0,    ed0 ) ;
	h.fill( hephi,   eph ) ;
	h.fill( heomega, eom ) ;
	h.fill( hez0,    ez0 ) ;
	h.fill( hetanL,  etL ) ;
	h.fill( hept,    ept ) ; 


   
	h.fill( hdd0,    dd0 ) ;
	h.fill( hdphi,   dph ) ;
	h.fill( hdomega, dom ) ;
	h.fill( hdz0,    dz0 ) ;
	h.fill( hdtanL,  dtL ) ;
	h.fill( hdpt,    dpt ) ;
	
	h.fill( hpd0,    dd0 / ed0 ) ;
	h.fill( hpphi,   dph / eph ) ;
	h.fill( hpomega, dom / eom ) ;
	h.fill( hpz0,    dz0 / ez0 ) ;
	h.fill( hptanL,  dtL / etL ) ;
	h.fill( hppt,    dpt / ept ) ;

	h.fill( hdptp2,    dpt / (pt*pt) ) ;


	h.fill( hpt_f , ptmcp ) ;
	h.fill( hcosth_f , costhmcp ) ;
	h.fill( hacth_f , std::abs( costhmcp) ) ;

	
      }
    }
    
    if( ! foundTrack )  {
      mcpTrksNotFound->push_back( trm ) ;

      streamlog_out( DEBUG3 )  << " ========= track not found " << lcshort( trm ) << std::endl ;
    }
  }

  streamlog_out( DEBUG4 )  << " ===== could not find " << mcpTrksNotFound->size()  <<  " tracks of " << mcpTracks->size()  
			   << std::endl ;

  //----------------------------------------------------------------------  
  
  _nEvt ++ ;

  clock_t end = clock () ; 
  
  streamlog_out( DEBUG )  << "---  histogramming time: " 
  			  <<  double( end - start ) / double(CLOCKS_PER_SEC) << std::endl  ;

}




/*************************************************************************************************/
void TrackEfficiency::check( LCEvent * evt ) { 
  /*************************************************************************************************/

}
//====================================================================================




void TrackEfficiency::end(){ 
  
  streamlog_out( MESSAGE )  << " processed " << _nEvt << " events in " << _nRun << " runs "
			    << std::endl ;
}

//====================================================================================================

