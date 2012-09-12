#include "TrackCheckMCTruth.h"


#include "marlin/AIDAProcessor.h"

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
#include "UTIL/BitSet32.h"

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


TrackCheckMCTruth aTrackCheckMCTruth ;

/** helper method to create a track collections and add it to the event */
inline IMPL::LCCollectionVec* newTrkCol(const std::string& name, LCEvent * evt ){

  IMPL::LCCollectionVec* col = new IMPL::LCCollectionVec( LCIO::TRACK ) ;  

  col->setFlag( UTIL::make_bitset32(  LCIO::TRBIT_HITS )  ) ;


  evt->addCollection( col , name ) ;

  return col ;
}
//----------------------------------------------------------------
/** helper method to get the collection from the event */
inline LCCollection* getCollection(  const std::string& name, LCEvent * evt ){

  LCCollection* col = 0 ;
  try{   col = evt->getCollection( name )  ; 
  } catch( lcio::DataNotAvailableException& e) { 
    streamlog_out( DEBUG4 ) <<  " input collection not in event : " << name << "  !!!  " << std::endl ;  
  } 
  return col ;
}

//----------------------------------------------------------------

// helper enum defining histogram index in vector 
namespace TrackCheckMCTruthHistos {
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
using namespace TrackCheckMCTruthHistos ;


//======================================================================================================

#define APPLY_CUT( LEVEL, Cut, Exp )  if( (Exp) == false ) { Cut = false ; \
    streamlog_out( LEVEL )<< "  ***** failed cut:  [ " <<  #Exp		\
			  <<  " ] in evt: " << evt->getEventNumber()	\
			  << " run: "  << evt->getRunNumber()   << std::endl ; }

//======================================================================================================
//    streamlog_out( LEVEL )
//======================================================================================================


TrackCheckMCTruth::TrackCheckMCTruth() : Processor("TrackCheckMCTruth"),
					 _nRun(0), _nEvt(0)  {
  
  // modify processor description
  _description = "TrackCheckMCTruth: analysis plots for Tracks" ;
  
  
  
  registerInputCollection( LCIO::MCPARTICLE,
			   "MCParticlesCollection" , 
			   "Name of the input collections with MCParticles"  ,
			   _mcpColName ,
			   std::string("MCParticlesSkimmed" ) ) ;

  StringVec sthVec ;
  sthVec.push_back( "TPCCollection" );

  registerInputCollections( LCIO::SIMTRACKERHIT,
			    "SimTrackerHitCollection" , 
			    "Name of the input collections with SimTrackerHits"  ,
			    _sthColNames ,
			    sthVec );
  
  registerInputCollection( LCIO::LCRELATION,
			   "TrackMCPRelation" , 
			   "Name of the input collections with LCRelation Tracks to MCParticles"  ,
			   _relTrkMCPColName ,
			   std::string("TrackRelation" ) ) ;
  
  registerInputCollection( LCIO::LCRELATION,
			   "MCPTrackRelation" , 
			   "Name of the input collections with LCRelation MCParticles to tracks"  ,
			   _relMCPTrkColName ,
			   std::string("MCPTrkRelation" ) ) ;
  
  registerInputCollection( LCIO::TRACK,
			   "TrackCollection" , 
			   "Name of the input collections with  Tracks"  ,
			   _trkColName ,
			   std::string("ClupatraTracks" ) ) ;


  registerOutputCollection( LCIO::TRACK,
			    "SplitTrackCollection" , 
			    "Name of the output collections with split Tracks"  ,
			    _splitColName ,
			    std::string("ClupatraSplitTracks" ) ) ;

  registerOutputCollection( LCIO::TRACK,
			    "MergedTrackCollection" , 
			    "Name of the output collections with merged Tracks"  ,
			    _mergedColName ,
			    std::string("ClupatraMergedTracks" ) ) ;

  registerProcessorParameter( "HitFractionSplit" , 
			      "maximal fractions of hits that are allowed to be missing for a split track" ,
			      _hitFractionSplit ,
			      float( 0.05 )  );

  registerProcessorParameter( "HitFractionMerged" , 
			      "maximal fractions of hits that are allowed to be missing for a merged track" ,
			      _hitFractionMerged ,
			      float( 0.05 )  );

  FloatVec ptRange ;
  ptRange.push_back( 0. ) ;
  ptRange.push_back( 1000. ) ;

  registerProcessorParameter( "PtRange" , 
			      "Min and max value of pt range [GeV]"  ,
			      _ptRange ,
			      ptRange ) ;
  
}


void TrackCheckMCTruth::init() { 

  // usually a good idea to
  printParameters() ;


  if( _ptRange.size() != 2 ) {

    throw Exception("TrackCheckMCTruth::init()  processor parameter 'ptRange' needs to have two values (min and max pt in GeV)" ) ;
  }


  _nRun = 0 ;
  _nEvt = 0 ;
}

void TrackCheckMCTruth::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 

void TrackCheckMCTruth::processEvent( LCEvent * evt ) { 


  clock_t start =  clock() ; 

  static const double bField = 3.5 ;   // should get this from Gear  ....

  static const double alpha =  2.99792458e-4 * bField ;  ;
  // p_t[GeV] = alpha / omega[mm] ;


  //----------------------------------------------------------------------
  Histograms h(_h1) ;

  if( isFirstEvent() ){

    // this creates a directory for this processor ....
    AIDAProcessor::histogramFactory( this ) ;
    
    _h1.resize( TrackCheckMCTruthHistos::Size ) ;
    
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

    h.create( hacth_t,    "hacth_t"    , "|cos theta| - true    ", 10 , 0. , 1. ) ; 
    h.create( hacth_f,    "hacth_f"    , "|cos theta| - found   ", 10 , 0. , 1. ) ; 

    //  void create(int idx, const char* n, const char* t,  int nBin=100, double min=0., double max=0. ){

    // variable  bin size histogram for pt
    static const int nBins = 9 ;
    double bins[nBins+1] = { 0.1, 0.2, 0.5 , 1.0 , 2., 5.0 , 10. , 20. , 50. , 100. } ;

    h.create( hpt_t, "hpt_t", "pt - true ", nBins , bins ) ;
    h.create( hpt_f, "hpt_f", "pt - found ", nBins , bins ) ;
    
  }
  //=========================================================================================================================


  LCCollection* relCol = getCollection( _relTrkMCPColName, evt ) ;
  if( ! relCol ) return ;
  
  LCCollection* relInverseCol = getCollection( _relMCPTrkColName, evt ) ;
  if( ! relInverseCol ) return ;
  
  LCCollection* trkCol = getCollection( _trkColName , evt ) ;
  if( ! trkCol ) return ;


  if( trkCol->getNumberOfElements() < 1 ) {
    
    streamlog_out( ERROR ) << " not tracks in collection " << _trkColName << " evt : " << evt->getEventNumber()
			   << " run : " << evt->getRunNumber() << std::endl ;
    return ;
  }

  LCCollectionVec* splitCol = newTrkCol( _splitColName, evt ) ;
  LCCollectionVec* mergedCol = newTrkCol( _mergedColName, evt ) ;
  splitCol->setSubset() ;
  mergedCol->setSubset() ;



  LCRelationNavigator nav( relCol  ) ;
  
  LCRelationNavigator navInverse( relInverseCol  ) ;

  //---------------------------------------------------------------------------------
  // create collections with low, ok and high weight of the truth relation
  //

  LCCollectionVec* laWCol = newTrkCol( "ClupatraTracksLargeWeight" , evt ) ;
  LCCollectionVec* okWCol = newTrkCol( "ClupatraTracksOkWeight"    , evt ) ;
  LCCollectionVec* loWCol = newTrkCol( "ClupatraTracksLowWeight"   , evt ) ;
  
  laWCol->setSubset() ;
  okWCol->setSubset() ;
  loWCol->setSubset() ;
  
  for( LCIterator<LCRelation> it( relCol) ; LCRelation* rel = it.next() ; ){ 
    
    float weight = rel->getWeight() ;
    Track* trk   = (Track*)rel->getFrom() ;
    
    
    if( weight <= 0.9 )
      loWCol->addElement( trk ) ;
    else if( weight > 0.9  && weight <=  1.0 )
      okWCol->addElement( trk ) ;
    else if( weight > 1. )
      laWCol->addElement( trk ) ;
  }

  //---------------------------------------------------------------------------------
  // create collections with low, ok and high weight of the inverse truth relation
  //

  LCCollectionVec* laWMCPCol = newTrkCol( "ClupatraTracksLargeMCPWeight" , evt ) ;
  LCCollectionVec* okWMCPCol = newTrkCol( "ClupatraTracksOkMCPWeight"    , evt ) ;
  LCCollectionVec* loWMCPCol = newTrkCol( "ClupatraTracksLowMCPWeight"   , evt ) ;
  
  laWMCPCol->setSubset() ;
  okWMCPCol->setSubset() ;
  loWMCPCol->setSubset() ;
  
  for( LCIterator<LCRelation> it( relInverseCol) ; LCRelation* rel = it.next() ; ){ 
    
    float weight = rel->getWeight() ;
    Track* trk   = (Track*)rel->getTo() ;
    
    
    if( weight <= 0.7 )
      loWMCPCol->addElement( trk ) ;
    else if( weight > 0.7  && weight <=  1.0 )
      okWMCPCol->addElement( trk ) ;
    else if( weight > 1. )
      laWMCPCol->addElement( trk ) ;
  }
  //---------------------------------------------------------------------------------


  LCIterator<MCParticle> mcpIt( evt, _mcpColName ) ;
  
  LCCollectionVec* mcpTracks = new LCCollectionVec( LCIO::MCPARTICLE )  ;
  mcpTracks->setSubset( true ) ;
  mcpTracks->reserve(  mcpIt.size() ) ;
  
  std::string name = "MCParticleTracks_" ;
  name += this->name() ;
  
  evt->addCollection( mcpTracks , name  ) ;

  LCCollectionVec* mcpTrksFound = new LCCollectionVec( LCIO::MCPARTICLE )  ;
  mcpTrksFound->setSubset( true ) ;
  mcpTrksFound->reserve(  mcpIt.size() ) ;

  name = "MCParticleTracksFound_" ;
  name += this->name() ;
  evt->addCollection( mcpTrksFound ,  name  ) ;


  //------ count sim hits from every MCParticle

  typedef std::map< MCParticle* , int > HITMAP ;
  HITMAP hitMap ;

  for( unsigned i=0 ,N = _sthColNames.size() ; i<N ; ++i ){
   
    try {
      LCIterator<SimTrackerHit> sthIt( evt, _sthColNames[i] ) ;
      
      while( SimTrackerHit* sth = sthIt.next()  ){
	
	hitMap[ sth->getMCParticle() ]++ ;
      }
    } catch( lcio::DataNotAvailableException& e){ 

      streamlog_out( DEBUG4 ) << " missing input collection: " << _sthColNames[i] << std::endl ;
    }
  }



  //================ create denominator for efficency ========================

  while( MCParticle* mcp = mcpIt.next()  ){
    
    bool cut = true ;
    APPLY_CUT( DEBUG2, cut,  std::abs( mcp->getCharge() )  > 0.5  ) ;
    
    APPLY_CUT( DEBUG2, cut,  mcp->getGeneratorStatus() != 2   ) ;   // no documentation lines

    gear::Vector3D v( mcp->getVertex()[0], mcp->getVertex()[1], mcp->getVertex()[2] );
    gear::Vector3D e( mcp->getEndpoint()[0], mcp->getEndpoint()[1], mcp->getEndpoint()[2] );
    gear::Vector3D p( mcp->getMomentum()[0], mcp->getMomentum()[1], mcp->getMomentum()[2] );

    APPLY_CUT( DEBUG2, cut, v.r() < 100.   ) ;   // start at IP+/-10cm

    // ==== vzeros =========
    // APPLY_CUT( DEBUG, cut, v.rho() > 100.   ) ;   // non prompt track
    // bool isVzero = false ;
    // if(   mcp->getParents().size() > 0  ) {
    //   isVzero =  ( std::abs( mcp->getParents()[0]->getCharge() )  < 0.1   ) ;
    // }
    // APPLY_CUT( DEBUG, cut,  isVzero    ) ;   // start at IP+/-10cm
    // ==== vzeros =========

    APPLY_CUT( DEBUG2, cut, e.rho()==0.0  || e.rho() > 400.   ) ; // end at rho > 40 cm

    APPLY_CUT( DEBUG2, cut, p.rho() > 0.1  ) ; // pt> 100 MeV

    APPLY_CUT( DEBUG2, cut, std::abs( cos( p.theta() ) )  < 0.99  ) ; //  | cos( theta ) | > 0.99

    APPLY_CUT( DEBUG2, cut, hitMap[ mcp ]  > 5 ) ; //  require at least 10 TPC hits (particle actually made it to the TPC)

    //....

    if( cut ) 
      mcpTracks->push_back( mcp ) ;
  }
  
  //================================================================

  // name = "MCParticleTracks_" ;
  // name += this->name() ;
  // LCIterator<MCParticle> trmIt( evt, name ) ;

  LCIterator<MCParticle> trmIt( mcpTracks );

  while( MCParticle* trm = trmIt.next()  ){
    
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
    
    const EVENT::LCObjectVec& trkV = nav.getRelatedFromObjects( trm ) ;

    std::vector<int> splitTrackIndices ;
    splitTrackIndices.reserve( trkV.size()  ) ;

    if( trkV.size() >  0 ){
      
      const FloatVec& wV = nav.getRelatedFromWeights( trm ) ;
      double wMax = 0.0 ;
      int iMax = 0 ;
      
      for(unsigned i=0 ; i<wV.size() ; ++i ) {
	
	if( wV[i] > wMax ){
	  wMax =  wV[i] ;
	  iMax = i ;
	}  

	if( wV[i] > 0.5 ) { // if we have a track segment that has a weight >0.5 it could be part of a split track....
	  splitTrackIndices.push_back( i ) ;
	}
      } 
      
      // if(  trkV.size() > 1 ){ // if we have more than one track it is a split track....
      // 	std::copy( trkV.begin() , trkV.end() , std::back_inserter( *splitCol ) ) ;
      // }
      if( splitTrackIndices.size() > 1 ){

	for(unsigned i=0,N=splitTrackIndices.size() ; i<N ; ++i ) {
	  splitCol->addElement(  trkV[i ] ) ;
	}
      }


      Track* tr = (Track*)( trkV.at(iMax)  ) ; 

 
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
      
      //      APPLY_CUT( DEBUG4, cut,  std::abs( dph )  <  (3.*eph)  ) ;
      //      APPLY_CUT( DEBUG4, cut,  std::abs( dom )  <  (5.*eom)  ) ;              // require 5 sigma on omega only for now
      //      APPLY_CUT( DEBUG4, cut,  std::abs( dtL  ) <  (3.*etL)  ) ;

      if( ! cut ) {
	streamlog_out(DEBUG3) << " phi : " << dph << "  -  " <<  3.*eph << std::endl ;
	streamlog_out(DEBUG3) << " ome : " << dom << "  -  " <<  3.*eom << std::endl ;
	streamlog_out(DEBUG3) << " tanL: " << dtL << "  -  " <<  3.*etL << std::endl ;
	streamlog_out(DEBUG3) << " pt: " << pt <<  " dpt: " << dpt << " ept : " << ept << " theta " <<  180. * atan( 1. / tL ) / M_PI  << std::endl ;
      }

      if( cut == true ){
	
	mcpTrksFound->push_back( trm ) ;
	
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
  }

  streamlog_out( DEBUG4 )  << " ===== found " << mcpTrksFound->size()  <<  " tracks of " << mcpTracks->size()  
			   << std::endl ;

  //----------------------------------------------------------------------  
  
  _nEvt ++ ;

  clock_t end = clock () ; 
  
  streamlog_out( DEBUG )  << "---  histogramming time: " 
  			  <<  double( end - start ) / double(CLOCKS_PER_SEC) << std::endl  ;

}




/*************************************************************************************************/
void TrackCheckMCTruth::check( LCEvent * evt ) { 
  /*************************************************************************************************/

  }
//====================================================================================




void TrackCheckMCTruth::end(){ 
  
  streamlog_out( MESSAGE )  << " processed " << _nEvt << " events in " << _nRun << " runs "
			    << std::endl ;
}

//====================================================================================================

