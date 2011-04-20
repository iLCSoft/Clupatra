#include "TrackEfficiency.h"


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
#include "UTIL/Operators.h"
#include "UTIL/LCRelationNavigator.h"
#include "UTIL/LCTypedVector.h"

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

#include "LCIterator.h"

using namespace lcio ;
using namespace marlin ;


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

  void fill( int idx , double val, double weight=1.0 ){  _h->at( idx )->Fill( val , weight ) ; }

protected:

  std::vector<TH1*>* _h;
};



}
using namespace TrackEfficiencyHistos ;


//======================================================================================================

#define APPLY_CUT( LEVEL, Cut, Exp )  if( (Exp) == false ) { Cut = false ; \
  streamlog_out( LEVEL ) << "  ***** failed cut:  [ " <<  #Exp	    \
  <<  " ] in evt: " << evt->getEventNumber()			    \
  << " run: "  << evt->getRunNumber()   << std::endl ; }

//======================================================================================================

//======================================================================================================


TrackEfficiency::TrackEfficiency() : Processor("TrackEfficiency") {
  
  // modify processor description
  _description = "TrackEfficiency: analysis plots for Tracks" ;
  
  
  
  registerInputCollection( LCIO::MCPARTICLE,
			   "MCParticlesCollection" , 
			   "Name of the input collections with MCParticles"  ,
			   _mcpColName ,
			   std::string("MCParticlesSkimmed" ) ) ;

  registerInputCollection( LCIO::LCRELATION,
			   "TrackMCPRelation" , 
			   "Name of the input collections with LCRelation Tracks to MCParticles"  ,
			   _relColName ,
			   std::string("TrackRelation" ) ) ;
  
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

  static const double ptMin = _ptRange[0] ;
  static const double ptMax = _ptRange[1] ;

  static double ptMean = ( ptMax + ptMin ) / 2. ;


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
    double ptBinL  = 0. ; 
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

  }
  //=========================================================================================================================


  LCRelationNavigator nav( evt->getCollection( _relColName )  ) ;
  

  LCIterator<MCParticle> mcpIt( evt, _mcpColName ) ;

  LCCollectionVec* mcpTracks = new LCCollectionVec( LCIO::MCPARTICLE )  ;
  mcpTracks->setSubset( true ) ;
  mcpTracks->reserve(  mcpIt.size() ) ;
  evt->addCollection( mcpTracks , "MCParticleTracks"   ) ;

  LCCollectionVec* mcpTrksFound = new LCCollectionVec( LCIO::MCPARTICLE )  ;
  mcpTrksFound->setSubset( true ) ;
  mcpTrksFound->reserve(  mcpIt.size() ) ;
  evt->addCollection( mcpTrksFound , "MCParticleTracksFound"   ) ;



  while( MCParticle* mcp = mcpIt.next()  ){
    
    // streamlog_out( DEBUG ) << " ... huhu " << std::endl ;

    bool cut = true ;
    APPLY_CUT( DEBUG, cut,  std::abs( mcp->getCharge() )  > 0.5  ) ;
    
    APPLY_CUT( DEBUG, cut,  mcp->getGeneratorStatus() != 2   ) ;   // no documentation lines

    gear::Vector3D v( mcp->getVertex()[0], mcp->getVertex()[1], mcp->getVertex()[2] );
    gear::Vector3D e( mcp->getEndpoint()[0], mcp->getEndpoint()[1], mcp->getEndpoint()[2] );
    gear::Vector3D p( mcp->getMomentum()[0], mcp->getMomentum()[1], mcp->getMomentum()[2] );

    APPLY_CUT( DEBUG, cut, v.r() < 100.   ) ;   // start at IP+/-10cm

    APPLY_CUT( DEBUG, cut, e.rho()==0.0  || e.rho() > 400.   ) ; // end at rho > 40 cm

    APPLY_CUT( DEBUG, cut, p.rho() > 0.1  ) ; // pt> 100 MeV

    APPLY_CUT( DEBUG, cut, std::abs( cos( p.theta() ) )  < 0.99  ) ; //  | cos( theta ) | > 0.99

    //....

    if( cut ) 
      mcpTracks->push_back( mcp ) ;
  }
  
 

  LCIterator<MCParticle> trmIt( evt,  "MCParticleTracks" ) ;

  while( MCParticle* trm = trmIt.next()  ){
    
    double pxmcp =  trm->getMomentum()[0]  ;
    double pymcp =  trm->getMomentum()[1]  ;
    double pzmcp =  trm->getMomentum()[2]  ;
    
    double pmcp  = sqrt( pxmcp*pxmcp + pymcp*pymcp + pzmcp*pzmcp ) ;
    double ptmcp = sqrt( pxmcp*pxmcp + pymcp*pymcp ) ;
    
    double thmcp  = atan2( ptmcp , pzmcp ) ;

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
    
    const EVENT::LCObjectVec& trkV = nav.getRelatedFromObjects( trm ) ;
    
    if( trkV.size() >  0 ){
      
      //      LCUTIL::LCTypedVector<Track*> trks( trkV ) ;
      
      const FloatVec& wV = nav.getRelatedFromWeights( trm ) ;
      double wMax = 0.0 ;
      int iMax = 0 ;
      for(unsigned i=0 ; i<wV.size() ; ++i ) {
	  
	  if( wV[i] > wMax ){
	    wMax =  wV[i] ;
	    iMax = i ;
	  }  
	} 

      Track* tr = dynamic_cast<Track*>( trkV.at(iMax)  ) ; 

 
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
      APPLY_CUT( DEBUG4, cut,  std::abs( dom )  <  (5.*eom)  ) ;              // require 5 sigma on omega only for now
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
void TrackEfficiency::check( LCEvent * evt ) { 
  /*************************************************************************************************/

  }
//====================================================================================




void TrackEfficiency::end(){ 
  
  streamlog_out( MESSAGE )  << " processed " << _nEvt << " events in " << _nRun << " runs "
			    << std::endl ;
}

//====================================================================================================

