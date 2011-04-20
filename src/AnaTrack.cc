#include "AnaTrack.h"


#include "marlin/AIDAProcessor.h"

//#include <time.h>
#include <vector>
#include <map>
#include <algorithm>
#include <math.h>
#include <cmath>

//---- LCIO ---
#include "EVENT/LCCollection.h"
#include "EVENT/Track.h"
#include "UTIL/Operators.h"
#include "UTIL/LCRelationNavigator.h"

//---- GEAR ----
// #include "marlin/Global.h"
// #include "gear/GEAR.h"
// #include "gear/TPCParameters.h"
// #include "gear/PadRowLayout2D.h"
// #include "gear/BField.h"

//---- ROOT -----
#include "TH1F.h" 

#include "LCIterator.h"

using namespace lcio ;
using namespace marlin ;


AnaTrack aAnaTrack ;

// helper enum defining histogram index in vector 
namespace AnaTrackHistos {
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
}
using namespace AnaTrackHistos ;


//======================================================================================================

#define APPLY_CUT( Cut, Exp ) Cut &= ( Exp ) ;				\
  if( ! Cut )  streamlog_out( DEBUG4 ) << "  ***** failed cut:  [ " <<  #Exp \
                                       <<  " ] in evt: " << evt->getEventNumber() \
				       << " run: "  << evt->getRunNumber()   << std::endl ;

//======================================================================================================

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


//======================================================================================================


AnaTrack::AnaTrack() : Processor("AnaTrack") {
  
  // modify processor description
  _description = "AnaTrack: analysis plots for Tracks" ;
  
  
  
  registerInputCollection( LCIO::TRACK,
			   "TrackCollection" , 
			   "Name of the input collections with Tracks"  ,
			   _trkColName ,
			   std::string("KalTestTracks" ) ) ;

  registerInputCollection( LCIO::LCRELATION,
			   "TrackMCPRelation" , 
			   "Name of the input collections with LCRelation Tracks to MCParticles"  ,
			   _relColName ,
			   std::string("TrackRelation" ) ) ;
  
  FloatVec ptRange ;
  ptRange.push_back( 45. ) ;
  ptRange.push_back( 55. ) ;

  registerProcessorParameter( "PtRange" , 
			      "Min and max value of pt range [GeV]"  ,
			      _ptRange ,
			      ptRange ) ;
  
}


void AnaTrack::init() { 

  // usually a good idea to
  printParameters() ;


  if( _ptRange.size() != 2 ) {

    throw Exception("AnaTrack::init()  processor parameter 'ptRange' needs to have two values (min and max pt in GeV)" ) ;
  }


  _nRun = 0 ;
  _nEvt = 0 ;
}

void AnaTrack::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 

void AnaTrack::processEvent( LCEvent * evt ) { 


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
    
    _h1.resize( AnaTrackHistos::Size ) ;
    
    static const int nBin = 50 ;

    //  define binning for pt and omega
    double ptBinL  = ptMin - ( ptMax - ptMin ) / 10.  ;
    double ptBinH  = ptMax + ( ptMax - ptMin ) / 10.  ;
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

  LCRelationNavigator nav( evt->getCollection( _relColName )  ) ;


  LCIterator<Track> trIt( evt, _trkColName ) ;

  while( Track* tr = trIt.next()  ){

    double d0 = tr->getD0() ;
    double ph = tr->getPhi() ;
    double om = tr->getOmega() ;
    double z0 = tr->getZ0() ;
    double tL = tr->getTanLambda() ;
    double pt = alpha / std::abs( om ) ;

    double ed0 = sqrt( tr->getCovMatrix()[0]  );
    double eph = sqrt( tr->getCovMatrix()[2]  );
    double eom = sqrt( tr->getCovMatrix()[5]  );
    double ez0 = sqrt( tr->getCovMatrix()[9]  );
    double etL = sqrt( tr->getCovMatrix()[14] );

    double ept = alpha * eom / ( om * om )  ; // ept = al/om^2 * eom 
    
    double d0mcp = 0 ;
    double phmcp = 0 ;
    double ommcp = 0 ; 
    double z0mcp = 0 ; 
    double tLmcp = 0 ; 
    double ptmcp = 0 ;

    const EVENT::LCObjectVec& mcpV = nav.getRelatedToObjects( tr ) ;
    
    if( mcpV.size() >  0 ){

      MCParticle* mcp = dynamic_cast<MCParticle*>( mcpV.at(0)  ) ; 

      double px =  mcp->getMomentum()[0]  ;
      double py =  mcp->getMomentum()[1]  ;
      double pz =  mcp->getMomentum()[2]  ;

      //      double pmcp  = sqrt( px*px + py*py + pz*pz ) ;

      ptmcp = sqrt( px*px + py*py ) ;

      d0mcp = 0. ;  //fixme: only true for prompt tracks
      phmcp = atan2( py, px) ;

      ommcp = alpha / ptmcp  ;
      if(  mcp->getCharge() < 0. ) 
	ommcp = -ommcp ;

      z0mcp = 0 ;  //fixme: only true for prompt tracks
      tLmcp = pz / ptmcp  ; 


    } else {
      
      streamlog_out( DEBUG4 ) << " no MCParticle found for track in evt: " << evt->getEventNumber() 
			      << " run : " << evt->getRunNumber()  << std::endl 
			      <<  tr->id() 
			      << std::endl ;
    }
    
    bool cut = true ;

    APPLY_CUT( cut,  mcpV.size() == 1 ) ;

    APPLY_CUT( cut,  tr->getChi2() / tr->getNdf() < 3. ) ;

    APPLY_CUT( cut, (ptMin < ptmcp && ptmcp < ptMax ) ) ;


    if( cut == true ){

      double dd0 = d0 - d0mcp ;
      double dph = ph - phmcp ; 
      if( std::abs( dph )  > M_PI )
	dph = 2.* M_PI - std::abs( dph ) ;

      double dom = om - ommcp ; 
      double dz0 = z0 - z0mcp ; 
      double dtL = tL - tLmcp ; 
      double dpt = pt - ptmcp ; 

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

      h.fill( hd0mcp,    d0mcp ) ;
      h.fill( hphimcp,   phmcp ) ;
      h.fill( homegamcp, ommcp ) ;
      h.fill( hz0mcp,    z0mcp ) ;
      h.fill( htanLmcp,  tLmcp ) ;
      h.fill( hptmcp,    ptmcp ) ;

   
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


  //----------------------------------------------------------------------  
  
  _nEvt ++ ;

  clock_t end = clock () ; 
  
  // streamlog_out( DEBUG )  << "---  histogramming time: " 
  // 			  <<  double( end - start ) / double(CLOCKS_PER_SEC) << std::endl  ;

}




/*************************************************************************************************/
void AnaTrack::check( LCEvent * evt ) { 
  /*************************************************************************************************/

  }
//====================================================================================




void AnaTrack::end(){ 
  
  streamlog_out( MESSAGE )  << " processed " << _nEvt << " events in " << _nRun << " runs "
			    << std::endl ;
}

//====================================================================================================

