#include "clupatra.h"
#include "gear/GEAR.h"

namespace clupatra{


  void addHitsAndFilter( GCluster* clu, GHitListVector& hLV , double dChi2Max, double chi2Cut, unsigned maxStep, bool backward) {
    
    Chi2_RPhi_Z_Hit ch2rzh ;

    KalTrack* trk =  clu->ext<ClusterInfo>()->track ;
    
    unsigned step = 0 ;
    
    while( step < maxStep + 1 ) {
      
      gear::Vector3D xv ;
      int layer ;
      
      bool hitAdded = false ;
      
      trk->findNextXingPoint(  xv , layer , step , backward ) ;
      
      streamlog_out( DEBUG2 ) <<  "  -- addHitsAndFilter(): searching in leftover hits for cluster : " << std::endl 
			      <<  "  omega : " <<  trk->getOmega()  
			      <<  "  Step : " << step 
			      <<  "  at layer: "   << layer         
			      <<  "  next xing point : " <<  xv  ;
      
      if( layer > -1 ) { // found a crossing point 
	
 	GHitList& hLL = hLV.at( layer ) ;
	
	double ch2Min = 1.e99 ;
	GHit* bestHit = 0 ;
	
 	for( GHitList::const_iterator ih = hLL.begin() ; ih != hLL.end() ; ++ih ){    

	  GHit* hit = *ih ;

	  double ch2 = ch2rzh( hit->first , xv )  ;
	  
 	  if( ch2 < ch2Min ){
	    
 	    ch2Min = ch2 ;
 	    bestHit = hit ;
 	  }
 	}
	
 	if( bestHit != 0 ){
	  
 	  Vector3D hPos(  bestHit->first->lcioHit->getPosition() ) ;
	  
 	  bestHit->first->chi2Residual = ch2Min ;
	  
	  if( ch2Min  < chi2Cut ) { 
	    
	    streamlog_out( DEBUG2 ) <<   " ---- assigning left over hit : " << hPos << " <-> " << xv
				   <<   " dist: " <<  (  hPos - xv ).r()
				   <<   " chi2: " <<  ch2Min 
				   <<   "  hit errors :  rphi=" <<  sqrt( bestHit->first->lcioHit->getCovMatrix()[0] 
									  + bestHit->first->lcioHit->getCovMatrix()[2] ) 
				   <<	 "  z= " <<  sqrt( bestHit->first->lcioHit->getCovMatrix()[5] )
				   << std::endl ;
	    

	    double deltaChi2 =  trk->testDeltaChi2(  bestHit->first->fitHit ) ;
	    
	    bestHit->first->deltaChi2 = deltaChi2 ;
	    
	    if(   deltaChi2 < dChi2Max   &&  trk->addAndFilter( bestHit->first->fitHit )  ){
	      
	      hitAdded = true ;
	      
	      hLL.remove(  bestHit ) ;
	      clu->addHit( bestHit ) ;
	      
	      backward = false ; 
	      // after we add the first (backward) hit, we have effectively switched the 
	      // direction of the track fit, i.e. the next search is 'forward'

	      streamlog_out( DEBUG2 ) <<   " ---- track state filtered with new hit ! ------- " << std::endl ;
	    }
	  } // chi2Cut 
	} // bestHit
      } // layer > -1 	

      if( hitAdded ){     
	step = 1 ;
      } else {  	
	++step  ;
      }
      
    } // while step < maxStep
  

  }
  //------------------------------------------------------------------------------------------------------------
  
  
  
}//namespace
