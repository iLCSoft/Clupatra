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
   //-----------------------------------------------------------------

  /** Find the nearest hits in the previous and next layers - (at most maxStep layers appart - NOT YET).
   */
  void  findNearestHits( GCluster& clu,  int maxLayerID, int maxStep){
    
    GHitListVector hitsInLayer( maxLayerID )  ; 
    addToHitListVector(  clu.begin(), clu.end(), hitsInLayer ) ;
    
    for(GCluster::iterator it= clu.begin() ; it != clu.end() ; ++it ){
      
      GHit* hit0 = *it ; 
      gear::Vector3D& pos0 =  hit0->first->pos ;
      
      int layer = hit0->first->layerID ;
      
      int l=0 ;
      if( (l=layer+1) <  maxLayerID ) {
	
	GHitList& hL = hitsInLayer[ l ] ;
	
	double minDist2 = 1.e99 ;
	GHit*   bestHit = 0 ;
	
	for( GHitList::iterator iH = hL.begin() ; iH != hL.end() ; ++iH ) {
	  
	  GHit* hit1 = *iH ; 
	  gear::Vector3D& pos1 = hit1->first->pos ;
	  
	  double dist2 = ( pos0 - pos1 ).r2() ;
	  
	  if( dist2 < minDist2 ){
	    minDist2 = dist2 ;
	    bestHit = hit1 ;
	  }
	}
	if( bestHit ){
	  
	  hit0->first->nNNHit =  bestHit->first ;
	  hit0->first->nDist  =  minDist2 ;
	}
      }
      
      if( (l=layer-1) > 0  ) {
	
	GHitList& hL = hitsInLayer[ l ] ;
	
	double minDist2 = 1.e99 ;
	GHit* bestHit = 0 ;
	
	for( GHitList::iterator iH = hL.begin() ; iH != hL.end() ; ++iH ) {
	  
	  GHit* hit1 = *iH ; 
	  gear::Vector3D&  pos1 =  hit1->first->pos ;
	  
	  double dist2 = ( pos0 - pos1 ).r2() ;
	  
	  if( dist2 < minDist2 ){
	    minDist2 = dist2 ;
	    bestHit = hit1 ;
	  }
	}
	if( bestHit ){
	  
	  hit0->first->pNNHit =  bestHit->first ;
	  hit0->first->pDist  =  minDist2 ;
	}
      }
      
    }
  }


  //------------------------------------------------------------------------------------------------------------------------- 
  //std::pair<GCluster*, GCluster* > create_two_clusters( GCluster& clu,  GHitVec& hV,  int maxLayerID){
  void create_two_clusters( const GHitVec& hV, GClusterVec& cluVec,  int maxLayerID) {

    GHitListVector hitsInLayer( maxLayerID )  ; 
    //    addToHitListVector(  clu.begin(), clu.end(), hitsInLayer ) ;
    addToHitListVector(  hV.begin(), hV.end(), hitsInLayer ) ;
    

    GCluster* clu0 = new GCluster ;
    GCluster* clu1 = new GCluster ;

    cluVec.push_back( clu0 ) ;
    cluVec.push_back( clu1 ) ;

    streamlog_out(  DEBUG3 ) << " create_two_clusters  --- called ! " << std::endl ;

    gear::Vector3D lastDiffVec(0.,0.,0.) ;
    
    for( int l=maxLayerID-1 ; l >= 0 ; --l){

      GHitList& hL = hitsInLayer[ l ] ;
	
      streamlog_out(  DEBUG3 ) << " create_two_clusters  --- layer " << l  <<  " size: " << hL.size() << std::endl ;

      if( hL.size() != 2 ){
	
	// if there is not exactly two hits in the layer, we simply copy the unassigned hits to the leftover hit vector
	// std::copy( hL.begin(),  hL.end() , std::back_inserter( hV ) ) ;
	continue ;
      }

      GHitList::iterator iH = hL.begin() ;

      GHit* h0 = *iH++ ;
      GHit* h1 = *iH   ;
      
      gear::Vector3D& p0 = h0->first->pos ;
      gear::Vector3D& p1 = h1->first->pos ;


      streamlog_out(  DEBUG3 ) << " create_two_clusters  ---  layer " << l 
			       << "  h0 : " << p0 
			       << "  h1 : " << p1 ;
      //			       << std::endl ;
      
      if( clu0->empty() ){ // first hit pair
	
	streamlog_out(  DEBUG3 ) << " create_two_clusters  --- initialize clusters " << std::endl ;
	
	clu0->addHit( h0 ) ;
	clu1->addHit( h1 ) ;
	
	lastDiffVec = p1 - p0 ;
	
	continue ;  
      }

      gear::Vector3D d = p1 - p0 ;
      
      float s0 =  ( lastDiffVec + d ).r() ;
      float s1 =  ( lastDiffVec - d ).r() ;
      
      if( s0 > s1 ){  // same orientation, i.e. h0 in this layer belongs to h0 in first layer
	
	streamlog_out(  DEBUG3 ) << " create_two_clusters  ---   same orientation " << std::endl ;
	clu0->addHit( h0 ) ;
	clu1->addHit( h1 ) ;
	
      } else{                // oposite orientation, i.e. h1 in this layer belongs to h0 in first layer

	streamlog_out(  DEBUG3 ) << " create_two_clusters  ---  oposite orientation " << std::endl ;
	clu0->addHit( h1 ) ;
	clu1->addHit( h0 ) ;
      }

    }

    // clu.freeHits() ;

    streamlog_out(  DEBUG3 ) << " create_two_clusters  --- clu0 " << clu0->size() 
			     <<  " clu1 " << clu1->size() << std::endl ;

    // return std::make_pair( clu0 , clu1 ) ;

    return ;
  }

  //------------------------------------------------------------------------------------------------------------------------- 
 
  
}//namespace
