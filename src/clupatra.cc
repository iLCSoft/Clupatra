#include "clupatra.h"
#include "gear/GEAR.h"

#include <set>


namespace clupatra{


  struct CDot{
    CDot( int i , int j , double d) : I(i) , J(j) , Dot( d ) {}
    int I ;
    int J ;
    double Dot ;
  };

  bool sort_CDot( const CDot& c0 ,const CDot& c1 ){
    return c0.Dot > c1.Dot ;
  }

  // inline double maxIdx( double a, double b ){
  //   return ( a > b ) ? 0 : 1  ;
  // }
  // inline double maxIdx( double a, double b, double c){
  //   //    return ( a > b ) ?  maxIdx( a, c ) : max( b , c ) ;
  //   if( a > b ) {    return ( a > c ) ?  0 : 2 ;
  //   } else      {    return ( b > c ) ?  1 : 2 ;
  //   }
  // }

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

    streamlog_out(  DEBUG ) << " create_two_clusters  --- called ! " << std::endl ;

    gear::Vector3D lastDiffVec(0.,0.,0.) ;
    
    for( int l=maxLayerID-1 ; l >= 0 ; --l){

      GHitList& hL = hitsInLayer[ l ] ;
	
      streamlog_out(  DEBUG ) << " create_two_clusters  --- layer " << l  <<  " size: " << hL.size() << std::endl ;

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


      streamlog_out(  DEBUG ) << " create_two_clusters  ---  layer " << l 
			      << "  h0 : " << p0 
			      << "  h1 : " << p1 ;
      //			       << std::endl ;
      
      if( clu0->empty() ){ // first hit pair
	
	streamlog_out(  DEBUG ) << " create_two_clusters  --- initialize clusters " << std::endl ;
	
	clu0->addHit( h0 ) ;
	clu1->addHit( h1 ) ;
	
	lastDiffVec = p1 - p0 ;
	
	continue ;  
      }

      gear::Vector3D d = p1 - p0 ;
      
      float s0 =  ( lastDiffVec + d ).r() ;
      float s1 =  ( lastDiffVec - d ).r() ;
      
      if( s0 > s1 ){  // same orientation, i.e. h0 in this layer belongs to h0 in first layer
	
	streamlog_out(  DEBUG ) << " create_two_clusters  ---   same orientation " << std::endl ;
	clu0->addHit( h0 ) ;
	clu1->addHit( h1 ) ;
	
      } else{                // oposite orientation, i.e. h1 in this layer belongs to h0 in first layer

	streamlog_out(  DEBUG3 ) << " create_two_clusters  ---  oposite orientation " << std::endl ;
	clu0->addHit( h1 ) ;
	clu1->addHit( h0 ) ;
      }

    }

    // clu.freeHits() ;

    streamlog_out(  DEBUG1 ) << " create_two_clusters  --- clu0 " << clu0->size() 
			     <<  " clu1 " << clu1->size() << std::endl ;

    // return std::make_pair( clu0 , clu1 ) ;

    return ;
  }

  //------------------------------------------------------------------------------------------------------------------------- 

  void create_three_clusters( const GHitVec& hV, GClusterVec& cluVec,  int maxLayerID) {
    
    GHitListVector hitsInLayer( maxLayerID )  ; 
    //    addToHitListVector(  clu.begin(), clu.end(), hitsInLayer ) ;
    addToHitListVector(  hV.begin(), hV.end(), hitsInLayer ) ;
    
    GCluster* clu[3] ;

    gear::Vector3D lastp[3] ;
    gear::Vector3D cluDir[3] ;

    clu[0] = new GCluster ;
    clu[1] = new GCluster ;
    clu[2] = new GCluster ;
    
    cluVec.push_back( clu[0] ) ;
    cluVec.push_back( clu[1] ) ;
    cluVec.push_back( clu[2] ) ;

       
    for( int l=maxLayerID-1 ; l >= 0 ; --l){

      GHitList& hL = hitsInLayer[ l ] ;
	
      streamlog_out(  DEBUG ) << " create_three_clusters  --- layer " << l  <<  " size: " << hL.size() << std::endl ;

      if( hL.size() != 3 ){
	
	// if there is not exactly two hits in the layer, we simply copy the unassigned hits to the leftover hit vector
	// std::copy( hL.begin(),  hL.end() , std::back_inserter( hV ) ) ;
	continue ;
      }

      GHitList::iterator iH = hL.begin() ;

      GHit* h[3] ;
      h[0] = *iH++ ;
      h[1] = *iH++   ;
      h[2] = *iH   ;
      
      gear::Vector3D p[3] ;
      p[0] = h[0]->first->pos ;
      p[1] = h[1]->first->pos ;
      p[2] = h[2]->first->pos ;
      
      
      streamlog_out(  DEBUG ) << " create_three_clusters  ---  layer " << l 
			      << "  h0 : " << p[0] 
			      << "  h1 : " << p[1] 
			      << "  h2 : " << p[2] ;
      //			       << std::endl ;
      
      if( clu[0]->empty() ){ // first hit triplet
	
	streamlog_out(  DEBUG ) << " create_three_clusters  --- initialize clusters " << std::endl ;
	
	clu[0]->addHit( h[0] ) ;
	clu[1]->addHit( h[1] ) ;
	clu[2]->addHit( h[2] ) ;
	
	lastp[0] =  p[0] ;
	lastp[1] =  p[1] ;
	lastp[2] =  p[2] ;
	
	lastp[0] = ( 1. / p[0].r() ) * p[0] ;
	lastp[1] = ( 1. / p[1].r() ) * p[1] ;
	lastp[2] = ( 1. / p[2].r() ) * p[2] ;
	
	continue ;  
      }


      // unit vector in direction of current hits
      gear::Vector3D pu[3] ;
      pu[0] = ( 1. / p[0].r() ) * p[0] ;
      pu[1] = ( 1. / p[1].r() ) * p[1] ;
      pu[2] = ( 1. / p[2].r() ) * p[2] ;
      

      std::list< CDot > cDot ; 

      // create list of dot products between last hit in cluster and current hit 
      // cos( angle between  hits as seen from IP ) 
      for( int i = 0 ; i < 3 ; ++i ){	         
	for( int j = 0 ; j < 3 ; ++j ){
	  
	  // // unit vector in direction of last hit and 
	  // gear::Vector3D pu = - ( p[0] - lastp[0] ) ;
	  // pu[1] = - ( p[1] - lastp[1] ) ;
	  // pu[2] = - ( p[2] - lastp[2] ) ;
	  // pu[0] =  (  1. / pu[0].r() ) * pu[0]  ;
	  // pu[1] =  (  1. / pu[1].r() ) * pu[1]  ;
	  // pu[2] =  (  1. / pu[2].r() ) * pu[2]  ;
	  // cDot.push_back( CDot( i , j , cluDir[i].dot( pu[ j ] ) )) ;
	  
	  cDot.push_back( CDot( i , j , lastp[i].dot( pu[ j ] ) )) ;
	  
	}   
      }   
      
      // sort dot products in descending order 
      cDot.sort( sort_CDot ) ;

      // assign hits to clusters with largest dot product ( smallest angle )

      std::set<int> cluIdx ;
      std::set<int> hitIdx ;

      for( std::list<CDot>::iterator it = cDot.begin() ; it != cDot.end() ; ++ it ) {
	
	streamlog_out(  DEBUG ) << " I : " << it->I  
				<< " J : " << it->J 
				<< " d : " << it->Dot 
				<< std::endl ;
	
	int i =  it->I ;
	int j =  it->J ;

	// ignore clusters or hits that hvae already been assigned
	if( cluIdx.find( i ) == cluIdx.end()  && hitIdx.find( j ) == hitIdx.end() ){

	  cluIdx.insert( i ) ;
	  hitIdx.insert( j ) ;

	  clu[ i ]->addHit( h[ j ] ) ;

	  // cluDir[i] = - ( p[j] - lastp[i] ) ;
	  // cluDir[i] =  (  1. / cluDir[i].r() ) * cluDir[i]  ;

	  lastp[i] =  p[j] ;

	  streamlog_out(  DEBUG ) << " **** adding to cluster : " << it->I  
				  << " hit  : " << it->J 
				  << " d : " << it->Dot 
				  << std::endl ;
	}
	
      }

      // lastp[0] =  pu[0] ;
      // lastp[1] =  pu[1] ;
      // lastp[2] =  pu[2] ;
      // lastp[0] =  p[0] ;
      // lastp[1] =  p[1] ;
      // lastp[2] =  p[2] ;
      
      
    }


    streamlog_out(  DEBUG ) << " create_three_clusters  --- clu[0] " << clu[0]->size() 
			    <<  " clu[1] " << clu[1]->size() 
			    <<  " clu[2] " << clu[2]->size() 
			    << std::endl ;


    // need to fit the cluster and re-assign false hits ....



    return ;
  }
  //------------------------------------------------------------------------------------------------------------------------- 
 
  
}//namespace
