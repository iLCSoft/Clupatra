/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#ifndef NNClusterer_h
#define NNClusterer_h 1

#include <list>
#include <vector>
#include <map>

#include "LCRTRelations.h"

/** Nearest neighbour type clusering for arbitrary types.
 *
 *  @author F.Gaede (DESY)
 *  @version $Id: NNClusterer.h 2468 2011-09-06 07:55:04Z gaede $
 */
namespace nnclu {
  
  
  
  // forward declaration:
  template <class U>
  class Cluster ;


  /** Wrapper class for elements that are clustered, holding a pointer to the actual 
   *  object "->first"  and a pointer to the cluster this obejct belongs to "->second". 
   *  
   *  @see Cluster
   *  @author F.Gaede (DESY)
   *  @version $Id:$
   */
  template <class T>
  class Element : public  std::pair< T*, Cluster<T>* >{
  
    typedef T value_type ;
    typedef Cluster<T> cluster_type ;

  public:
  
    /** Default c'tor takes a pointer to the original element type object. The optional index can be used to 
     *  code nearest neighbour bins, e.g. in z-coordinate to speed up the clustering process.
     */
    Element(T* element, int index0 = 0 ) : Index0( index0 )   {  
      Pair::first = element ; 
      Pair::second = 0 ;
    }
  
    /** C'tor that also takes a pointer to the cluster this element belongs to - in case seed elements/clusters are used.
     */
    Element(T* element ,  Cluster<T>* cl , int index0 = 0) : Index0( index0 ) {
      Pair::first =  element ;
      Pair::second = cl ;
    }
  
    /** Index that can be used to code nearest neighbour bins, e.g. in z-coordinate
     *  to speed up the clustering process.
     */
    int Index0 ;
  
  protected:
    typedef std::pair< T*, Cluster<T>* > Pair ;

    /** Don't allow default c'tor w/o element */
    Element() ;
  } ;


  /** Helper class that creates an Elements for an objects of type T.
   */
  template <class T>
  struct MakeElement{  
    Element<T>  operator()( T* t) { return new Element<T>( t ) ;  }    
  } ;


  /** Extension of std::vector that allows to take ownership of objects pointed to and
   *  delete these when going out of scope.
   */
  template <class T> 
  class PtrVector : public std::vector<T*> {
    typedef std::vector<T*> vec ;
    bool _isOwner ;
  public:
    PtrVector() : _isOwner( false ) {}
    ~PtrVector() {
      if( _isOwner ) 
        for( typename vec::iterator i = vec::begin(),end = vec::end(); i != end ; delete *i++ ) ; //++i ) delete *i ;
    }
    void setOwner( bool val=true ) { _isOwner = val ; }
  };
  

  /** Extension of std::list that allows to take ownership of objects pointed to and
   *  delete these when going out of scope.
   */
  template <class T> 
  class PtrList : public std::list<T*> {
    typedef std::list<T*> vec ;
    bool _isOwner ;
  public:
    PtrList() : _isOwner( false ) {}
    ~PtrList() {
      if( _isOwner ) 
        for( typename vec::iterator i = vec::begin(),end = vec::end(); i != end ; delete *i++ ) ; //++i ) delete *i ;
    }
    void setOwner( bool val=true ) { _isOwner = val ; }
  };
  
  

  /** Templated class for generic clusters  of Elements that are clustered with
   *  an NN-like clustering algorithm. Effectively this is just a list of elements.
   * 
   *  @see Element
   *  @author F.Gaede (DESY)
   *  @version $Id:$
   */
  template <class T >
  class Cluster : public std::list< Element<T> * >, public lcrtrel::LCRTRelations {
  
  public :
    typedef Element<T> element_type ; 
  
    int ID ; //DEBUG
  
    Cluster()  {}
  
    /** C'tor that takes the first element */
    Cluster( Element<T>* element)  {
      static int SID ;  //DEBUG
      ID = SID++ ;      //DEBUG
      addElement( element ) ;
    }
  
    /** Add a element to this cluster - updates the element's pointer to cluster */
    void addElement( Element<T>* element ) {
    
      element->second = this ;
      push_back( element ) ;
    }

    // /** Remove all elements from the cluster and reset the cluster association, i.e. elements can be
    //  *  used for another clustering procedure.
    //  */
    // template <class Out>
    // void takeElements(Out result){
    //   typename Cluster<T>::iterator it = this->begin() ;
    //   while( it !=  this->end() ){
    //     (*it)->second = 0 ;
    //     result++ = *it ;
    //     it = this->erase(it) ;
    //   }
    // }


    /** Free all elements, ie. reset their cluster association - the elements are still in the list !.
     */
    void freeElements(){

      for( typename Cluster<T>::iterator it = this->begin(), end = this->end() ; it != end ; it++ ){
        (*it)->second = 0 ;
      }
    
      // typename Cluster<T>::iterator it = this->begin() ;
      // while( it !=  this->end() ){
      //   (*it)->second = 0 ;
      //   it = this->erase(it) ;
      // }
    }
  
    /** Merges all elements from the other cluster cl into this cluster */
    void mergeClusters( Cluster<T>* cl ) {
    
      for( typename Cluster<T>::iterator it = cl->begin(), end = cl->end() ; it != end ; it++ ){
        (*it)->second = this  ;
      }
      this->merge( *cl ) ;
    }
  
    /** D'tor frees all remaining elements that still belong to this cluster */
    ~Cluster()  {
    
      //  typename Cluster<T>::iterator it = this->begin() ;
      //  while( it !=  this->end()  )
    
      for( typename Cluster<T>::iterator it = this->begin() , end =  this->end()  ;  it != end ; ++it ){
      
        typename Cluster<T>::value_type h = *it ; 
      
        if( h != 0 && h->second == this )
          h->second = 0 ;
      
      }
    }
  } ;


  template <class T>
  /** Main class for a nearest neighbour type clustering. 
   * 
   *  @author F.Gaede (DESY)
   *  @version $Id:$
   */
  class NNClusterer{

  public:
    typedef T value_type ;
    typedef Cluster<T> cluster_type ; 
    typedef Element<T> element_type ; 
    typedef PtrVector< element_type >    element_vector ;
    typedef PtrVector< cluster_type >    cluster_vector ;
    typedef PtrList< element_type >    element_list ;
    typedef PtrList< cluster_type >    cluster_list ;
    
    /** Simple nearest neighbour (NN) clustering algorithm. Users have to provide an input iterator of
     *  Element objects and an output iterator for the clusters found. The predicate has to have 
     *  a method with the following signature: bool operator()( const Element<T>*, const Element<T>*).
     *  All pairs of elements for which this method returns 'true' will be merged into one output cluster 
     *  - all other pairs of elements will be in different clusters.
     */
  
    template <class In, class Out, class Pred > 
    void cluster( In first, In last, Out result, Pred& pred , const unsigned minSize=1) {
    
      // typedef typename Pred::element_type ElementType ;
    
      cluster_vector tmp ; 
      tmp.reserve( 1024 ) ;
      
      while( first != last ) {
      
        for( In other = first+1 ;   other != last ; other ++ ) {
        
          if( pred( (*first) , (*other) ) ) {
          
            if( (*first)->second == 0 && (*other)->second == 0 ) {  // no cluster exists
            
              //             Cluster<ElementType >* cl = new Cluster<ElementType >( (*first) ) ;
              cluster_type* cl = new cluster_type( (*first) ) ;
            
              cl->addElement( (*other) ) ;
            
              tmp.push_back( cl ) ;
            
            }
            else if( (*first)->second != 0 && (*other)->second != 0 ) { // two clusters
            
              if(  (*first)->second != (*other)->second )  // don't call merge on identical clusters
                (*first)->second->mergeClusters( (*other)->second ) ;
            
            } else {  // one cluster exists
            
              if( (*first)->second != 0 ) {
              
                (*first)->second->addElement( (*other)  ) ;
              
              } else {                           
              
                (*other)->second->addElement( (*first)  ) ;
              }
            }
          
          } // dCut 
          //       ++j ;
        }
        //     ++i ;
        ++first ;
      }
    
      // remove empty clusters 
      //   std::remove_copy_if( tmp.begin() , tmp.end() , result , &empty_list< Cluster<ElementType > > ) ;
    
      for( typename cluster_vector::iterator i = tmp.begin(); i !=  tmp.end() ; i++ ){
      
        if( (*i)->size() > minSize-1 ) {
        
          result++ = *i ;
        } 
        else {  
        
          delete *i ; 
        }
      }
    }


  };



  /**Splits a list into two based on a predicate. The new list will 
   * hold all elements for which Pred is true which are in turn removed 
   * from the original list.
   */
  template <class List, class Out, class Pred > 
  void split_list( List& list, Out result, Pred pred) {
    
    typename List::iterator it = list.begin() ;
    
    while(  it != list.end() ){
      
      if( pred( *it ) ){

        result++ = *it ;
        it = list.erase( it ) ;

      } else
        ++it ;
    }
  }
  
  

  // /** Same as above - but requires the elements to be sorted in index0 (only compare neighbouring bins in index0). */
  // template <class In, class Out, class Pred > 
  // void cluster_sorted( In first, In last, Out result, Pred* pred , const unsigned minSize=1) {
  
  //   typedef typename In::value_type ElementPtr ;
  //   typedef typename Pred::element_type ElementType ;
  
  //   typedef std::vector< Cluster<ElementType >* >  ClusterList ;
  
  //   ClusterList tmp ; 
  //   tmp.reserve( 1024 ) ;
  
  //   //   int i(0),j(0) ;

  //   while( first != last ) {

  //     //     j=i+1 ;

 
  //     for( In other = first+1 ;   other != last ; other ++ ) {

      
  //       if( notInRange<-1,1>(   (*first)->Index0 - (*other)->Index0  )   ) 
  //          break ;


  //       if( pred->mergeElements( (*first) , (*other) ) ) {
	
  //         if( (*first)->second == 0 && (*other)->second == 0 ) {  // no cluster exists
	  
  //           Cluster<ElementType >* cl = new Cluster<ElementType >( (*first) ) ;

  //           cl->addElement( (*other) ) ;
	  
  //           tmp.push_back( cl ) ;

  //         }
  //         else if( (*first)->second != 0 && (*other)->second != 0 ) { // two clusters
	  
  //           if(  (*first)->second != (*other)->second )  // don't call merge on identical clusters
  //             (*first)->second->mergeClusters( (*other)->second ) ;
	  
  //         } else {  // one cluster exists
	  
  //           if( (*first)->second != 0 ) {
	    
  //             (*first)->second->addElement( (*other)  ) ;
	    
  //           } else {                           
	    
  //             (*other)->second->addElement( (*first)  ) ;
  //           }
  //         }
	
  //       } // dCut 
  //       //       ++j ;
  //     }
  //     //     ++i ;
  //     ++first ;
  //   }

  //   // remove empty clusters 
  //   //   std::remove_copy_if( tmp.begin() , tmp.end() , result , &empty_list< Cluster<ElementType > > ) ;

  //   for( typename ClusterList::iterator i = tmp.begin(); i !=  tmp.end() ; i++ ){

  //     if( (*i)->size() > minSize-1 ) {

  //       result++ = *i ;
  //     } 
  //     else {  

  //       delete *i ; 
  //     }
  //   }
  // }





  // /** Helper method that copies all element pointers from an LCIO collection that fullfill the predicate to
  //  *  a ElementVec. The predicate can either be a bool funtion or functor that takes a T*, e.g.
  //  *  @see EnergyCut
  //  */
  // template <class T, class Pred> 
  // void addToElementVec(ElementVec<T>& v, EVENT::LCCollection* col, Pred pred ){

  //   for( int i=0 ; i < col->getNumberOfElements() ; i++ ){

  //     T* element = dynamic_cast<T*>( col->getElementAt( i) ) ;

  //     if( pred( element ) ){

  //       v.push_back( new Element<T>( element ) ) ;
  //     }
  //   }
  // } 
  // /** Same as addToElementVec(ElementVec<T>& v, LCCollection* col, Pred pred ) except that
  //  *  an additional order function/functor can be given that defines the index of the element, e.g.
  //  *  @see ZIndex.
  //  */
  // template <class T, class Pred, class Order> 
  // void addToElementVec(ElementVec<T>& v, EVENT::LCCollection* col, Pred pred , Order order ){

  //   for( int i=0 ; i < col->getNumberOfElements() ; i++ ){

  //     T* element = dynamic_cast<T*>( col->getElementAt( i) ) ;

  //     if( pred( element ) ){

  //       v.push_back( new Element<T>( element , order(element) ) ) ;
  //     }
  //   }
  // } 

  // /** Same as addToElementVec(ElementVec<T>& v, LCCollection* col, Pred pred ) except that
  //  *  no LCIO collection but an iterator pair is given. If no order index is needed, use NullIndex.
  //  */
  // template <class T, class Pred, class Order, class In> 
  // void addToElementVec(ElementVec<T>& v, In first, In last, Pred pred , Order order ){

  //   while( first != last ){
    
  //     T* element = *first++  ;

  //     if( pred( element) ){
  //       v.push_back( new Element<T>( element , order(element) ) ) ;
  //     }
  //   }
  // } 
  // // /** Helper method that copies all element pointers from a container of LCIO objects that fullfill the predicate to
  // //  *  a ElementVec. The predicate can either be a bool funtion or functor that takes a T*, e.g.
  // //  *  @see EnergyCut
  // //  */
  // // template <class T, class IT, class Pred> 
  // // void copyToElementVec(ElementVec<T>& v, IT start, IT end,  Pred pred ){
  // //   for(IT it = start; it != end ; ++it ){
  // //     T* element = (T*) *(it) ;
  // //     if( pred( element ) ){
  // //       v.push_back( new Element<T>( element ) ) ;
  // //     }
  // //   } 
  // // }



  // /** Helper vector of Cluster<T> taking care of memory management.
  //  * It is the users responsibility that various instances of this class do not
  //  * contain duplicate entries and that no elements are removed from the list
  //  * w/o being deleted or added to another list.
  //  */
  // template <class T> 
  // class cluster_vector : public std::list< Cluster<T>* > {
  //   typedef std::list< Cluster<T>* > List ;
  // public:
  //   ~cluster_vector() {
  //     //static int c=0 ;
  //     for( typename cluster_vector::iterator i = List::begin() ; i != List::end() ; i++) {
      
  // //       std::cout << " --- deleting cluster at : " << std::hex << *i << ", " 
  // //                 << std::dec<< c++ 
  // //                 << " - size: " << (*i)->size()
  // //                 << " ID: " <<  (*i)->ID 
  // //                 << std::endl ;
  //       delete *i ;
  //     }
  //   }

  //   void clear() {
  // //     for( typename cluster_vector::iterator i = List::begin() ; i != List::end() ; i++) {
  // //       delete *i ;
  // //     }
  //     std::list< Cluster<T>* >::clear() ;
  //   }
  // };


  // /** Simple predicate class for NN clustering. Requires 
  //  *  PosType* ElementClass::getPosition(), e.g for CalorimeterElements use: <br>
  //  *  NNDistance<CalorimeterElement,float> dist( myDistCut ) ; 
  //  */
  // template <class ElementClass, typename PosType > 
  // class NNDistance{
  // public:

  //   /** Required typedef for cluster algorithm 
  //    */
  //   typedef ElementClass element_type ;

  //   /** C'tor takes merge distance */
  //   NNDistance(float dCut) : _dCutSquared( dCut*dCut ) , _dCut(dCut)  {} 


  //   /** Merge condition: true if distance  is less than dCut given in the C'tor.*/ 
  //   inline bool mergeElements( Element<ElementClass>* h0, Element<ElementClass>* h1){
    
  //     if( std::abs( h0->Index0 - h1->Index0 ) > 1 ) return false ;

  //     const PosType* pos0 =  h0->first->getPosition() ;
  //     const PosType* pos1 =  h1->first->getPosition() ;
    
  //     return 
  //       ( pos0[0] - pos1[0] ) * ( pos0[0] - pos1[0] ) +
  //       ( pos0[1] - pos1[1] ) * ( pos0[1] - pos1[1] ) +
  //       ( pos0[2] - pos1[2] ) * ( pos0[2] - pos1[2] )   
  //       < _dCutSquared ;
  //   }
  
  // protected:
  //   NNDistance() ;
  //   float _dCutSquared ;
  //   float _dCut ;
  // } ;


  // /** Simple predicate class for applying an energy cut to the objects of type T.
  //  *  Requires float/double T::getEnergy().
  //  */
  // template <class T>
  // class EnergyCut{
  // public:
  //   EnergyCut( double eCut ) : _eCut( eCut ) {}  

  //   inline bool operator() (T* element) {  return element->getEnergy() > _eCut ; }

  // protected:

  //   EnergyCut() {} ;
  //   double _eCut ;
  // } ;


  // /** Simple predicate class for computing an index from N bins of the z-coordinate of LCObjects
  //  *  that have a float/double* getPostion() method.
  //  */
  // template <class T, int N>
  // class ZIndex{
  // public:
  //   /** C'tor takes zmin and zmax - NB index can be negative and larger than N */
  //   ZIndex( float zmin , float zmax ) : _zmin( zmin ), _zmax( zmax ) {}  

  //   inline int operator() (T* element) {  

  //     return (int) std::floor( ( element->getPosition()[2] - _zmin ) / ( _zmax - _zmin ) * N ) ; 
  //   }

  // protected:

  //   ZIndex() {} ;
  //   float _zmin ;
  //   float _zmax ;
  // } ;

  // /** Simple predicate class for 'computing' a NULL index.
  //  */
  // //template <class T>
  // struct NullIndex{
  //   inline int operator() (void* element) {  
  //     return 0 ;
  //   }
  // };

  // struct AllwaysTrue{
  //   inline int operator() (void* element) {  return true ; }  
  // };

  // /** Helper class that creates an lcio::Cluster from a generic cluster with element types that have a 
  //  *  getPosition() and a getEnergy() method.
  //  */
  // template <class T>
  // struct LCIOCluster{

  //   inline lcio::Cluster* operator() (Cluster<T>* c) {  
    
  //     IMPL::ClusterImpl* clu = new IMPL::ClusterImpl ;
    
  //     unsigned n = c->size() ;
  //     unsigned i=0 ;

  //     float a[n], x[n], y[n], z[n] ;

  //     for( typename Cluster<T>::iterator hi = c->begin(); hi != c->end() ; hi++) {
      
  //       T* element = (*hi)->first ;  

  //       a[i] = element->getEnergy() ;
  //       x[i] = element->getPosition()[0] ;
  //       y[i] = element->getPosition()[1] ;
  //       z[i] = element->getPosition()[2] ;

  //       clu->addElement( element , a[i] ) ;
      
  //       ++i ;
  //     }
    
  //     ClusternShapes cs( n,a,x,y,z) ;

  //     clu->setEnergy( cs.getTotalAmplitude()  ) ;
  //     clu->setPosition( cs.getCenterOfGravity() ) ;

  //     // direction of cluster's PCA
  //     float* d = cs.getEigenVecInertia() ;

  //     Hep3Vector v( d[0], d[1], d[2] ) ;
    
  //     clu->setITheta( v.theta() )  ;
  //     clu->setIPhi(   v.phi() ) ;
  
  //     std::vector<float> param(5) ;

  //     //     float* ev = cs.getEigenValInertia() ;
  //     //     param[0] = ev[0] ;
  //     //     param[1] = ev[1] ;
  //     //     param[2] = ev[2] ;

  //     param[0] = cs.getElipsoid_r1() ;
  //     param[1] = cs.getElipsoid_r2() ;
  //     param[2] = cs.getElipsoid_r3() ;
  //     param[3] = cs.getElipsoid_vol() ;
  //     param[4] = cs.getWidth() ;

  //     clu->setShape( param ) ;

  //     return clu ;

  //   }

  // } ;

}


#endif


