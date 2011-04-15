// -*- C++ -*-
#ifndef LCIterator_include
#define LCIterator_include

#include <string>
#include <sstream>
#include "EVENT/LCEvent.h"
#include "EVENT/LCCollection.h"


/** Test class for LCCollectionIterator - should save some boiler plate code....
 */

template <class T>
class LCIterator{

  LCIterator<T>() {} 

public:
  LCIterator<T>( EVENT::LCEvent* evt, const std::string& name ) : _i(0) {
    
    _col = evt->getCollection( name ) ;
    
    _n = _col->getNumberOfElements() ;

    if( _n > 0 ){
      T* t = dynamic_cast<T*>(  _col->getElementAt(0)  );
      if( t == 0 ){
	std::stringstream s ;
	s << " invalid iterator type  : " << typeid( t ).name() << " for collection " <<  name  << std::endl ; 
	throw lcio::Exception( s.str() ) ;
      }
    }

  }
  
  T* next(){

    if( _i < _n ) 
      return (T*) _col->getElementAt( _i++ ) ;
    else
      return 0 ;
  }

  int size() { return _n ; }

protected:
  int _n, _i ;
  EVENT::LCCollection* _col ;
} ;



#endif
