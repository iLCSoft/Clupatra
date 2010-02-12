// -*- C++ -*-
#ifndef LCIterator_include
#define LCIterator_include

#include <string>
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
  }
  
  T* next(){

    if( _i < _n ) 
      return (T*) _col->getElementAt( _i++ ) ;
    else
      return 0 ;
  }

protected:
  int _n, _i ;
  EVENT::LCCollection* _col ;
} ;



#endif
