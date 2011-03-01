#ifndef ConfigFlags_h
#define ConfigFlags_h

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <exception>


class ConfigFlags ;

inline std::ostream& operator<<(std::ostream& os, const ConfigFlags& cf) ;
  
class ConfigFlags{
  
  friend  std::ostream& operator<<(std::ostream& os, const ConfigFlags& flags) ;

public:

  /** Helper class that holds a fixed number of boolean properties for configuration.
   *  Should be used with const indices defined outside the class (eg static const int in struct).  
   */
  ConfigFlags(unsigned size) : _flags( size ) ,
			       _names( size ) {
  }
  


  ~ConfigFlags(){}
  

  void registerOption( unsigned index, const std::string& name, bool defaultValue=false ){

    if( index > _flags.size() - 1 ) {

      std::stringstream ss ;
      ss << "ConfigFlags - index " << name << " out of range : " << index << " size: " <<  _flags.size() << std::endl ;

      throw std::out_of_range( ss.str() ) ;
    }

    _names[ index ] = name ;
    _flags[ index ] = defaultValue ;
  }
  
  bool option(unsigned index) const { 
    return _flags.at( index ) ; 
  }
  
  bool operator[](unsigned index) const { 
    return _flags.at( index ) ; 
  }
  
  void setOption(unsigned index , bool val) { 
    _flags[ index ] = val ; 
  }
  
  
  std::string& optionName(unsigned index) {
    return _names.at( index ) ;
  }
  
  protected:
  std::vector<bool> _flags ;
  std::vector<std::string> _names ;

};

inline std::ostream& operator<<(std::ostream& os, const ConfigFlags& cf)  {
  
  for( unsigned i=0 ; i < cf._flags.size() ; ++i){

    os << "  option: " << cf._names[i] <<  "\t: " << cf._flags[i] << std::endl ; 
  }
      
  return os ;
} 



#endif
  
