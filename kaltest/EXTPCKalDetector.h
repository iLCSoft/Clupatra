#ifndef __EXTPCDETECTOR__
#define __EXTPCDETECTOR__

#include "EXVKalDetector.h"

class TNode;

namespace gear{
  class TPCParameters ;
}


class EXTPCKalDetector : public EXVKalDetector {
public:
  EXTPCKalDetector(Int_t m = 100);
  
  /** Initialize the TPC from GEAR */
  EXTPCKalDetector(const gear::TPCParameters& tpcParam );
  
  ~EXTPCKalDetector();
  
  static Double_t GetVdrift() { return fgVdrift; }
  
  using EXVKalDetector::Draw;
  void  Draw(Int_t color, const Char_t *opt = "");
  
private:
  TNode *fNodePtr;            // node pointer
   static Double_t fgVdrift;   // drift velocity

   //ClassDef(EXTPCKalDetector,1)   // Sample hit class
};

#endif
