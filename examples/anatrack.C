#include <string>
#include <iostream>
#include "TFile.h"
#include "TDirectory.h"


//plot6hists(TFile* fi, char* title, char* h1,char* h2,char* h3,char* h4, char* h5, char* h6=0){
plot6hists(TFile* fi,const char* dir, const std::string& fileName, char* title, char* h1,char* h2,char* h3,char* h4, char* h5, char* h6=0){
 
  c1 = new TCanvas(title,title,-5);
  c1->Clear();
  c1->Divide(3,2);
  
  int pad = 1;
  
  std::string path = std::string("/") + std::string( dir ) + std::string("/" ) ;
  //  std::string path( "/tracks_mup_50gev/" ) ;

  std::string hName ;
  TH1D* hist = 0 ;

  c1->cd(pad++); hName = path + std::string(h1) ;  hist = (TH1D*) fi->Get( hName.c_str() ) ; if(hist)  hist->Draw() ; 
  c1->cd(pad++); hName = path + std::string(h2) ;  hist = (TH1D*) fi->Get( hName.c_str() ) ; if(hist)  hist->Draw() ; 
  c1->cd(pad++); hName = path + std::string(h3) ;  hist = (TH1D*) fi->Get( hName.c_str() ) ; if(hist)  hist->Draw() ; 
  c1->cd(pad++); hName = path + std::string(h4) ;  hist = (TH1D*) fi->Get( hName.c_str() ) ; if(hist)  hist->Draw() ; 
  c1->cd(pad++); hName = path + std::string(h5) ;  hist = (TH1D*) fi->Get( hName.c_str() ) ; if(hist)  hist->Draw() ; 


  if( h6 != 0 ) {

    c1->cd(pad++); hName = path + std::string(h6) ;  hist = (TH1D*) fi->Get( hName.c_str() ) ; if(hist)  hist->Draw() ; 
  }

  // std::string fileName( title ) ;
  // fileName += ".pdf" ;
  c1->Print( fileName.c_str() ) ;

}



void anatrack(const char* FILEN, const char* dirName) {

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  //  gStyle->SetOptFit(111);
  gStyle->SetOptStat(111111);
  gStyle->SetStatFontSize(.05);
  gStyle->SetStatW(.5);
  gStyle->SetStatH(.15);
  // style->SetStatX(0.90);
  // style->SetStatY(0.90);
  gStyle->SetTitleFont(32,"xyz");
  gStyle->SetTitleSize(0.07,"xyz");
  gStyle->SetTitleOffset(1.1,"yz");
  gStyle->SetTitleOffset(1.0,"x");


  
  std::string fN =  std::string(dirName) + ".pdf" ;

  TFile* f = new TFile(FILEN) ; 
  f->cd( dirName  ) ;

  const char* dn = dirName ;
  
  plot6hists( f, dn , fN+"(" ,"track_params",       "hphi",    "hd0" ,    "homega" ,    "hz0",    "htanL"    , "hpt"    ) ;
  plot6hists( f, dn , fN     ,"track_params_mcp",   "hphimcp", "hd0mcp" , "homegamcp" , "hz0mcp", "htanLmcp" , "hptmcp" ) ;
  plot6hists( f, dn , fN     ,"track_params_error", "hephi",   "hed0" ,   "heomega" ,   "hez0",   "hetanL"   , "hept"   ) ;
  plot6hists( f, dn , fN     ,"track_params_delta", "hdphi",   "hdd0" ,   "hdomega" ,   "hdz0",   "hdtanL"   , "hdpt"   ) ;
  plot6hists( f, dn , fN+")" ,"track_params_pulls", "hpphi",   "hpd0" ,   "hpomega" ,   "hpz0",   "hptanL"   , "hppt"   ) ;

}





