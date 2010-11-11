#include <string>
#include <iostream>
#include "TFile.h"
#include "TDirectory.h"

plot6hists(TFile* fi, char* title, char* h1,char* h2,char* h3,char* h4, char* h5, char* h6=0){
 
  c1 = new TCanvas(title,title,-5);
  c1->Clear();
  c1->Divide(3,2);
  
  int pad = 1;
  
  std::string path("MyAnaTrack/") ;
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

  std::string fileName( title ) ;
  fileName += ".pdf" ;
  c1->Print( fileName.c_str() ) ;

}


void anatrack(const char* FILEN) {

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptFit(111);
  gStyle->SetOptStat(111111);
  
  
  TFile* f = new TFile(FILEN) ; 
  f->cd("MyAnaTrack") ;


  plot6hists( f , "track_params",       "hphi",    "hd0" ,    "homega" ,    "hz0",    "htanL"    , "hpt"  ) ;
  plot6hists( f , "track_params_mcp",   "hphimcp", "hd0mcp" , "homegamcp" , "hz0mcp", "htanLmcp"          ) ;
  plot6hists( f , "track_params_error", "hephi",   "hed0" ,   "heomega" ,   "hez0",   "hetanL"   , "hept" ) ;
  plot6hists( f , "track_params_delta", "hdphi",   "hdd0" ,   "hdomega" ,   "hdz0",   "hdtanL"   , "hdpt" ) ;
  plot6hists( f , "track_params_pulls", "hpphi",   "hpd0" ,   "hpomega" ,   "hpz0",   "hptanL"   , "hppt" ) ;

}





