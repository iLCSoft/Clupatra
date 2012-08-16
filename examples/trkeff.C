
void trkeff(const char* FILEN ) {

  float minEff = 0.8 ;
  float maxEff = 1.01 ;

  float minpt = 0.2 ;
  float maxpt = 500. ;

  float minct = 0.0 ;
  float maxct = 1. ;

   const char* title = "TPC track finding efficiency - WW @ 1000 GeV" ;

  //   const char* title = "TPC track finding efficiency - tau pairs @ 500 GeV" ;
  //const char* title = "TPC track finding efficiency - ttbar @ 500 GeV" ;
  //const char* title = "TPC track finding efficiency - tracks from v-zeros (in ttbar @ 500 GeV)" ;

  char* dirClupa = "clupaeff" ;
  //  char* dirTPC   = "marlintrkeff" ;
  char* dirLDC   = "marlintrkeff" ;

  
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  //  gStyle->SetOptFit(111);
  gStyle->SetOptStat(0000);
  //  gStyle->SetOptStat(111111);
  gStyle->SetStatFontSize(.05);
  gStyle->SetStatW(.5);
  gStyle->SetStatH(.15);
  // style->SetStatX(0.90);
  // style->SetStatY(0.90);

  gStyle->SetTitleFont(32,"xyz");
  gStyle->SetTitleSize(0.07,"xyz");
  gStyle->SetTitleOffset(1.1,"yz");
  gStyle->SetTitleOffset(1.0,"x");
  
  // gStyle->SetMarkerStyle( kFullCircle ) ;
  // gStyle->SetMarkerStyle( kOpenCircle ) ;
  gStyle->SetMarkerStyle( kDot ) ;


  gROOT->ForceStyle();


  TFile* f = new TFile( FILEN ) ;


  //---- create some graphs with asymmetric efficiency errors -----------

  f->cd( dirClupa ) ;

  TGraphAsymmErrors* gpt_clupa = new TGraphAsymmErrors ;
  gpt_clupa->Divide( hpt_f , hpt_t , "v" ) ;
  gpt_clupa->SetName("gpt_clupa" ) ;

  // TGraphAsymmErrors* gcosth_clupa = new TGraphAsymmErrors ;
  // gcosth_clupa->Divide( hcosth_f , hcosth_t , "v" ) ;

  TGraphAsymmErrors* gacth_clupa = new TGraphAsymmErrors ;
  gacth_clupa->Divide( hacth_f , hacth_t , "v" ) ;
  gacth_clupa->SetName("gacth_clupa" ) ;


//   f->cd( dirTPC ) ; 
//   TGraphAsymmErrors* gpt_tpc = new TGraphAsymmErrors ;
//   gpt_tpc->Divide( hpt_f , hpt_t , "v" ) ;

//   // TGraphAsymmErrors* gcosth_tpc = new TGraphAsymmErrors ;
//   // gcosth_tpc->Divide( hcosth_f , hcosth_t , "v" ) ;

//   TGraphAsymmErrors* gacth_tpc = new TGraphAsymmErrors ;
//   gacth_tpc->Divide( hacth_f , hacth_t , "v" ) ;
//   gacth_tpc->SetName("gacth_tpc" ) ;



  f->cd( dirLDC ) ; 
  TGraphAsymmErrors* gpt_ldc = new TGraphAsymmErrors ;
  gpt_ldc->Divide( hpt_f , hpt_t , "v" ) ;

  // TGraphAsymmErrors* gcosth_ldc = new TGraphAsymmErrors ;
  // gcosth_ldc->Divide( hcosth_f , hcosth_t , "v" ) ;

  TGraphAsymmErrors* gacth_ldc = new TGraphAsymmErrors ;
  gacth_ldc->Divide( hacth_f , hacth_t , "v" ) ;
  gacth_ldc->SetName("gacth_ldc" ) ;



  // ----------- create histograms that defnie the axis ranges of the plot ----

  hpt_clupa = new TH2F( "hpt_clupa", title ,10, minpt, maxpt, 10, minEff , maxEff  ) ;

  hacth_clupa = new TH2F( "hcosth_clupa", "" , 10, minct, maxct, 10, minEff , maxEff  ) ;
  // hpt_clupa = new TH2F("hpt_clupa","track finding efficiency vs. pt",10, minpt, maxpt, 10, minEff , maxEff  ) ;
  // hacth_clupa = new TH2F("hcosth_clupa","track finding efficiency vs. |cos(theta)|",10, minct, maxct, 10, minEff , maxEff  ) ;

  hpt_clupa->GetYaxis()->SetTitle( "#epsilon_{t_find}" );
  hpt_clupa->GetYaxis()->SetTitleSize( 0.04 ) ;
  hpt_clupa->GetYaxis()->SetTitleOffset( 0.5 );

  hpt_clupa->GetXaxis()->SetTitle( "p_{t}/GeV" );
  hpt_clupa->GetXaxis()->SetTitleSize( 0.04 ) ;
  hpt_clupa->GetXaxis()->SetTitleOffset( 1.2 );

  hacth_clupa->GetYaxis()->SetTitle( "#epsilon_{t_find}" );
  hacth_clupa->GetYaxis()->SetTitleSize( 0.04 ) ;
  hacth_clupa->GetYaxis()->SetTitleOffset( 0.5 );

  hacth_clupa->GetXaxis()->SetTitle( "|cos(#theta)|" );
  hacth_clupa->GetXaxis()->SetTitleSize( 0.05 ) ;
  hacth_clupa->GetXaxis()->SetTitleOffset( 1. );
  //------------------------------
 


  c1 = new TCanvas("c_trkeff", title ,-5);
  c1->Clear();
  c1->Divide(1,2);
  
  int pad = 1;
  
  //-------------------------------------------------------------------------------------
  c1->cd( pad++ ) ;

  gPad->SetLogx() ;

  c1->GetPad(1)->SetTitle("track finding efficiency vs. pt") ;
   
  hpt_clupa->Draw() ;

  gpt_clupa->SetMarkerStyle( kOpenCircle) ;
  gpt_clupa->SetMarkerColor( kRed) ;
  gpt_clupa->SetLineColor( kRed ) ;
  gpt_clupa->Draw("P") ;

//   gpt_tpc->SetMarkerStyle( kOpenTriangleUp ) ;
//   gpt_tpc->SetMarkerColor( kBlue) ;
//   gpt_tpc->SetLineColor( kBlue ) ;
//   gpt_tpc->Draw("P") ;

  gpt_ldc->SetMarkerStyle( kDot) ;
  gpt_ldc->SetMarkerColor( kBlack) ;
  gpt_ldc->SetLineColor( kBlack ) ;
  gpt_ldc->Draw("P") ;


  //========================================================================================
  c1->cd( pad++ ) ;

  c_trkeff->GetPad(2)->SetTitle("track finding efficiency vs. |cos(theta)|") ;

  hacth_clupa->Draw() ;
  gacth_clupa->SetMarkerStyle( kOpenCircle) ;
  gacth_clupa->SetMarkerColor( kRed) ;
  gacth_clupa->SetLineColor( kRed ) ;
  gacth_clupa->Draw("P") ;

//   gacth_tpc->SetMarkerStyle( kOpenTriangleUp ) ;
//   gacth_tpc->SetMarkerColor( kBlue) ;
//   gacth_tpc->SetLineColor( kBlue ) ;
//   gacth_tpc->Draw("P") ;

  gacth_ldc->SetMarkerStyle( kDot) ;
  gacth_ldc->SetMarkerColor( kBlack) ;
  gacth_ldc->SetLineColor( kBlack ) ;
  gacth_ldc->Draw("P") ;


  // TLegend* tl = ((TPad*)c1->GetPad(2))->BuildLegend() ;
  //  tl->Draw()

  TLegend* leg = new TLegend(0.5, 0.2, .8, .5);
  leg->AddEntry( "gacth_clupa" ,  "Clupatra (TPC only)" , "lep") ;
//   leg->AddEntry( "gacth_tpc"   ,  "LEPTracking (TPC only)", "lep") ;
//  leg->AddEntry( "gacth_ldc"   ,  "MarlinTrk (incl VTX/SIT,FTD)", "lep") ;
  leg->Draw();


  // --------------- save to pdf file ----------------------------------

  std::string pdfFile( std::string( FILEN ) + std::string( ".pdf" ) ) ;

  c1->Print( pdfFile.c_str() ) ;
 
}



