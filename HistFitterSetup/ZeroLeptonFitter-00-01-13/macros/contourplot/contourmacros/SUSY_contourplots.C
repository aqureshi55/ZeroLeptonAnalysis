#include "contourmacros/CombinationGlob.C"
#include "TColor.h"
#include <TStyle.h>
#include "TPave.h"

#include <TFile.h>
#include <TTree.h>
#include <Riostream.h>
#include <TGraph.h>
#include "contourmacros/GG_direct_8TeVobs.C"
#include "contourmacros/SS_direct_8TeVobs.C"
#include "contourmacros/SS_onestepCC_8TeVobs.C"
#include "contourmacros/GG_onestepCCx12_8TeVobs.C"
#include "contourmacros/GG_onestep_mlsp60_8TeVobs.C"
#include "Run1result/SG_pMSSM/SG_pMSSM_0_Run1.C"
#include "Run1result/SG_pMSSM/SG_pMSSM_395_Run1.C"
#include "Run1result/SG_pMSSM/SG_pMSSM_695_Run1.C"
#include "TROOT.h"
#include <algorithm>
#include <iostream>
#include <sstream>
#include "summary_harvest_tree_description.h"

const bool doSmooth=true;
TFile* outfile(0);

void initialize() {
  gSystem->Load("libSusyFitter.so");
}

TH2F* FixAndSetBorders( const TH2& hist, const char* name, const char* title, Double_t val );

TGraph* ContourGraph( TH2F* hist);

TH2F* linearsmooth(const TH2& hist, const char* name, const char* title);

void modifyPointToRemoveLowMLine(double& x, double& y, double frameXmin,  double frameXmax, double frameYmin, TString Grid){
  float xmax = (Grid=="SM_GG_N2") ? 900. : frameXmax/2.;
  //float ymax = (Grid=="SM_GG_N2") ? frameXmin-30. : frameXmin;
  float ymax = frameXmin-30.;
  if( y<ymax && x<xmax && frameXmin<x && frameYmin<y ){
    y = frameYmin-300;
    if( x<=frameXmin )  x=frameXmin-300;
  } 
};

void modifyPoint(double& x, double& y, TString Grid){
  if( Grid=="GG_onestep_mlsp60" ){
    if( y>1. ) y=1. ;
  }else if(Grid=="SG_pMSSM_0"){
    if(x<1800) x=-10000;
    if(y<2000) y=-10000;
    if(y<2500&&x<2200) y=-10000;
  }else if(Grid=="SG_pMSSM_395"){
    ;
  }else if(Grid=="SG_pMSSM_695"){
    if(x<1500) x=-1000;
    if(y<1800) y=-1000;
    if(y<2500&&x<2000) y=-1000;
  }else if(Grid=="SG_pMSSM_995"||Grid=="SG_pMSSM_all"){
    ;
  }else if(Grid!="long_lived_gl"){
    if( x<y ) y=x;
  }
}

  void
DummyLegendExpected(TLegend* leg, TString what,  Int_t fillColor, Int_t fillStyle, Int_t lineColor, Int_t lineStyle, Int_t lineWidth)
{
  TGraph* gr = new TGraph();
  gr->SetFillColor(fillColor);
  gr->SetFillStyle(fillStyle);
  gr->SetLineColor(lineColor);
  gr->SetLineStyle(lineStyle);
  gr->SetLineWidth(lineWidth);
  leg->AddEntry(gr,what,"LF");
}

  void
DummyLegendBestSR(TLegend* leg, TString what,  Int_t fillColor, Int_t fillStyle, Int_t lineColor, Int_t lineStyle, Int_t lineWidth)
{
  TGraph* gr = new TGraph();
  gr->SetFillColor(fillColor);
  gr->SetFillStyle(fillStyle);
  gr->SetLineColor(lineColor);
  gr->SetLineStyle(lineStyle);
  gr->SetLineWidth(lineWidth);
  leg->AddEntry(gr,what,"LF");
}

  void
DummyLegendLine(TLegend* leg, TString what, Int_t lineColor, Int_t lineStyle, Int_t lineWidth)
{
  TGraph* gr = new TGraph();
  gr->SetFillColor(0);
  gr->SetFillStyle(0);
  gr->SetLineColor(lineColor);
  gr->SetLineStyle(lineStyle);
  gr->SetLineWidth(lineWidth);
  leg->AddEntry(gr,what,"L");
}

  TGraph*
ContourGraph( TH2F* hist)
{

  std::cout<<" ContourGraph in mySUSY_contorplots.C : hist address : "<<hist<<std::endl;
  TH2F* h = (TH2F*)hist->Clone();
  TCanvas* ctmp=new TCanvas("ctmp","ctmp",800,600);
  cout << "==> Will dumb histogram: " << h->GetName() << " into a graph" <<endl;
  //h->SetContour( 1 );
  h->SetContour( 2 );
  double pval = CombinationGlob::cl_percent[1];
  double signif = TMath::NormQuantile(1-pval);
  std::cout<<" contour : significance = "<<signif<<std::endl;
  h->SetContourLevel( 0, signif );
  h->SetContourLevel( 1, (signif<5)?5:99 );
  cout << "Draw contour histogram " << h->GetName() <<" : "<<h->GetEntries()<<endl;
  h->Draw("CONT LIST");
  h->SetDirectory(0);
  cout << "Canvas update"<<endl;
  ctmp->Update();
  TObjArray *contours(0);
  contours = (TObjArray*) gROOT->GetListOfSpecials()->FindObject("contours");
  if(!contours){
    std::cout<<"Error!! Failed to get contours from gROOT"<<std::endl;
    gROOT->ls();
    exit(1);
  }
  cout << "Retrieve contours as TGraph"<<endl;
  Int_t ncontours     = contours->GetSize();
  TList *list = (TList*)contours->At(0);
  Int_t number_of_lists = list->GetSize();
  cout << "Retrieve first graph of contour"<<endl;
  TGraph* gr(0);
  gr = (TGraph*)list->At(0)->Clone("gr");
  if(!gr){
    std::cout<<"Error!! Failed to get first graph of contour."<<std::endl;
    exit(1);
  }
  TGraph* grTmp;
  for (int k = 0 ; k<number_of_lists ; k++){
    cout << "Retrieve "<<k<<"th/"<<number_of_lists<<" graph of contour"<<endl;
    grTmp = 0;
    grTmp = (TGraph*)list->At(k)->Clone((TString)("grTmp"+std::to_string(k)) );
    if( !grTmp ){
      std::cout<<"Error!! Failed to get "<<k<<" graph of contour."<<std::endl;
      exit(1);
    }
    Int_t N = gr->GetN();
    Int_t N_tmp = grTmp->GetN();
    if(N < N_tmp){
      std::cout << "Set "<<k<<"th graph as contour."<<std::endl;
      gr = grTmp;
    }
    //
    //    mg->Add((TGraph*)list->At(k));
  }
  std::cout << "Finish retrieving graph of contour"<<std::endl;

  gr->SetName(hist->GetName());
  int N = gr->GetN();
  double x0, y0;

  //  for(int j=0; j<N; j++) {
  //    gr->GetPoint(j,x0,y0);
  //    cout << j << " : " << x0 << " : "<<y0 << endl;
  //  }
  //  //  gr->SetMarkerSize(2.0);
  //  gr->SetMarkerSize(2.0);
  //  gr->SetMarkerStyle(21);

  //  gr->Draw("LP");


  //  cout << "Generated graph " << gr << " with name " << gr->GetName() << endl;
  return gr;

}


TGraph*
DrawExpectedBand( TGraph* gr1,  TGraph* gr2, Int_t fillColor, Int_t fillStyle, Int_t cut = -99999 , TString Grid="", TCanvas* c=0, bool removeLowMLine=false )
{

  //  TGraph* gr1 = new TGraph( *graph1 );
  //  TGraph* gr2 = new TGraph( *graph2 );

  int number_of_bins = max(gr1->GetN(),gr2->GetN());

  const Int_t gr1N = gr1->GetN();
  const Int_t gr2N = gr2->GetN();
  //cout<<" gr1 / gr2 npoint = "<<gr1N<<" / "<<gr2N<<endl;

  const Int_t N = number_of_bins;
  Double_t x1[N], y1[N], x2[N], y2[N];

  Double_t xx0, yy0;
  int npoint=0;

  for(int j=0; j<gr1N; j++) {
    gr1->GetPoint(npoint,xx0,yy0);
    // Modify points forming low mass line
    if(removeLowMLine&&c){
      TH1* hframe=(TH1*)c->FindObject("frame");
      double xmin = hframe->GetXaxis()->GetXmin();
      double xmax = hframe->GetXaxis()->GetXmax();
      double ymin = hframe->GetYaxis()->GetXmin();
      modifyPointToRemoveLowMLine(xx0,yy0,xmin,xmax,ymin,Grid);
    }
    // Modify points above the forbidden line 
    modifyPoint(xx0,yy0,Grid);
    x1[npoint] = xx0;
    y1[npoint] = yy0;
    /*
       if(Grid.BeginsWith("SG_pMSSM_")&&(xx0==-1000||yy0==-1000)){
       if(npoint==0) std::cout<<" Error !! there is not previous point!!"<<std::endl;
       x1[npoint] = x1[npoint-1];
       y1[npoint] = y1[npoint-1];
       }
       */
    cout<<" gr1 i="<<j<<" : (x,y)=("<<x1[npoint]<<","<<y1[npoint]<<")"<<endl;

    npoint++;
  }
  int lastpoint=npoint-1;
  while (npoint < N) {
    x1[npoint] = x1[lastpoint];
    y1[npoint] = y1[lastpoint];
    cout<<" gr1 i="<<npoint<<" : (x,y)=("<<x1[npoint]<<","<<y1[npoint]<<")"<<endl;
    npoint++;
  }

  Double_t xx1, yy1;
  npoint=0;

  for(int j=0; j<gr2N; j++) {
    gr2->GetPoint(npoint,xx1,yy1);
    // Moidfy points forming low mass line
    if(removeLowMLine&&c){
      TH1* hframe=(TH1*)c->FindObject("frame");
      double xmin = hframe->GetXaxis()->GetXmin();
      double xmax = hframe->GetXaxis()->GetXmax();
      double ymin = hframe->GetYaxis()->GetXmin();
      modifyPointToRemoveLowMLine(xx1,yy1,xmin,xmax,ymin,Grid);
    }
    // Modify points above the forbidden line 
    modifyPoint(xx1,yy1,Grid);

    x2[npoint] = xx1;
    y2[npoint] = yy1;
    /*
       if(Grid.BeginsWith("SG_pMSSM_")&&(xx1==-1000||yy1==-1000)){
       if(npoint==0) std::cout<<" Error !! there is not previous point!!"<<std::endl;
       x2[npoint] = x2[npoint-1];
       y2[npoint] = y2[npoint-1];
       }
       */
    cout<<" gr2 i="<<j<<" : (x,y)=("<<x2[npoint]<<","<<y2[npoint]<<")"<<endl;

    npoint++;
  }
  lastpoint=npoint-1;
  while (npoint < N) {
    x2[npoint] = x2[lastpoint];
    y2[npoint] = y2[lastpoint];
    cout<<" gr2 i="<<npoint<<" : (x,y)=("<<x2[npoint]<<","<<y2[npoint]<<")"<<endl;
    npoint++;
  }


  TGraph *grshade = new TGraphAsymmErrors(2*N);

  for (int i=0;i<N;i++) {
    if (x1[i] > cut)
      grshade->SetPoint(i,x1[i],y1[i]);
    if (x2[N-i-1] > cut)
      grshade->SetPoint(N+i,x2[N-i-1],y2[N-i-1]);
  }

  // Apply the cut in the shade plot if there is something that doesn't look good...
  int Nshade = grshade->GetN();
  double x0, y0;
  double x00, y00;

  for(int j=0; j<Nshade; j++) {
    grshade->GetPoint(j,x0,y0);
    if ((x0 != 0) && (y0 != 0)) {
      x00 = x0;
      y00 = y0;
      break;
    }
  }

  for(int j=0; j<Nshade; j++) {
    grshade->GetPoint(j,x0,y0);
    if ((x0 == 0) && (y0 == 0)) 
      grshade->SetPoint(j,x00,y00);
  }

  for(int j=0; j<Nshade; j++) {
    grshade->GetPoint(j,x0,y0);
    cout<<" gShade i="<<j<<" : (x,y)=("<<x0<<","<<y0<<")"<<endl;
  }

  // Now draw the plot... 
  grshade->SetFillStyle(fillStyle);
  grshade->SetFillColor(fillColor);
  grshade->SetLineColor(fillColor);
  grshade->SetMarkerStyle(21);
  /*
     grshade->GetXaxis()->SetRangeUser(200.,2000.);
     grshade->GetXaxis()->SetLimits(200.,2000.);
     grshade->SetMinimum(0.);
     grshade->SetMaximum(1400.);
     */

  ///* to check point
  /*
     grshade->SetMarkerStyle(20);
     grshade->SetMarkerSize(1.);
  //grshade->Draw("FC same");
  */
  grshade->Draw("FL same");

  return grshade;
}


  void
MirrorBorders( TH2& hist )
{
  int numx = hist.GetNbinsX();
  int numy = hist.GetNbinsY();

  Float_t val;
  // corner points
  hist.SetBinContent(0,0,hist.GetBinContent(1,1));
  hist.SetBinContent(numx+1,numy+1,hist.GetBinContent(numx,numy));
  hist.SetBinContent(numx+1,0,hist.GetBinContent(numx,1));
  hist.SetBinContent(0,numy+1,hist.GetBinContent(1,numy));

  for(int i=1; i<=numx; i++){
    hist.SetBinContent(i,0,	   hist.GetBinContent(i,1));
    hist.SetBinContent(i,numy+1, hist.GetBinContent(i,numy));
  }
  for(int i=1; i<=numy; i++) {
    hist.SetBinContent(0,i,      hist.GetBinContent(1,i));
    hist.SetBinContent(numx+1,i, hist.GetBinContent(numx,i));
  }
}


  TH2F*
AddBorders( const TH2& hist, const char* name=0, const char* title=0)
{
  int nbinsx = hist.GetNbinsX();
  int nbinsy = hist.GetNbinsY();

  double xbinwidth = ( hist.GetXaxis()->GetBinCenter(nbinsx) - hist.GetXaxis()->GetBinCenter(1) ) / double(nbinsx-1);
  double ybinwidth = ( hist.GetYaxis()->GetBinCenter(nbinsy) - hist.GetYaxis()->GetBinCenter(1) ) / double(nbinsy-1);

  double xmin = hist.GetXaxis()->GetBinCenter(0) - xbinwidth/2. ;
  double xmax = hist.GetXaxis()->GetBinCenter(nbinsx+1) + xbinwidth/2. ;
  double ymin = hist.GetYaxis()->GetBinCenter(0) - ybinwidth/2. ;
  double ymax = hist.GetYaxis()->GetBinCenter(nbinsy+1) + ybinwidth/2. ;

  TH2F* hist2 = new TH2F(name, title, nbinsx+2, xmin, xmax, nbinsy+2, ymin, ymax);

  for (Int_t ibin1=0; ibin1 <= hist.GetNbinsX()+1; ibin1++) {
    for (Int_t ibin2=0; ibin2 <= hist.GetNbinsY()+1; ibin2++)
      hist2->SetBinContent( ibin1+1, ibin2+1, hist.GetBinContent(ibin1,ibin2) );
  }

  return hist2;
}


void SetBorders( TH2 &hist, Double_t val=0 )
{
  int numx = hist.GetNbinsX();
  int numy = hist.GetNbinsY();

  for(int i=0; i <= numx+1 ; i++){
    hist.SetBinContent(i,0,val);
    hist.SetBinContent(i,numy+1,val);
  }
  for(int i=0; i <= numy+1 ; i++) {
    hist.SetBinContent(0,i,val);
    hist.SetBinContent(numx+1,i,val);
  }
}


  TH2F* 
FixAndSetBorders( const TH2& hist, const char* name=0, const char* title=0, Double_t val=0 )
{
  TH2F* hist0 = (TH2F*) hist.Clone(); // histogram we can modify

  MirrorBorders( *hist0 );    // mirror values of border bins into overflow bins

  TH2F* hist1 = AddBorders( *hist0, "hist1", "hist1" );   
  // add new border of bins around original histogram,
  // ... so 'overflow' bins become normal bins
  SetBorders( *hist1, val );                              
  // set overflow bins to value 1

  TH2F* histX = AddBorders( *hist1, "histX", "histX" );   
  // add new border of bins around original histogram,
  // ... so 'overflow' bins become normal bins

  TH2F* hist3 = (TH2F*)histX->Clone();
  hist3->SetName( name!=0 ? name : "hist3" );
  hist3->SetTitle( title!=0 ? title : "hist3" );

  delete hist0; delete hist1; delete histX;
  return hist3; // this can be used for filled contour histograms
}


  void 
DrawContourSameColor( TLegend *leg, TH2F* hist, Int_t nsigma, TString color, Bool_t second=kFALSE, TH2F* inverse=0, Bool_t linesOnly=kFALSE, Bool_t isnobs=kFALSE )
{
  if (nsigma < 1 || nsigma > 3) {
    cout << "*** Error in CombinationGlob::DrawContour: nsigma out of range: " << nsigma 
      << "==> abort" << endl;
    exit(1);
  }
  nsigma--; // used as array index

  Int_t lcol_sigma;
  Int_t fcol_sigma[3];

  if( color == "pink" ){
    lcol_sigma    = CombinationGlob::c_VDarkPink;
    fcol_sigma[0] = CombinationGlob::c_LightPink;
    fcol_sigma[1] = CombinationGlob::c_LightPink;
    fcol_sigma[2] = CombinationGlob::c_LightPink;
  }
  else if( color == "green" ){ // HF
    lcol_sigma    = CombinationGlob::c_VDarkGreen;
    fcol_sigma[0] = CombinationGlob::c_DarkGreen;
    fcol_sigma[1] = CombinationGlob::c_LightGreen;
    fcol_sigma[2] = CombinationGlob::c_VLightGreen;
  } 
  else if( color == "yellow" ){
    lcol_sigma    = CombinationGlob::c_VDarkYellow;
    fcol_sigma[0] = CombinationGlob::c_DarkYellow;
    fcol_sigma[1] = CombinationGlob::c_DarkYellow;
    fcol_sigma[2] = CombinationGlob::c_White; //c_DarkYellow;
  }
  else if( color == "orange" ){
    lcol_sigma    = CombinationGlob::c_VDarkOrange;
    fcol_sigma[0] = CombinationGlob::c_DarkOrange;
    fcol_sigma[1] = CombinationGlob::c_LightOrange; // c_DarkOrange
    fcol_sigma[2] = CombinationGlob::c_VLightOrange;
  }
  else if( color == "gray" ){
    lcol_sigma    = CombinationGlob::c_VDarkGray;
    fcol_sigma[0] = CombinationGlob::c_LightGray;
    fcol_sigma[1] = CombinationGlob::c_LightGray;
    fcol_sigma[2] = CombinationGlob::c_LightGray;
  }
  else if( color == "blue" ){
    lcol_sigma    = CombinationGlob::c_DarkBlueT1;
    fcol_sigma[0] = CombinationGlob::c_BlueT5;
    fcol_sigma[1] = CombinationGlob::c_BlueT3;
    fcol_sigma[2] = CombinationGlob::c_White;  //CombinationGlob::c_BlueT2;
  }

  // contour plot
  TH2F* h = new TH2F( *hist );
  h->SetContour( 1 );
  double pval = CombinationGlob::cl_percent[1];
  double signif = TMath::NormQuantile(1-pval);
  double dnsigma = double(nsigma)-1.;
  double dsignif = signif + dnsigma;
  h->SetContourLevel( 0, dsignif );

  if( !second ){
    h->SetFillColor( fcol_sigma[nsigma] );

    if (!linesOnly) h->Draw( "samecont0" );
  }

  /*  
  // inverse contours plot, needed in case of contours within contours
  if (inverse!=0) {
  TH2F* g = new TH2F( *inverse );
  g->SetContour( 1 );
  g->SetContourLevel( 0, CombinationGlob::cl_percent[nsigma] );
  if( !second ){
  g->SetFillColor( 0 );
  if (!linesOnly) g->Draw( "samecont0" );
  }
  }
  */

  h->SetLineColor( nsigma==1? 4 : lcol_sigma );
  if (isnobs)h->SetLineColor( nsigma==1? 1 : lcol_sigma );
  //h->SetLineStyle( 4 );
  h->SetLineWidth( 2 );
  h->Draw( "samecont3" );

  if (linesOnly&&!isnobs)
    if(nsigma==1){ leg->AddEntry(h,"expected 95% C.L. exclusion ","l");}
  if (isnobs)
    if(nsigma==1){ leg->AddEntry(h,"observed 95% C.L. exclusion","l");}  
  if (!linesOnly) {
    if(nsigma==0){ leg->AddEntry(h,"expected 68% C.L. exclusion","l"); }
    if(nsigma==2){ leg->AddEntry(h,"expected 99% C.L. exclusion","l");}
  }


}


  void 
DrawContourSameColorDisc( TLegend *leg, TH2F* hist, Double_t nsigma, TString color, Bool_t second=kFALSE, TH2F* inverse=0, Bool_t linesOnly=kFALSE )
{
  if (nsigma < 0.5 || nsigma > 10.5 ) {
    cout << "*** Error in CombinationGlob::DrawContour: nsigma out of range: " << nsigma 
      << "==> abort" << endl;
    exit(1);
  }

  Int_t lcol_sigma;
  Int_t fcol_sigma[3];

  if( color == "pink" ){
    lcol_sigma    = CombinationGlob::c_DarkPink;
    fcol_sigma[0] = CombinationGlob::c_VLightPink;
    fcol_sigma[1] = CombinationGlob::c_VLightPink;
    fcol_sigma[2] = CombinationGlob::c_VLightPink;
  }
  else if( color == "green" ){ // HF
    lcol_sigma    = CombinationGlob::c_VDarkGreen;
    fcol_sigma[0] = CombinationGlob::c_DarkGreen;
    fcol_sigma[1] = CombinationGlob::c_LightGreen;
    fcol_sigma[2] = CombinationGlob::c_VLightGreen;
  } 
  else if( color == "yellow" ){
    lcol_sigma    = CombinationGlob::c_VDarkYellow;
    fcol_sigma[0] = CombinationGlob::c_DarkYellow;
    fcol_sigma[1] = CombinationGlob::c_DarkYellow;
    fcol_sigma[2] = CombinationGlob::c_White; //c_DarkYellow;
  }
  else if( color == "orange" ){
    lcol_sigma    = CombinationGlob::c_VDarkOrange;
    fcol_sigma[0] = CombinationGlob::c_DarkOrange;
    fcol_sigma[1] = CombinationGlob::c_LightOrange; // c_DarkOrange
    fcol_sigma[2] = CombinationGlob::c_VLightOrange;
  }
  else if( color == "gray" ){
    lcol_sigma    = CombinationGlob::c_VDarkGray;
    fcol_sigma[0] = CombinationGlob::c_LightGray;
    fcol_sigma[1] = CombinationGlob::c_LightGray;
    fcol_sigma[2] = CombinationGlob::c_LightGray;
  }
  else if( color == "blue" ){
    lcol_sigma    = CombinationGlob::c_DarkBlueT1;
    fcol_sigma[0] = CombinationGlob::c_LightBlue;
    fcol_sigma[1] = CombinationGlob::c_LightBlue;
    fcol_sigma[2] = CombinationGlob::c_LightBlue;
  }

  // contour plot
  TH2F* h = new TH2F( *hist );
  h->SetContour( 1 );
  double dsignif = double (nsigma);
  h->SetContourLevel( 0, dsignif );

  //Int_t mycolor = (nsigma==3   ? 0 : 2);
  Int_t mycolor = (nsigma==2 ? 1 : 2);

  h->SetFillStyle(3003);

  if( !second ){
    h->SetFillColor( fcol_sigma[mycolor] );
    if (!linesOnly) h->Draw( "samecont0" );
  }

  h->SetLineColor( (nsigma==3) ? lcol_sigma : lcol_sigma );

  h->SetLineStyle( nsigma==3 || nsigma==2 ? 1 : 2 );
  h->SetLineWidth( nsigma==3 || nsigma==2 ? 2 : 1 );


  h->Draw( "samecont3" );

  if(nsigma==3)   { leg->AddEntry(h,"3 #sigma evidence","l"); }
  if(nsigma==6)   { leg->AddEntry(h,"N (int) #sigma","l"); }
  if(nsigma==2)   { leg->AddEntry(h,"2 #sigma evidence","l"); }
}




  void
DrawContourMassLine(TH2F* hist, Double_t mass, int color=14 )
{

  // contour plot
  TH2F* h = new TH2F( *hist );

  //  Double_t contours[5] = {500, 1000, 1500, 2000, 2500}
  h->SetContour( 1 );
  //h->SetContour( 5, contours )
  //  h->SetContourLevel( 0, contours );
  h->SetContourLevel( 0, mass );

  h->SetLineColor( color );
  h->SetLineStyle( 7 );
  h->SetLineWidth( 1 );
  h->Draw( "samecont3" );

}


void 
DrawContourLine95( TLegend *leg, TH2F* hist, const TString& text="", Int_t linecolor=CombinationGlob::c_VDarkGray, Int_t linestyle=2, Int_t linewidth=2 , TString Grid="", TCanvas* c=0, bool removeLowMLine=false, bool save=true )
{

  // contour plot
  TH2F* h = new TH2F( *hist );
  h->SetContour( 1 );
  double pval = CombinationGlob::cl_percent[1];
  double signif = TMath::NormQuantile(1-pval);
  //cout << "signif: " <<signif << endl;
  h->SetContourLevel( 0, signif );

  h->SetLineColor( linecolor );
  h->SetLineWidth( linewidth );
  h->SetLineStyle( linestyle );
  //h->Draw( "samecont3" );

  //if (!text.IsNull()) leg->AddEntry(h,text.Data(),"l");
  if (text!="") leg->AddEntry(h,text.Data(),"l");

  // Get TGraph from TH2
  TString histname = h->GetName();
  TGraph* g = (TGraph*) ContourGraph( h )->Clone("g"+histname);
  g->SetLineColor( linecolor );
  g->SetLineWidth( linewidth );
  g->SetLineStyle( linestyle );
  // Modify points
  double x=0.,y=0.;
  for( int n=0 ; n<g->GetN() ; n++ ){
    x=0., y=0.;
    g->GetPoint(n,x,y);
    cout<<" graph : "<<histname<<" (x,y)=("<<x<<","<<y<<")"<<endl;
    // Modify points forming low mass line
    if(removeLowMLine){
      TH1* hframe=(TH1*)c->FindObject("frame");
      double xmin = hframe->GetXaxis()->GetXmin();
      double xmax = hframe->GetXaxis()->GetXmax();
      double ymin = hframe->GetYaxis()->GetXmin();
      modifyPointToRemoveLowMLine(x,y,xmin,xmax,ymin,Grid);
    }
    // Modify points above the forbidden
    modifyPoint(x,y, Grid);

    g->SetPoint(n,x,y);
    cout<<"          -->> set point (x,y)=("<<x<<","<<y<<")"<<endl;
  }

  // Draw
  if(c) c->cd();
  g->Draw("same");

  // Save TH2&TGraph to root file
  if(save){
    TString outfileName = Grid+"_contours.root" ;
    if(!outfile) outfile = new TFile(outfileName, "recreate");
    if( histname!="contour_esigxspls"&&histname!="contour_esigxsmls" ){
      outfile->cd();
      g->Write();
      h->Write();
    }
  }

}

void 
DrawContourLine95( TLegend *leg, TGraph* graph, const TString& text="", Int_t linecolor=CombinationGlob::c_VDarkGray, Int_t linestyle=2, Int_t linewidth=2 , TString Grid="", TCanvas* c=0, bool removeLowMLine=false, bool save=true )
{

  graph->SetLineColor( linecolor );
  graph->SetLineWidth( linewidth );
  graph->SetLineStyle( linestyle );
  // Modify points
  double x=0.,y=0.;
  for( int n=0 ; n<graph->GetN() ; n++ ){
    x=0., y=0.;
    graph->GetPoint(n,x,y);
    cout<<" graph : "<<graph->GetName()<<" (x,y)=("<<x<<","<<y<<")"<<endl;
    // Modify points forming low mass line
    if(removeLowMLine){
      TH1* hframe=(TH1*)c->FindObject("frame");
      double xmin = hframe->GetXaxis()->GetXmin();
      double xmax = hframe->GetXaxis()->GetXmax();
      double ymin = hframe->GetYaxis()->GetXmin();
      modifyPointToRemoveLowMLine(x,y,xmin,xmax,ymin,Grid);
    }
    // Modify points above the forbidden
    modifyPoint(x,y,Grid);

    graph->SetPoint(n,x,y);
    cout<<"          -->> set point (x,y)=("<<x<<","<<y<<")"<<endl;
  }

  // Draw
  if(c) c->cd();
  graph->Draw( "samel" );

  //if (!text.IsNull()) leg->AddEntry(h,text.Data(),"l");
  if (text!="") leg->AddEntry(graph,text.Data(),"l");

  TString outfileName = Grid+"_contours.root" ;
  TString name = graph->GetName();
  if(save){
    if(!outfile) outfile = new TFile(outfileName, "recreate");
    outfile->cd();
    graph->Write();
  }

}



TGraphAsymmErrors*
GetContourLine95( TH2F* hist, TString Grid="" )
{
  //cout<<"input contour : "<<hist->GetName()<<endl;
  // contour plot
  TH2F* h = new TH2F( *hist );
  h->SetContour( 1 );
  double pval = CombinationGlob::cl_percent[1];
  double signif = TMath::NormQuantile(1-pval);
  //cout << "signif: " <<signif << endl;
  h->SetContourLevel( 0, signif );
  TString histname=h->GetName();

  TGraphAsymmErrors* g = (TGraphAsymmErrors*) ContourGraph( h )->Clone("g"+histname);
  return g;
}


void
DrawContourLine3sigma( TLegend *leg, TH2F* hist, const TString& text="", Int_t linecolor=CombinationGlob::c_VDarkGray, Int_t linestyle=2, Int_t linewidth=2 )
{
  // contour plot
  TH2F* h = new TH2F( *hist );
  h->SetContour( 1 );
  double pval = (1.-0.9973)*0.5; // one-sided
  double signif = TMath::NormQuantile(1-pval);
  cout <<"DrawContourLine3sigma: pval="<<pval<<", "<<signif<<"sigma for "<<text<<endl;
  h->SetContourLevel( 0, signif );

  h->SetLineColor( linecolor );
  h->SetLineWidth( linewidth );
  h->SetLineStyle( linestyle );
  h->Draw( "samecont3" );

  if (!text.IsNull()) leg->AddEntry(h,text.Data(),"l");
}

void
DrawContourLine2sigma( TLegend *leg, TH2F* hist, const TString& text="", Int_t linecolor=CombinationGlob::c_VDarkGray, Int_t linestyle=2, Int_t linewidth=2 )
{
  // contour plot
  TH2F* h = new TH2F( *hist );
  h->SetContour( 1 );
  double pval = (1.-0.9545)*0.5; // one-sided
  double signif = TMath::NormQuantile(1-pval);
  cout <<"DrawContourLine3sigma: pval="<<pval<<", "<<signif<<"sigma for "<<text<<endl;
  h->SetContourLevel( 0, signif );

  h->SetLineColor( linecolor );
  h->SetLineWidth( linewidth );
  h->SetLineStyle( linestyle );
  h->Draw( "samecont3" );

  if (!text.IsNull()) leg->AddEntry(h,text.Data(),"l");
}

void
DrawContourLine99( TLegend *leg, TH2F* hist, const TString& text="", Int_t linecolor=CombinationGlob::c_VDarkGray, Int_t linestyle=2, Int_t linewidth=2 )
{
  // contour plot
  TH2F* h = new TH2F( *hist );
  h->SetContour( 1 );
  double pval = CombinationGlob::cl_percent[2];
  double signif = TMath::NormQuantile(1-pval);

  h->SetContourLevel( 0, signif );

  h->SetLineColor( linecolor );
  h->SetLineWidth( linewidth );
  h->SetLineStyle( linestyle );
  h->Draw( "samecont3" );

  if (!text.IsNull()) leg->AddEntry(h,text.Data(),"l");
}


void
DrawContourLine68( TLegend *leg, TH2F* hist, const TString& text="", Int_t linecolor=CombinationGlob::c_VDarkGray, Int_t linestyle=2 )
{
  // contour plot
  TH2F* h = new TH2F( *hist );
  h->SetContour( 1 );
  double pval = CombinationGlob::cl_percent[0];
  double signif = TMath::NormQuantile(1-pval);

  h->SetContourLevel( 0, signif );

  h->SetLineColor( linecolor );
  h->SetLineWidth( 2 );
  h->SetLineStyle( linestyle );
  h->Draw( "samecont3" );

  if (!text.IsNull()) leg->AddEntry(h,text.Data(),"l");
}

TH2F* linearsmooth(const TH2& hist, const char* name, const char* title) {
  int nbinsx = hist.GetNbinsX();
  int nbinsy = hist.GetNbinsY();

  double xbinwidth = ( hist.GetXaxis()->GetBinCenter(nbinsx) - hist.GetXaxis()->GetBinCenter(1) ) / double(nbinsx-1);
  double ybinwidth = ( hist.GetYaxis()->GetBinCenter(nbinsy) - hist.GetYaxis()->GetBinCenter(1) ) / double(nbinsy-1);

  int nbinsxsm = 2*nbinsx - 1 ;
  int nbinsysm = 2*nbinsy - 1 ;

  double xmin = hist.GetXaxis()->GetBinCenter(1) - xbinwidth/4. ;
  double xmax = hist.GetXaxis()->GetBinCenter(nbinsx) + xbinwidth/4. ;
  double ymin = hist.GetYaxis()->GetBinCenter(1) - ybinwidth/4. ;
  double ymax = hist.GetYaxis()->GetBinCenter(nbinsy) + ybinwidth/4. ;

  TH2F* hist2 = new TH2F(name, title, nbinsxsm, xmin, xmax, nbinsysm, ymin, ymax);

  for (Int_t ibin1=1; ibin1 < hist.GetNbinsX(); ibin1++) {
    for (Int_t ibin2=1; ibin2 < hist.GetNbinsY(); ibin2++) {
      float f00 = hist.GetBinContent(ibin1,ibin2);
      float f10 = hist.GetBinContent(ibin1+1,ibin2);
      float f01 = hist.GetBinContent(ibin1,ibin2+1);
      float f11 = hist.GetBinContent(ibin1+1,ibin2+1);

      for (Int_t i=0; i<=2; ++i)
        for (Int_t j=0; j<=2; ++j) {
          float x = i*0.5; float y = j*0.5;
          float val = (1-x)*(1-y)*f00 + x*(1-y)*f10 + (1-x)*y*f01 + x*y*f11 ;
          Int_t jbin1 = 2*ibin1 - 1 + i;
          Int_t jbin2 = 2*ibin2 - 1 + j;
          hist2->SetBinContent(jbin1,jbin2,val);
        }
    }
  }

  return hist2; // caller owns histogram
}

void GetSRName(TString txtfile, std::map<int,TString>& SRNames){

  TString infoFileName = txtfile + "_infoFile";
  std::ifstream ifs(infoFileName);
  if( !ifs ){
    std::cout<<" File open error!! : "<<infoFileName<<std::endl;
    std::cout<<" Could not get SR name!!"<<std::endl;
    return;
  }
  while(!ifs.eof()){
    TString colon, SRName;
    int  index;
    ifs>>index>>colon>>SRName;
    SRNames.insert(std::make_pair(index,SRName));
  }
  ifs.close();

  return;
}


TString GetSRName(int fID, TString txtfile){

  TString text = "non" ;
  std::map<int, TString> SRNames;
  GetSRName(txtfile, SRNames);
  std::map<int, TString>::iterator find = SRNames.find(fID);
  std::map<int, TString>::iterator end  = SRNames.end();
  if( find!=end ){
    text = (*find).second;
  }else{
    std::cout<<" Failded to find fID="<<fID<<" SR name."<<std::endl;
  }
  return text;
}


void Show_SR(TString oredList,  TCanvas *can, TH1* frame, float xlow, float xhigh, float ylow, float yhigh, TString Grid, TLegend* leg)
{
  std::stringstream ss; TString str;
  Int_t c_myRed      = TColor::GetColor("#aa000");

  can->cd();

  TLatex lat;
  //lat.SetTextAlign( 11 );
  lat.SetTextColor( 12 );
  lat.SetTextFont( 42 );

  cout << "Draw signal region labels." << endl;
  //gROOT->ProcessLine(".L summary_harvest_tree_description.h+");
  gSystem->Load("libSusyFitter.so");

  TString txtfile=oredList;
  txtfile.ReplaceAll(".root","");

  TTree* tree = harvesttree( txtfile!=0 ? txtfile : "" );
  if (tree==0) { 
    cout << "Cannot open list file. Exit." << endl;
    return;
  }

  Char_t fID[3];
  Float_t m0; 
  Float_t m12; 

  TBranch *b_m0;
  TBranch *b_m12;
  TBranch *b_fID;

  TString m_x = "m0";
  TString m_y = "m12";
  if(Grid=="GG_onestepCC"){
    m_x = "mgluino";
    m_y = "mlsp";
  }else if(Grid=="SS_onestepCC"){
    m_x = "msquark";
    m_y = "mlsp";
  }else if (Grid=="SS_direct" || Grid=="GG_direct"){
    m_x = "m0";
    m_y = "m12";
  }else if (Grid=="SM_GG_N2"){
    m_x = "mgluino";
    m_y = "mlsp2";
  }else if(Grid=="SS_onestepCC_mlsp60"){
    m_x = "msquark";
    m_y = "mchargino";
  }else if(Grid=="GG_onestepCC_mlsp60"){
    m_x = "mgluino";
    m_y = "mchargino";
  }else if(Grid=="SS_onestepN2C_mlsp60"){
    m_x = "msquark";
    m_y = "mneutralino2";
  }else if(Grid=="GG_onestepN2C_mlsp60"){
    m_x = "mgluino";
    m_y = "mneutralino2";
  }else if(Grid=="long_lived_gl"){
    m_x = "lifetime";
    m_y = "mgluino";
  }else if (Grid.BeginsWith("SG_pMSSM")){
    m_x = "mgluino";
    m_y = "msquark";
  }

  tree->SetBranchAddress(m_x, &m0, &b_m0);
  tree->SetBranchAddress(m_y, &m12, &b_m12);
  tree->SetBranchAddress("fID",  &(fID[0]),  &b_fID);

  for( Int_t i = 0; i < tree->GetEntries(); i++ ){

    tree->GetEntry( i );

    if(Grid=="long_lived_gl")    m0  = std::log10(m0);

    int   _m0 = (int)   m0;
    float _m12= (float) m12;
    // 1% length of canvas frame width
    float dX = (xhigh-xlow)/100;
    float dY = (yhigh-ylow)/100;

    int fid = std::atoi((TString)fID);
    TString mySR = GetSRName(fid,txtfile);
    mySR.ReplaceAll("SRJigsawSR","");
    mySR.ReplaceAll("SR","");
    cout <<"BestSR for "<< m0 << " " << m12 << " is: " << fID <<" "<<mySR<<" xmax="<<frame->GetXaxis()->GetXmax()<< endl;
    if( Grid.EndsWith("_mlsp60") ){
      float y = (_m12-60.) / (_m0-60.) ;
      lat.SetTextSize( 0.02 );
      if( (_m0 >= xlow) && (y <= yhigh) && (m0<=frame->GetXaxis()->GetXmax()) ){
        float xlat=m0-2.5*dX;
        float ylat=y+dY;
        if( m0-dX*2.5 < xlow+dX  ) xlat=xlow+dX;
        if( ylow-0.2<y && y<ylow ) ylat=ylow+dY;
        lat.DrawLatex(xlat, ylat, mySR.Data());
        ss.str(""); ss<<"("<<m0<<","<<y<<")"; str=ss.str();
      }
    }else{
      if( (_m0 >= xlow) && (_m12 >= ylow)&& (_m12 <= yhigh) && (m0<=frame->GetXaxis()->GetXmax()) ){

        // ignore points overlappting with legend
        if(leg){
          float x1 = leg->GetX1()*(xhigh-xlow)+xlow;
          float y1 = leg->GetY1()*(yhigh-ylow)+ylow;
          std::cout<<" legend x1="<<x1<<", y1="<<y1<<" dX="<<dX<<" dY="<<dY<<std::endl;
          float marginX=5*dX;
          float marginY=0.;
          if(Grid.BeginsWith("SG_pMSSM_")) marginY = 20*dX;
          if(_m0>x1-marginX&&_m12>y1-marginY) continue;
        }

        lat.SetTextSize( 0.02 );
        if     ( m0==frame->GetXaxis()->GetXmax() ) lat.DrawLatex(m0-dX*6. , m12+dY,   mySR.Data());
        else if( m12==yhigh                       ) lat.DrawLatex(m0-dX*2.5, m12-dY*3, mySR.Data());
        else                                        lat.DrawLatex(m0-dX*2.5, m12+dY,   mySR.Data());
        //ss.str(""); ss<<"("<<m0<<","<<m12<<")"; str=ss.str();
        //lat.DrawLatex(m0, m12, str.Data());
        //lat.DrawLatex(m0, m12, "+");
      }

      //}
  }
}
}

TGraphAsymmErrors* GetLimitGraph(TString filename, TString graphName="observed"){
  TFile* file(0);
  file = new TFile(filename);
  if( !file || !file->IsOpen() ){
    cout<<" file open error!! : "<<filename<<endl;
    return 0;
  }
  cout<<" get "<<graphName<<" from "<<filename<<endl;
  TGraphAsymmErrors* g(0);
  g = (TGraphAsymmErrors*)file->Get(graphName);
  if( !g ){
    cout<<" failed to get "<<graphName<<" from "<<filename<<endl;
    file->ls();
  }
  return (TGraphAsymmErrors*)g;
};



void SUSY_contourplots(
    TString fname0="", TString fname1="", TString fname2="", 
    const char* prefix="",
    TString lumi = 5,
    TString Grid="GG_direct",
    TString outdir="plots",
    bool showSR=false,
    bool unblind=false,
    bool showAllSRExp=false,
    bool drawOldLimits=true,
    bool useTGraph=false,
    int discexcl = 1,
    TString SRLabel = "",
    bool removeLowMLine=false)
{

  initialize();
  cout<<"Grid = "<<Grid<<endl;

  if(useTGraph) showAllSRExp=false;

  //set style and remove existing canvas'
  CombinationGlob::Initialize();

  TString yMass="mLSP";
  TString yParticleTex="m_{#tilde{#chi}_{1}^{0}}";
  if( Grid=="SM_GG_N2" ){
    yMass="m2ndLSP";
    yParticleTex="m_{#tilde{#chi}_{2}^{0}}";
  }
  if( Grid=="GG_onestepCC_mlsp60" ){
    yMass="x";
    yParticleTex="x=#Delta m(#tilde{#chi}_{1}^{#pm}, #tilde{#chi}_{1}^{0})/#Delta m(#tilde{g}, #tilde{#chi}_{1}^{0})";
  }
  if( Grid=="SS_onestepCC_mlsp60" ){
    yMass="x";
    yParticleTex="x=#Delta m(#tilde{#chi}_{1}^{#pm}, #tilde{#chi}_{1}^{0})/#Delta m(#tilde{g}, #tilde{#chi}_{1}^{0})";
  }
  if( Grid=="GG_onestepN2C_mlsp60" ){
    yMass="x";
    yParticleTex="x=#Delta m(#tilde{#chi}_{2}^{0}, #tilde{#chi}_{1}^{0})/#Delta m(#tilde{g}, #tilde{#chi}_{1}^{0})";
  }
  if( Grid=="SS_onestepN2C_mlsp60" ){
    yMass="x";
    yParticleTex="x=#Delta m(#tilde{#chi}_{2}^{0}, #tilde{#chi}_{1}^{0})/#Delta m(#tilde{q}, #tilde{#chi}_{1}^{0})";
  }
  if( Grid=="long_lived_gl" ){
    yMass="mgluino";
    yParticleTex="m_{#tilde{g}}";
  }
  if( Grid.BeginsWith("SG_pMSSM") ){
    yMass="msquark";
    yParticleTex="m_{#tilde{q}}";
  }

  // for SG_pMSSM_all
  std::vector<TGraph*> gobscls0all;
  std::vector<TGraph*> gobscls1all;
  std::vector<TGraph*> gobscls2all;
  std::vector<TGraph*> g0_expclsall;
  std::vector<TGraph*> g0_exp1suall;
  std::vector<TGraph*> g0_exp1sdall;
  std::vector<TGraphAsymmErrors*> gRun1all;

  int mlsps[3]={0,695,995};
  if(Grid=="SG_pMSSM_all"){
    // get TGraph from output root file of each SG_pMSSM_xxx

    TFile * file0 = new TFile("SG_pMSSM_0_contours.root");
    TFile * file1 = new TFile("SG_pMSSM_695_contours.root");
    TFile * file2 = new TFile("SG_pMSSM_995_contours.root");
    std::vector<TFile*> files;
    files.push_back(file0);
    files.push_back(file1);
    files.push_back(file2);

    for( int i=0; i<files.size() ; i++ ){
      gobscls0all .push_back( (TGraph*)files[i]->Get("gobserved"  ) );
      gobscls1all .push_back( (TGraph*)files[i]->Get("gobserved_1"  ) );
      gobscls2all .push_back( (TGraph*)files[i]->Get("gobserved_2"  ) );
      g0_expclsall.push_back( (TGraph*)files[i]->Get("gexpected"    ) );
      g0_exp1suall.push_back( (TGraph*)files[i]->Get("gexpected_1su") );
      g0_exp1sdall.push_back( (TGraph*)files[i]->Get("gexpected_1sd") );
      gRun1all    .push_back( (TGraphAsymmErrors*)files[i]->Get("gRun1") );
    }

    fname0 = fname0.ReplaceAll("all","0");
    fname1 = fname1.ReplaceAll("all","0");
    fname2 = fname2.ReplaceAll("all","0");
    cout << "file name0 = "<<fname0 << endl;
    cout << "file name1 = "<<fname1 << endl;
    cout << "file name2 = "<<fname2 << endl;
  } // end of SG_pMSSM_all setting



  cout << "--- Plotting mGluino versus "+yMass+" " << endl;

  // Get all SR names
  std::map<int, TString> SRNamesMap;
  TString harvestFileName=fname0;
  harvestFileName.ReplaceAll(".root","");
  GetSRName(harvestFileName,SRNamesMap);
  std::vector<TString> SRNames;
  for( auto itr = SRNamesMap.begin() ; itr!=SRNamesMap.end() ; ++itr ){
    SRNames.push_back((*itr).second);
  }

  // Open files and retrieve histogram
  TFile* f0;
  TFile* f1;
  TFile* f2;

  TString dirname=lumi+"ifb";

  // nominal file
  cout << "--- Reading root base file0: " << fname0 << endl;
  f0 = TFile::Open( fname0, "READ" );
  //f0->ls();
  // up file
  cout << "--- Reading root base file1: " << fname1 << endl;
  f1 = TFile::Open( fname1, "READ" );
  //f1->ls();
  // down file
  cout << "--- Reading root base file2: " << fname2 << endl;
  f2 = TFile::Open( fname2, "READ" );
  //f2->ls();

  if(!f0 || !f1 || !f2) cout << "Cannot open!"<< endl;

  TH2F* histobscls0(0);
  TH2F* histobscls1(0);
  TH2F* histobscls2(0);
  TH2F* hist0_expcls(0);
  TH2F* hist0_exp1su(0);
  TH2F* hist0_exp1sd(0);
  TGraph* gobscls0(0);
  TGraph* gobscls1(0);
  TGraph* gobscls2(0);
  TGraph* g0_expcls(0);
  TGraph* g0_exp1su(0);
  TGraph* g0_exp1sd(0);

  TString excl_p1exp;
  TString excl_p1obs;

  bool draw_obscls_error=true;

  if(!useTGraph){
    // Which variable to plot?
    excl_p1exp = "sigp1expclsf";
    excl_p1obs = "sigp1clsf";

    // hists for observed xsec 0/+/-1 lines
    cout << "--- Get observed limits " << endl;
    // nominal
    histobscls0 = (TH2F*)f0->Get(excl_p1obs);
    // +/- 1 sigma
    if(f1) histobscls1 = (TH2F*)f1->Get(excl_p1obs);  
    if(f2) histobscls2 = (TH2F*)f2->Get(excl_p1obs);
    cout << " Draw observed +/- 1sigma liens? : " << draw_obscls_error << endl;
    //TH2F* histe_esigxsp1s;                              
    //TH2F* histe_esigxsm1s;      
    //if(f1) histe_esigxsp1s = (TH2F*)f1->Get( excl_p1obs ); 
    //if(f2) histe_esigxsm1s = (TH2F*)f2->Get( excl_p1obs );

    // hists for expected xsec 0/+/-1 lines
    cout << "--- Get expected limits " << endl;
    // nominal
    hist0_expcls = (TH2F*)f0->Get(excl_p1exp);
    // +/- 1 sigma (from nominal file)
    hist0_exp1su = (TH2F*)f0->Get( "sigclsu1s" );
    if(!hist0_exp1su){
      std::cout<<"Error!! Failed to get sigclsu1s from "<<f0->GetName()<<std::endl;
    }
    hist0_exp1sd = (TH2F*)f0->Get( "sigclsd1s" );
    if(!hist0_exp1sd){
      std::cout<<"Error!! Failed to get sigclsd1s from "<<f0->GetName()<<std::endl;
    }
    draw_obscls_error=histobscls1&&histobscls2;
  }else{
    // Which variable to plot?
    TString excl_p1exp = "CLsexp_Contour_0";
    TString excl_p1obs = "CLs_Contour_0";

    TString excl_p1expu1s = "clsu1s_Contour_0";
    TString excl_p1expd1s = "clsd1s_Contour_0";

    // hists for observed xsec 0/+/-1 lines
    cout << "--- Get observed limits " << endl;
    // nominal
    gobscls0 = (TGraph*)f0->Get(excl_p1obs);
    // +/- 1 sigma
    if(f1) gobscls1 = (TGraph*)f1->Get(excl_p1obs);  
    if(f2) gobscls2 = (TGraph*)f2->Get(excl_p1obs);
    draw_obscls_error=gobscls1&&gobscls2;

    // hists for expected xsec 0/+/-1 lines
    cout << "--- Get expected limits " << endl;
    // nominal
    g0_expcls = (TGraph*)f0->Get(excl_p1exp);
    // +/- 1 sigma (from nominal file)
    g0_exp1su = (TGraph*)f0->Get(excl_p1expu1s);
    if(!g0_exp1su){
      std::cout<<"Error!! Failed to get "<<excl_p1expu1s<<" from "<<f0->GetName()<<std::endl;
      //exit(1);
    }
    g0_exp1sd = (TGraph*)f0->Get(excl_p1expd1s);
    if(!g0_exp1sd){
      std::cout<<"Error!! Failed to get "<<excl_p1expd1s<<" from "<<f0->GetName()<<std::endl;
      //exit(1);
    }
  }

  if(doSmooth&&!useTGraph){
    histobscls0->Smooth();
    if(draw_obscls_error){
      histobscls1->Smooth();
      histobscls2->Smooth();
    }
    //histe_esigxsp1s->Smooth();
    //histe_esigxsm1s->Smooth();
    hist0_expcls->Smooth();
    hist0_exp1su->Smooth();
    hist0_exp1sd->Smooth();
  }

  TH2F* contour_exp(0);
  TH2F* contour_exp1su(0);
  TH2F* contour_exp1sd(0);
  TGraph* gr_contour_exp1su(0);
  TGraph* gr_contour_exp1sd(0);
  TH2F* contour_obscls_1(0);
  TH2F* contour_obscls_2(0);
  TH2F* contour_obscls_0(0);
  if(!useTGraph){
    // get observed exclusion TH2 : nominal
    if (hist0_expcls!=0) { contour_exp     = FixAndSetBorders( *hist0_expcls, "expected", "contour_exp", 0 ); }

    // get expected exclusion TH2 : +/- 1 sigma
    cout << "--- FixAndSetBorders : " << hist0_exp1su->GetName() << endl;
    if (hist0_exp1su!=0) { contour_exp1su     = FixAndSetBorders( *hist0_exp1su, "expected_1su", "expected limit 95\% CL+1 sigma", 0 ); }
    cout << "--- FixAndSetBorders : " << hist0_exp1sd->GetName() << endl;
    if (hist0_exp1sd!=0) { contour_exp1sd     = FixAndSetBorders( *hist0_exp1sd, "expected_1sd", "Expected limit 95\% CL-1 sigma", 0 ); }     
    // get expected exclusion TGraph : +/- 1 sigma
    cout << "--- Retrieveing TGraph from TH2D : " << hist0_exp1su->GetName() << endl;
    //TGraph* gr_contour_exp1su = (TGraph*) ContourGraph( hist0_exp1su )->Clone("expected_1su");
    gr_contour_exp1su = (TGraph*) ContourGraph( contour_exp1su )->Clone("expected_1su");
    gr_contour_exp1su->SetTitle("Exptected limit 95\% CL+1 sigma");
    cout << "--- Retrieveing TGraph from TH2D : " << hist0_exp1sd->GetName() << endl;
    //TGraph* gr_contour_exp1sd = (TGraph*) ContourGraph( hist0_exp1sd )->Clone("expected_1sd");
    gr_contour_exp1sd = (TGraph*) ContourGraph( contour_exp1sd )->Clone("expected_1sd");
    gr_contour_exp1sd->SetTitle("Expected limit 95\% CL-1 sigma");

    // get observed exclusion TH2 : +/- 1 sigma
    /*
       cout << "--- FixAndSetBorders observed band : " << histe_esigxsp1s->GetName()<< endl;
       TH2F* contour_esigxsp1s
       = ( histe_esigxsp1s!=0 ? FixAndSetBorders( *histe_esigxsp1s, "observed_1su", "Observed limit 95\% CL+1 sigma", 0 ) : 0);                              
       cout << "--- FixAndSetBorders expected band : " << histe_esigxsm1s->GetName()<< endl;
       TH2F* contour_esigxsm1s                              
       = ( histe_esigxsm1s!=0 ? FixAndSetBorders( *histe_esigxsm1s, "observed_1sd", "Observed limit 95\% CL-1 sigma", 0 ) : 0);      
       */

    // set observed exclusion TH2 : +/- 1 sigma
    if (draw_obscls_error) { 
      cout << "--- FixAndSetBorders observed band : " << histobscls1->GetName()<< endl;
      contour_obscls_1     = FixAndSetBorders( *histobscls1, "observed_1", "contour_obscls_1", 0 );
      cout << "--- FixAndSetBorders observed band : " << histobscls2->GetName()<< endl;
      contour_obscls_2     = FixAndSetBorders( *histobscls2, "observed_2", "contour_obscls_2", 0 );
    }
    cout << "--- FixAndSetBordersobserved band : " << histobscls0->GetName()<< endl;
    histobscls0->Print();
    if (histobscls0!=0) { contour_obscls_0     = FixAndSetBorders( *histobscls0, "observed", "Observed limit 95\% CL", 0 ); }   

  }else{
    gr_contour_exp1su = g0_exp1su;
    gr_contour_exp1sd = g0_exp1sd;
  }


  gStyle->SetPaintTextFormat(".2g");

  // --- draw
  cout << "--- Drawing " << endl;
  // Create canvas
  TCanvas* c = new TCanvas( "c", "A scan of m_#tilde{g} versus "+yParticleTex, 0, 0, 
      CombinationGlob::StandardCanvas[0], CombinationGlob::StandardCanvas[1] );

  cout << "--- Set & Draw frame histogram " << endl;
  // Create and draw the frame
  float frameXmin=400, frameXmax=2000;
  float frameYmin=0  , frameYmax=1400;
  if( Grid=="SS_direct" ){
    frameXmin=400.; frameXmax=1800.;
    frameYmin=0.  ; frameYmax=1200.;
  }else if( Grid=="SS_onestepCC" ){
    frameXmin=200.; frameXmax=1300.;
    frameYmin=0.  ; frameYmax=900.;
  }else if( Grid=="GG_direct" ){
    frameXmin=400.; frameXmax=2400.;
    frameYmin=0.  ; frameYmax=1600.;
  }else if( Grid=="GG_onestepCC" ){
    frameXmin=600.; frameXmax=2300.;
    frameYmin=50. ; frameYmax=1500.;
  }else if( Grid=="SM_GG_N2" ){
    frameXmin=600.; frameXmax=2400.;
    frameYmin=100.; frameYmax=2500.;
  }else if( Grid=="GG_onestepCC_mlsp60" ){
    frameXmin=400.; frameXmax=2200.;
    frameYmin=0.  ; frameYmax=1.3;
  }else if( Grid=="SS_onestepCC_mlsp60" ){
    frameXmin=400.; frameXmax=2200.;
    frameYmin=0.  ; frameYmax=1.3;
  }else if( Grid=="GG_onestepN2C_mlsp60" ){
    frameXmin=800.; frameXmax=2200.;
    frameYmin=0.2 ; frameYmax=1.3;
  }else if( Grid=="SS_onestepN2C_mlsp60" ){
    frameXmin=800.; frameXmax=2100.;
    frameYmin=0.2 ; frameYmax=1.3;
  }else if( Grid=="long_lived_gl" ){
    //frameXmin=10; frameXmax=100000.;
    frameXmin=1; frameXmax=5.;
    frameYmin=600.  ; frameYmax=2800.;
  }else if( Grid.BeginsWith("SG_pMSSM_") ){
    frameXmin=800.; frameXmax=6000.;
    frameYmin=800.; frameYmax=6000.;
  }
  //TH2F *frame = new TH2F("frame", "m_#tilde{g} vs "+yParticleTex+" - ATLAS", 180, frameXmin, frameXmax, 200, frameYmin, frameYmax );
  TH2F *frame = new TH2F("frame", "", 180, frameXmin, frameXmax, 200, frameYmin, frameYmax );

  // Set common frame style
  gPad->SetTopMargin( 0.07  );
  gPad->SetBottomMargin( 0.120  );
  gPad->SetRightMargin( 0.12 );
  gPad->SetLeftMargin( 0.1  );

  //Palette
  TString STmass="";
  if       (Grid.BeginsWith("SS_")){
    STmass = "m_{#tilde{q}}";
  }else if (Grid=="long_lived_gl"){
    STmass = "lifetime";
  }else{
    STmass = "m_{#tilde{g}}";
  }

  if(Grid=="long_lived_gl"){
    frame->SetXTitle( "lifetime [ps]" );
    frame->SetYTitle( yParticleTex+" [GeV]" );
  }else{
    frame->SetXTitle( STmass+" [GeV]" );
    frame->SetYTitle( yParticleTex+" [GeV]" );
  }
  frame->SetZTitle( "X-Section" );
  frame->GetXaxis()->SetTitleOffset(1.15);
  frame->GetYaxis()->SetTitleOffset(1.15);
  frame->GetZaxis()->SetTitleOffset(1);

  frame->GetXaxis()->SetTitleFont( 42 );
  frame->GetYaxis()->SetTitleFont( 42 );
  frame->GetXaxis()->SetLabelFont( 42 );
  frame->GetYaxis()->SetLabelFont( 42 );

  frame->GetXaxis()->SetTitleSize( 0.04 );
  frame->GetYaxis()->SetTitleSize( 0.04 );
  frame->GetXaxis()->SetLabelSize( 0.035 );
  frame->GetYaxis()->SetLabelSize( 0.035 );
  frame->GetZaxis()->SetLabelSize( 0.015 );

  frame->Draw("axis same");

  // creat legend
  cout << "--- Create legend " << endl;
  TLegend *leg;
  if     (Grid=="SM_GG_N2")             leg = new TLegend(0.38,0.69,0.63,0.9);
  else if(Grid=="SS_onestepCC")         leg = new TLegend(0.40,0.69,0.63,0.9);
  else if(Grid.EndsWith("_mlsp60"))     leg = new TLegend(0.38,0.78,0.63,0.9);
  else if(Grid=="SG_pMSSM_all")         leg = new TLegend(0.40,0.53,0.7,0.81);
  else if(Grid.BeginsWith("SG_pMSSM_")) leg = new TLegend(0.45,0.58,0.7,0.81);
  else                                  leg = new TLegend(0.45,0.67,0.7,0.9);
  leg->SetTextSize( CombinationGlob::DescriptionTextSize );
  leg->SetFillStyle(0000); 
  leg->SetTextSize( 0.03 );
  leg->SetTextFont( 42 );

  Int_t c_myYellow   = TColor::GetColor("#ffe938"); 
  Int_t c_myRed      = TColor::GetColor("#aa000");
  Int_t c_myExp      = kBlue+2; //TColor::GetColor("#28373c");
  Int_t c_myGreen    = kGreen+2; //TColor::GetColor("#28373c");
  Int_t c_myMagenta  = kMagenta+1; //TColor::GetColor("#28373c");

  // color code for colz
  const Int_t ncontours=99;
  const Int_t NRGBs=4;
  Double_t stops[NRGBs] = {0.00, 0.33, 0.66, 1.00};
  Double_t red[NRGBs]   = {1.00, 238./255., 139./255., 0.00};
  Double_t green[NRGBs] = {1.00, 201./255., 90./255., 0.00};
  Double_t blue[NRGBs]  = {1.00, 0.00, 0.00, 0.00};

  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, ncontours);
  gStyle->SetNumberContours(ncontours);

  // plot p0 map in 2D 
  //contour_exp->Draw("colz same");
  TLatex latexz=TLatex();
  latexz.SetTextSize(0.035);
  latexz.SetTextFont( 42 );
  latexz.SetTextAngle(90.);
  //contour_exp->GetZaxis()->SetRangeUser(0.0, 7);
  //latexz.DrawText(2200, 20,"CL_s with best SR");

  TGraph *LEPchar = new TGraph(6);
  TGraph *targetpoint = new TGraph(3);
  TGraph *targetpoint_white = new TGraph(3);
  /*
     if (Grid=="GG_direct"){
     LEPchar->SetPoint(0,197.5,0);
     if(lumi=="8"){
     LEPchar->SetPoint(1,995,0);
     }else{
     LEPchar->SetPoint(1,1195,0);
     }
     LEPchar->SetPoint(2,650.,630.);
     LEPchar->SetPoint(3,630.,630.);
     LEPchar->SetPoint(4,197.5,197.5);
     LEPchar->SetPoint(5,197.5,0.);
     targetpoint->SetPoint(1,1100,500+20);
     targetpoint->SetPoint(2,750,650+20);
     targetpoint_white->SetPoint(1,1100,500+20);
     targetpoint_white->SetPoint(2,750,650+20);
     if(lumi=="8"){
     targetpoint->SetPoint(0,1400,0+20);
     targetpoint_white->SetPoint(0,1400,0+20);
     }else{
     targetpoint->SetPoint(0,1400,0+20);
     targetpoint_white->SetPoint(0,1400,0+20);
     }
     }else if(Grid=="SS_direct"){
     LEPchar->SetPoint(0,197.5,0);
     LEPchar->SetPoint(1,810,0);
     LEPchar->SetPoint(2,430.,430.);
     LEPchar->SetPoint(3,197.5,197.5);
     LEPchar->SetPoint(4,197.5,0.);     
     if(lumi=="8"){
     targetpoint->SetPoint(0,1200,0+20);
     targetpoint_white->SetPoint(0,1200,0+20);
     }else{
     targetpoint->SetPoint(0,1000,0+20);
     targetpoint_white->SetPoint(0,1000,0+20);
     }     
     targetpoint->SetPoint(1,800,400+20);
     targetpoint->SetPoint(2,425,375+20);
     targetpoint_white->SetPoint(1,800,400+20);
     targetpoint_white->SetPoint(2,425,375+20);
     }else if (Grid=="GG_onestepCC"){
     LEPchar->SetPoint(0,197.5,0);
     LEPchar->SetPoint(1,1195,0);
     LEPchar->SetPoint(3,400.,375.);
     LEPchar->SetPoint(2,400.,400.);
     LEPchar->SetPoint(4,197.5,197.5);
     LEPchar->SetPoint(5,197.5,0.);
     targetpoint->SetPoint(0,1545,25);
     targetpoint->SetPoint(1,1265,625);
     targetpoint->SetPoint(2,825,745);
     targetpoint_white->SetPoint(0,1545,25);
     targetpoint_white->SetPoint(1,1265,625);
     targetpoint_white->SetPoint(2,825,745);
     }
     */

  cout << "--- Draw expected band " << endl;
  c->cd();
  TString tmpGrid = Grid;
  TString labelPrefix="";
  if(Grid=="SG_pMSSM_all"){
    Grid="SG_pMSSM_0";
    labelPrefix = "m(#tilde{#chi_{1}^{0}})=0 GeV ";
  }
  TGraph* grshadeExp(0);
  if( gr_contour_exp1su && gr_contour_exp1sd ){
    grshadeExp =  (TGraph*) DrawExpectedBand( gr_contour_exp1su, gr_contour_exp1sd,   c_myYellow , 1001   , (Grid=="SG_pMSSM_695")?-9999:0 , Grid, c, removeLowMLine )->Clone();
  }
  if(!useTGraph){

    c->cd();
    DrawContourLine95( leg, contour_exp1su,         "", c_myYellow, 1, 2 , Grid, c, removeLowMLine );
    c->cd();
    DrawContourLine95( leg, contour_exp1sd,         "", c_myYellow, 1, 2 , Grid, c, removeLowMLine );

    /* duplicate observed +/- 1 sigma lines ?
       c->cd();
       if (contour_esigxsp1s!=0) DrawContourLine95( leg, contour_esigxsp1s, "", CombinationGlob::c_DarkRed, 3, 2 , Grid, c, removeLowMLine );
       c->cd();
       if (contour_esigxsm1s!=0) DrawContourLine95( leg, contour_esigxsm1s, "", CombinationGlob::c_DarkRed, 3, 2 , Grid, c, removeLowMLine );
       */

    if(unblind){
      c->cd();
      if (draw_obscls_error) DrawContourLine95( leg, contour_obscls_1, ""              , CombinationGlob::c_DarkRed, 3, 2 , Grid, c, removeLowMLine ); 
      c->cd();
      if (contour_obscls_0!=0&&draw_obscls_error) DrawContourLine95( leg, contour_obscls_0, labelPrefix+"Obs. limit (#pm1 #sigma_{theory}^{SUSY})", CombinationGlob::c_DarkRed, 1, 4 , Grid, c, removeLowMLine );
      if (contour_obscls_0!=0&&(!draw_obscls_error)) DrawContourLine95( leg, contour_obscls_0, labelPrefix+"Obs. limit", CombinationGlob::c_DarkRed, 1, 4 , Grid, c, removeLowMLine );
      c->cd();
      if (draw_obscls_error) DrawContourLine95( leg, contour_obscls_2, ""              , CombinationGlob::c_DarkRed, 3, 2 , Grid, c, removeLowMLine ); 
    }


    c->cd();
    DummyLegendExpected(leg, labelPrefix+"Exp. limits (#pm1 #sigma_{exp})",  c_myYellow, 1001, c_myExp, 6, 2);

    c->cd();
    DrawContourLine95( leg, contour_exp, "", c_myExp, 6, 2 , Grid, c, removeLowMLine );
  }else{
    c->cd();
    if(gr_contour_exp1su) DrawContourLine95( leg, gr_contour_exp1su,         "", c_myYellow, 1, 2 , Grid, c, removeLowMLine , false);
    c->cd();
    if(gr_contour_exp1sd) DrawContourLine95( leg, gr_contour_exp1sd,         "", c_myYellow, 1, 2 , Grid, c, removeLowMLine , false);


    if(unblind){
      c->cd();
      if (gobscls1!=0) DrawContourLine95( leg, gobscls1, ""              , CombinationGlob::c_DarkRed, 3, 2 , Grid, c, removeLowMLine , false); 
      c->cd();
      if (gobscls0!=0) DrawContourLine95( leg, gobscls0, labelPrefix+"Obs. limit (#pm1 #sigma_{theory}^{SUSY})", CombinationGlob::c_DarkRed, 1, 4 , Grid, c, removeLowMLine , false);
      c->cd();
      if (gobscls2!=0) DrawContourLine95( leg, gobscls2, ""              , CombinationGlob::c_DarkRed, 3, 2 , Grid, c, removeLowMLine , false); 
    }


    c->cd();
    DummyLegendExpected(leg, labelPrefix+"Exp. limits (#pm1 #sigma_{exp})",  c_myYellow, 1001, c_myExp, 6, 2);


    TFile* fMEff = TFile::Open( "output_MEff.root", "READ" );
    TFile* fRJR = TFile::Open( "output_RJR.root", "READ" );

    if(fMEff && fRJR){
      DrawContourLine95( leg, (TGraph*)fMEff->Get("CLs_Contour_0"),    "Obs. limits MEff", c_myMagenta,   1, 1 , Grid, c, removeLowMLine , false);
      DrawContourLine95( leg, (TGraph*)fMEff->Get("CLsexp_Contour_0"), "Exp. limits MEff", c_myMagenta,   7, 1 , Grid, c, removeLowMLine , false);
      DrawContourLine95( leg, (TGraph*)fRJR->Get("CLs_Contour_0"),     "Obs. limits RJR",  c_myGreen, 1, 1 , Grid, c, removeLowMLine , false);
      DrawContourLine95( leg, (TGraph*)fRJR->Get("CLsexp_Contour_0"),  "Exp. limits RJR",  c_myGreen, 7, 1 , Grid, c, removeLowMLine , false);
    }

    c->cd();
    DrawContourLine95( leg, g0_expcls, "", c_myExp, 6, 2 , Grid, c, removeLowMLine , false);


    if(unblind){
      c->cd();
      if (gobscls1!=0) DrawContourLine95( leg, gobscls1, ""              , CombinationGlob::c_DarkRed, 3, 2 , Grid, c, removeLowMLine , false); 
      c->cd();
      if (gobscls0!=0) DrawContourLine95( leg, gobscls0, "", CombinationGlob::c_DarkRed, 1, 4 , Grid, c , false);
      c->cd();
      if (gobscls2!=0) DrawContourLine95( leg, gobscls2, ""              , CombinationGlob::c_DarkRed, 3, 2 , Grid, c, removeLowMLine , false); 
    }


  }
  if(tmpGrid=="SG_pMSSM_all") Grid="SG_pMSSM_all";


  // draw other mlsp mass grids
  if( Grid=="SG_pMSSM_all" ){
    int colors[3]={CombinationGlob::c_DarkRed, CombinationGlob::c_LightBlue, CombinationGlob::c_DarkGreen};
    for( int i=1 ; i<gobscls0all.size() ; i++ ){
      c->cd();
      if(unblind){
        TString labelObs="m({#tilde{#chi_{1}^{0}}})="+std::to_string(mlsps[i])+" GeV, Obs. limits";
        if (gobscls0all[i]!=0) DrawContourLine95( leg, gobscls0all[i], labelObs, colors[i], 1, 4, Grid, c , false);
      }
      TString labelExp="m({#tilde{#chi_{1}^{0}}})="+std::to_string(mlsps[i])+" GeV, Exp. limits";
      if(g0_expclsall[i]!=0) DrawContourLine95( leg, g0_expclsall[i], labelExp, colors[i], 6, 2, Grid, c , false);
    }
  }

  if( drawOldLimits ){
    // plot 8 TeV exclusion:
    cout << "--- Get 8 TeV exclusionend " << endl;
    TGraphAsymmErrors *g(0);
    if(Grid=="GG_direct"){
      //g=GG_direct_8TeVobs();
      ;
    }else if (Grid=="SS_direct"){
      //g=SS_direct_8TeVobs();
      ;
    }else if (Grid=="SS_onestepCC"){
      // AW uncommented 
      //g=SS_onestepCC_8TeVobs();
      //TFile file("contourmacros/SS_onestep_x05_mergedHists_8TeVobs.root");
      //g=(TGraphAsymmErrors*)file.Get("combined_observed_0");
      g=GetLimitGraph("Run1result/SS_onestepCC_Run1.root","observed");
      //TFile file("Run1result/SS_onestepCC_Run1.root");
      //g=(TGraphAsymmErrors*)file.Get("observed");
      ;
    }else if (Grid=="GG_onestepCC"){
      //      g=GG_onestepCCx12_8TeVobs(); 
      ;
    }else if (Grid=="GG_onestepCC_mlsp60"){
      //g=GG_onestep_mlsp60_8TeVobs();
      //g=GetLimitGraph("Run1result/GG_onestep_mlsp60_Run1.root","combined_observed_0");
      ;
    }else if (Grid=="SG_pMSSM_0"||Grid=="SG_pMSSM_all"){
      g=SG_pMSSM_0_Run1();
    }else if (Grid=="SG_pMSSM_395"){
      g=SG_pMSSM_395_Run1();
    }else if (Grid=="SG_pMSSM_695"){
      g=SG_pMSSM_695_Run1();
    }
    c->cd();
    if(g){
      cout << "--- Draw graph : " << g->GetName()<< endl;
      g->SetName("gRun1");

      // Modify line
      double x=0.,y=0.;
      for( int n=0 ; n<g->GetN() ; n++ ){
        x=0., y=0.;
        g->GetPoint(n,x,y);
        // Modify points forming low mass line
        if(removeLowMLine) modifyPointToRemoveLowMLine(x,y,frameXmin,frameXmax,frameYmin,Grid);
        // Modify points above forbidden line
        if( Grid.EndsWith("_mlsp60") ){
          if(y>1.)  y=1.;
        }else if(Grid!="long_lived_gl"&&(!Grid.BeginsWith("SG_pMSSM_"))){
          if(y>x)   y=x;
        }
        g->SetPoint(n,x,y);
        cout<<"  g : -->> set point (x,y)=("<<x<<","<<y<<")"<<endl;
      }

      g->SetLineStyle(1);

      //g->SetFillColor(kBlue-10);
      g->SetFillColorAlpha(kGray,0.5);
      g->SetLineWidth(1);
      g->SetFillStyle(1001);
      //g->SetFillColor(kBlue-7);
      //g->SetLineWidth(2);
      //g->SetFillStyle(3490);

      g->SetLineColor(kGray+1);
      //g->SetLineColor(17);
      g->Draw("F same");
      g->Draw("L same");
      if (Grid=="SS_onestepCC")
	leg->AddEntry(g,"0L+1L combination obs. limit (20.3 fb^{-1}, 8 TeV "+labelPrefix+")","fl");
      else
	leg->AddEntry(g,"Observed limit (20.3 fb^{-1}, 8 TeV "+labelPrefix+")","fl");
      TString outfileName = Grid+"_contours.root" ;
      if(!outfile) outfile = new TFile(outfileName, "recreate");
      outfile->cd();
      g->Write();
    }


    //begin AW uncomment
    // plot 2015 exclusion:
    cout << "--- Get 2015 exclusionend " << endl;
    TGraphAsymmErrors *g2015(0);
    if(Grid=="GG_direct"){
      //g2015=GetLimitGraph("2015result/GG_direct_2015.root","observed");
    }else if (Grid=="SS_direct"){
      //      g2015=GetLimitGraph("2015result/SS_direct_2015.root","observed");
    }else if (Grid=="SS_onestepCC"){
    //g2015=GetLimitGraph("2015result/SS_onestepCC_2015.root","observed");
    }else if (Grid=="GG_onestepCC"){
    g2015=GetLimitGraph("2015result/GG_onestepCC_2015.root","observed");
    }else if (Grid=="SM_GG_N2"){
    //g2015=GetLimitGraph("2015result/SM_GG_N2_2015.root","observed");
    TFile* f2015  = new TFile("2015result/SM_GG_N2_2015_Z+MET.root");
    if(f2015){
    TH2D* hist2015 = (TH2D*)f2015->Get("sigp1clsf;2");
    if(doSmooth)  hist2015->Smooth();
    TH2F* contour_2015(0);
    if(hist2015) { contour_2015 = FixAndSetBorders( *hist2015, "", "", 0 ); }
    if(contour_2015) g2015 = GetContourLine95( contour_2015, Grid );
    }
    }
    c->cd();
    if(g2015){
    cout << "--- Draw graph : " << g2015->GetName()<< endl;
    g2015->SetName("g2015");
    // Modify line
    double x=0.,y=0.;
    for( int n=0 ; n<g2015->GetN() ; n++ ){
    x=0., y=0.;
    g2015->GetPoint(n,x,y);
    // Modify points forming low mass line
    if(removeLowMLine) modifyPointToRemoveLowMLine(x,y,frameXmin,frameXmax,frameYmin,Grid);
    // Modify points above forbidden line
    if( Grid.EndsWith("_mlsp60") ){
    if(y>1.)  y=1.;
    }else if(Grid!="long_lived_gl"){
    if(y>x)   y=x;
    }
    g2015->SetPoint(n,x,y);
    cout<<"  g2015 : -->> set point (x,y)=("<<x<<","<<y<<")"<<endl;
    }

    g2015->SetLineStyle(1);

    //g2015->SetFillColor(kBlue-10);
    g2015->SetFillColorAlpha(kAzure-4,0.5);
    g2015->SetLineColor(kAzure-5);
    g2015->SetLineWidth(1);
    g2015->SetFillStyle(1001);
    //g2015->SetFillColor(kBlue-7);
    //g2015->SetLineWidth(2);
    //g2015->SetFillStyle(3490);

    //g2015->SetLineColor(17);
    g2015->Draw("F same");
    g2015->Draw("L same");
    if(Grid=="SM_GG_N2"){
    leg->AddEntry(g2015,      "2L+MET analysis obs. limit (3.2 fb^{-1}, 2015)","fl");
    }else{
    leg->AddEntry(g2015,"Obs. limit (3.2 fb^{-1}, 2015)","fl");
    }
    TString outfileName = Grid+"_contours.root" ;
    if(!outfile) outfile = new TFile(outfileName, "recreate");
    outfile->cd();
    g2015->Write();
    }
    // end AW uncomment


    // plot 2016 exclusion:
    cout << "--- Get 2016 exclusionend " << endl;
    TGraphAsymmErrors *g2016(0);
    TGraphAsymmErrors *g2016ZMET(0);
    if(Grid=="GG_direct"){
      g2016=GetLimitGraph("2016result/GG_direct_contours.root","gobserved");
    }else if (Grid=="SS_direct"){
      g2016=GetLimitGraph("2016result/SS_direct_contours.root","gobserved");
    }else if (Grid=="SS_onestepCC"){
      //g2016=GetLimitGraph("2016result/SS_onestepCC_contours.root","gobserved");
    }else if (Grid=="GG_onestepCC"){
      //g2016=GetLimitGraph("2016result/GG_onestepCC_contours.root","gobserved");
    }else if (Grid=="SM_GG_N2"){
      g2016=GetLimitGraph("2016result/SM_GG_N2_contours.root","gobserved");
      g2016ZMET=GetLimitGraph("2016result/SM_GG_N2_contours_ZMET.root","contour_obscls");
    }
    c->cd();
    if(g2016){
      cout << "--- Draw graph : " << g2016->GetName()<< endl;
      g2016->SetName("g2016");
      // Modify line
      double x=0.,y=0.;
      for( int n=0 ; n<g2016->GetN() ; n++ ){
        x=0., y=0.;
        g2016->GetPoint(n,x,y);
        // Modify points forming low mass line
        if(removeLowMLine) modifyPointToRemoveLowMLine(x,y,frameXmin,frameXmax,frameYmin,Grid);
        // Modify points above forbidden line
        if( Grid.EndsWith("_mlsp60") ){
          if(y>1.)  y=1.;
        }else if(Grid!="long_lived_gl"){
          if(y>x)   y=x;
        }
        g2016->SetPoint(n,x,y);
        cout<<"  g2016 : -->> set point (x,y)=("<<x<<","<<y<<")"<<endl;
      }

      g2016->SetLineStyle(1);

      //g2016->SetFillColor(kBlue-10);
      g2016->SetFillColorAlpha(kAzure-4,0.5);
      g2016->SetLineColor(kAzure-5);
      g2016->SetLineWidth(1);
      g2016->SetFillStyle(1001);
      //g2016->SetFillColor(kBlue-7);
      //g2016->SetLineWidth(2);
      //g2016->SetFillStyle(3490);

      //g2016->SetLineColor(17);
      g2016->Draw("F same");
      g2016->Draw("L same");
      leg->AddEntry(g2016,"Obs. limit (13.3 fb^{-1}, 2016)","fl");
      TString outfileName = Grid+"_contours.root" ;
      if(!outfile) outfile = new TFile(outfileName, "recreate");
      outfile->cd();
      g2016->Write();
    }
    if(g2016ZMET){
      cout << "--- Draw graph : " << g2016ZMET->GetName()<< endl;
      g2016ZMET->SetName("g2016ZMET");
      // Modify line
      double x=0.,y=0.;
      for( int n=0 ; n<g2016ZMET->GetN() ; n++ ){
        x=0., y=0.;
        g2016ZMET->GetPoint(n,x,y);
        // Modify points forming low mass line
        if(removeLowMLine) modifyPointToRemoveLowMLine(x,y,frameXmin,frameXmax,frameYmin,Grid);
        if(y>x) y=x;
        g2016ZMET->SetPoint(n,x,y);
        cout<<"  g2016ZMET : -->> set point (x,y)=("<<x<<","<<y<<")"<<endl;
      }

      g2016ZMET->SetLineStyle(1);

      g2016ZMET->SetFillColorAlpha(kGray,0.5);
      g2016ZMET->SetLineColor(kGray+1);
      g2016ZMET->SetLineWidth(1);
      g2016ZMET->SetFillStyle(1001);

      g2016ZMET->Draw("F same");
      g2016ZMET->Draw("L same");
      leg->AddEntry(g2016ZMET,"Obs. limit (14.7 fb^{-1}, Z+MET analysis 2016)","fl");
      TString outfileName = Grid+"_contours.root" ;
      if(!outfile) outfile = new TFile(outfileName, "recreate");
      outfile->cd();
      g2016ZMET->Write();
    }



  } // end of drawOldLimits

  // Draw each SR's expected limits
  if( showAllSRExp ){
    int colors[]={kRed,kBlue,kGreen,kMagenta,kCyan,kPink-4,kAzure-4,kSpring-4,kRed+3,kBlue+3,kGreen+3,kOrange+7,kViolet+7,kYellow+2,kGray+2};
    //TLegend * legSREach=new TLegend(0.9,0.1,1.0,0.9);
    TLegend * legSREach=new TLegend(0.73,0.5,0.93,0.9);
    for( size_t i=0 ; i<SRNames.size() ; i++ ){
      TFile* file(0);
      TString fname0Each = fname0;
      TString SRName=SRNames.at(i);
      cout<<"SR : "<<SRName<<endl;
      fname0Each.ReplaceAll("combined",SRName);
      file = new TFile( fname0Each, "READ" );
      if( !file || !file->IsOpen() ){
        cout<<"Warning! Failed to open file : "<<fname0Each<<endl;
        continue;
      }
      TH2F* hist0_expclsEach = (TH2F*)file->Get(excl_p1exp);
      if( doSmooth ) hist0_expclsEach->Smooth();
      TH2F* contour_expEach(0);
      if (hist0_expclsEach!=0) { contour_expEach = FixAndSetBorders( *hist0_expclsEach, "expected", "contour_exp"+SRName, 0 ); }
      c->cd();
      DrawContourLine95( legSREach, contour_expEach, SRName, colors[i], 3, 2, Grid, c, false, false );
    }
    c->cd();
    legSREach->Draw();
    cout<<"finish drawing all of SR's expected limits"<<endl;
  }


  c->cd();
  // uncertainty lines in observed limit legend
  TLine* lineObsLimitUnc = new TLine();
  int nLegEntry = leg->GetNRows();
  float legX1=leg->GetX1();
  float legX2=leg->GetX2();
  float legY1=leg->GetY1();
  float legY2=leg->GetY2();
  float lineY =legY2-(legY2-legY1)*0.5/(float)nLegEntry;
  float deltaY=(legY2-legY1)*0.2/(float)nLegEntry;
  float marginX=(legX2-legX1) * leg->GetMargin();
  float Xspace =marginX*0.15;
  cout<<"nLegEntry="<<nLegEntry<<endl;
  cout<<"leg X="<<legX1<<"--"<<legX2<<endl;
  cout<<"leg Y="<<legY1<<"--"<<legY2<<endl;
  cout<<"marginX="<<marginX<<endl;

  if(draw_obscls_error){
    cout<<"line Y ="<<lineY<<endl;
    cout<<"deleta leg Y="<<deltaY<<endl;
    lineObsLimitUnc->SetLineStyle(3);
    lineObsLimitUnc->SetLineWidth(1);
    lineObsLimitUnc->SetLineColor(CombinationGlob::c_DarkRed);
    lineObsLimitUnc->DrawLineNDC(legX1+Xspace, lineY+deltaY, legX1+marginX-Xspace, lineY+deltaY ) ;                
    lineObsLimitUnc->DrawLineNDC(legX1+Xspace, lineY-deltaY, legX1+marginX-Xspace, lineY-deltaY ) ;                
    lineObsLimitUnc->DrawLine(0,0.5,1,0.5) ;                
  }

  LEPchar->SetFillColor(17);
  //LEPchar->Draw("FL same");

  targetpoint->SetMarkerSize(1.5);
  targetpoint->SetMarkerStyle(29);
  //targetpoint->Draw("P same");

  targetpoint_white->SetMarkerSize(1.5);
  targetpoint_white->SetMarkerStyle(30);
  targetpoint_white->SetMarkerColor(0);
  //targetpoint_white->Draw("P same");

  TObjArray* arr = fname0.Tokenize("/");
  TObjString* objstring = (TObjString*)arr->At( arr->GetEntries()-1 );
  TString outfile = objstring->GetString().ReplaceAll(".root","").ReplaceAll("_1_harvest_list","");

  // mass forbidden line
  TLine lineExcl;
  int maxframeMin=TMath::Max(frameXmin,frameYmin);
  int minframeMax=TMath::Min(frameXmax,frameYmax);
  if( Grid=="SM_GG_N2" ){
    lineExcl = ( TLine(maxframeMin,maxframeMin,minframeMax-600,minframeMax-600));
  }else if( Grid.EndsWith("_mlsp60") ){
    lineExcl = ( TLine(frameXmin,1.,frameXmax,1.));
  }else if( Grid=="GG_direct" ){
    lineExcl = ( TLine(maxframeMin,maxframeMin,1000,1000));
  }else{
    lineExcl = ( TLine(maxframeMin,maxframeMin,800,800));
  }
  lineExcl.SetLineStyle(9);
  lineExcl.SetLineWidth(1);
  lineExcl.SetLineColor(14);

  /*
  // cover above mass forbidden line
  TGraph* gCover = new TGraph();
  gCover->SetLineColor(0);
  gCover->SetMarkerColor(1);
  gCover->SetMarkerSize(1);
  gCover->SetMarkerStyle(20);
  gCover->SetFillColor(10);
  gCover->SetFillStyle(1001);
  int nPoint=0;
  gCover->SetPoint(nPoint++,frameXmin,frameYmax);
  gCover->SetPoint(nPoint++,maxframeMin,maxframeMin);
  gCover->SetPoint(nPoint++,minframeMax,minframeMax);
  gCover->SetPoint(nPoint++,frameXmin+(minframeMax-frameXmin)/2.,frameYmax);
  gCover->SetPoint(nPoint++,frameXmin,frameYmax);
  */

  // mass forbidden label
  float massforbiddenX=500;
  float massforbiddenY=580;
  float angle=45.;
  TString massforbiddenLabel = STmass+" < "+yParticleTex ; 
  if( Grid=="SM_GG_N2" ){
    massforbiddenX=700;
    massforbiddenY=750;
    angle=33;
  }else if(Grid=="SS_direct"){
    angle=40;
  }else if( Grid.EndsWith("_mlsp60") ){
    massforbiddenX=1800;
    massforbiddenY=1.05;
    massforbiddenLabel = "x=1";
    angle=0;
  }
  TLatex massforbidden = TLatex(massforbiddenX,massforbiddenY,massforbiddenLabel);
  massforbidden.SetTextSize(0.032);
  massforbidden.SetTextColor(14);
  massforbidden.SetTextAngle(angle);
  massforbidden.SetTextFont(42);

  if( Grid!="long_lived_gl" && !Grid.BeginsWith("SG_pMSSM_") ){
    //gCover->Draw("sameFL");
    lineExcl.Draw("same") ;
    massforbidden.Draw("same"); 
  }

  // Cosmetics for plotting 
  Float_t textSizeOffset = +0.000;
  Double_t xmax = frame->GetXaxis()->GetXmax();
  Double_t xmin = frame->GetXaxis()->GetXmin();
  Double_t ymax = frame->GetYaxis()->GetXmax();
  Double_t ymin = frame->GetYaxis()->GetXmin();
  Double_t dx   = xmax - xmin;
  Double_t dy   = ymax - ymin;

  float labAx=0.13, labAy=0.87; // label of "ATLAS"
  float lab2x=0.13, lab2y=0.74; // label of prefix(0-leptons,2-6jets)
  float lab3x=0.13, lab3y=0.7;// label of "All limits at 95% CL."
  if(Grid=="long_lived_gl"){
    // add log x-axis
    TGaxis* xaxis = new TGaxis(frameXmin,frameYmin,frameXmax,frameYmin,1e+1,1e+5,510,"G");
    xaxis->SetLabelSize(frame->GetXaxis()->GetLabelSize());
    xaxis->SetLabelFont(frame->GetXaxis()->GetLabelFont());
    xaxis->SetTitle(frame->GetXaxis()->GetTitle());
    xaxis->SetTitleSize(frame->GetXaxis()->GetTitleSize());
    xaxis->SetTitleFont(frame->GetXaxis()->GetTitleFont());
    xaxis->SetTitleOffset(frame->GetXaxis()->GetTitleOffset());
    xaxis->Draw();
    TGaxis* xaxisTop = new TGaxis(frameXmin,frameYmax,frameXmax,frameYmax,1e+1,1e+5,510,"G-");
    xaxisTop->SetLabelSize(0);
    xaxisTop->SetTitle("");
    xaxisTop->Draw();
    // remove non-log x-axis
    frame->GetXaxis()->SetTitleOffset(999);
    frame->GetXaxis()->SetLabelSize(0);
    frame->GetXaxis()->SetTickSize(0);
  }else if(Grid.EndsWith("_mlsp60")){
    lab2y=labAy-0.08;
    lab3y=lab2y-0.04;
  }else if(Grid.BeginsWith("SG_pMSSM_")){
    labAx=leg->GetX1();
    labAy=leg->GetY2()+0.06;
    lab2x=leg->GetX1()+0.21;
    lab2y=labAy;
    lab3x=lab2x;
    lab3y=lab2y-0.04;
  }

  TLatex *Leg0 = new TLatex( xmin, ymax + dy*0.025, "" );
  Leg0->SetTextAlign( 11 );
  Leg0->SetTextFont( 42 );
  Leg0->SetTextSize( CombinationGlob::DescriptionTextSize);
  Leg0->SetTextColor( 1 );
  Leg0->AppendPad();

  TString STdecayLabel="";
  TString STdecayLabel1="";
  TString STdecayLabel2="";
  TString mLSP=""; // only for SG_pMSSM
  if (Grid=="GG_direct"){
    STdecayLabel="#tilde{g}#tilde{g} production, #it{B}(#tilde{g} #rightarrow qq #tilde{#chi}_{1}^{0})=100%";
  }else if (Grid=="SS_direct"){
    STdecayLabel="#tilde{q}#tilde{q} production, #it{B}(#tilde{q} #rightarrow q #tilde{#chi}_{1}^{0})=100%";
  }else if (Grid=="SS_onestepCC"){
    STdecayLabel="#tilde{q}#tilde{q} production, #it{B}(#tilde{q} #rightarrow q #tilde{#chi}_{1}^{#pm} #rightarrow q W^{#pm} #tilde{#chi}_{1}^{0})=100%, m(#tilde{#chi}_{1}^{#pm})=(m(#tilde{q}) + m(#tilde{#chi}_{1}^{0}))/2";
  }else if (Grid=="SS_onestepCC_mlsp60"){
    STdecayLabel="#tilde{q}#tilde{q} production, #it{B}(#tilde{q} #rightarrow q #tilde{#chi}_{1}^{#pm} #rightarrow q W^{#pm} #tilde{#chi}_{1}^{0})=100%, m(#tilde{#chi}_{1}^{0})=60 GeV";
  }else if (Grid=="SS_onestepN2C_mlsp60"){
    STdecayLabel ="#tilde{q}#tilde{q} production, #it{B}(#tilde{q} #rightarrow q #tilde{#chi}_{1}^{#pm} #rightarrow q W^{#pm} #tilde{#chi}_{1}^{0})=50%, #it{B}(#tilde{q} #rightarrow q #tilde{#chi}_{2}^{0} #rightarrow q Z/h #tilde{#chi}_{1}^{0})=50%, m(#tilde{#chi}_{1}^{0})=60 GeV";
  }else if (Grid=="GG_onestepCC"){
    STdecayLabel="#tilde{g}#tilde{g} production, #it{B}(#tilde{g} #rightarrow qq #tilde{#chi}_{1}^{#pm} #rightarrow qq W^{#pm} #tilde{#chi}_{1}^{0})=100%, m(#tilde{#chi}_{1}^{#pm})=(m(#tilde{g}) + m(#tilde{#chi}_{1}^{0}))/2";
  }else if (Grid=="GG_onestepCC_mlsp60"){
    STdecayLabel="#tilde{g}#tilde{g} production, #it{B}(#tilde{g} #rightarrow qq #tilde{#chi}_{1}^{#pm} #rightarrow qq W^{#pm} #tilde{#chi}_{1}^{0})=100%, m(#tilde{#chi}_{1}^{0})=60 GeV";
  }else if (Grid=="GG_onestepN2C_mlsp60"){
    STdecayLabel ="#tilde{g}#tilde{g} production, #it{B}(#tilde{g} #rightarrow qq #tilde{#chi}_{1}^{#pm} #rightarrow qq W^{#pm} #tilde{#chi}_{1}^{0})=50%, #it{B}(#tilde{g} #rightarrow qq #tilde{#chi}_{2}^{0} #rightarrow qq Z/h #tilde{#chi}_{1}^{0})=50%, m(#tilde{#chi}_{1}^{0})=60 GeV";
  }else if (Grid=="SM_GG_N2"){
    STdecayLabel="#tilde{g}#tilde{g} production, #it{B}(#tilde{g} #rightarrow qq #tilde{#chi}_{2}^{0} #rightarrow qq Z #tilde{#chi}_{1}^{0})=100%, m(#tilde{#chi}_{1}^{0})=1 GeV";
  }else if (Grid=="long_lived_gl"){
    STdecayLabel="#tilde{g}#tilde{g} production, #it{B}(#tilde{g} #rightarrow qq #tilde{#chi}_{1}^{0})=100%, m(#tilde{#chi}_{1}^{0})=100 GeV";
  }else if (Grid.BeginsWith("SG_pMSSM_")){
    std::cout<<" string size = "<<sizeof("SG_pMSSM_")/sizeof(char)<<std::endl;
    mLSP = Grid; mLSP.Remove(0,sizeof("SG_pMSSM_")/sizeof(char)-1);
    if(mLSP=="all") STdecayLabel="Squark-gluino-neutralino model";
    else            STdecayLabel="Squark-gluino-neutralino model, m(#tilde{#chi}_{1}^{0})="+mLSP+" GeV";
  }
  TLatex* decayLabel = new TLatex();
  decayLabel->SetNDC();
  decayLabel->SetTextFont(42);
  decayLabel->SetTextColor(EColor::kBlack);
  decayLabel->SetTextSize( 0.0335 );
  // small font
  if(STdecayLabel1!=""&&STdecayLabel2!="") {
    decayLabel->DrawLatex(0.14, 0.962,STdecayLabel);
    decayLabel->SetTextSize( 0.024 );
    decayLabel->DrawLatex(0.30, 0.975,STdecayLabel1);
    decayLabel->DrawLatex(0.30, 0.945,STdecayLabel2);
  }else if(Grid.Contains("N2C")){
    decayLabel->DrawLatex(0.01, 0.962,STdecayLabel);
  }else{
    decayLabel->DrawLatex(0.14, 0.962,STdecayLabel);
  }

  TLatex* atlasLabel = new TLatex();
  atlasLabel->SetNDC();
  atlasLabel->SetTextFont(42);
  atlasLabel->SetTextColor(EColor::kBlack);
  atlasLabel->SetTextSize( 0.03 );
  atlasLabel->DrawLatex(labAx, labAy,"#bf{#it{ATLAS}}");
  atlasLabel->AppendPad();

  TLatex* progressLabel = new TLatex();
  progressLabel->SetNDC();
  progressLabel->SetTextFont(42);
  progressLabel->SetTextColor(EColor::kBlack);
  progressLabel->SetTextSize( 0.03 );
  progressLabel->DrawLatex(labAx+0.07, labAy,"Internal");
  progressLabel->AppendPad();

  TLatex *Leg1 = new TLatex();
  Leg1->SetNDC();
  Leg1->SetTextFont( 42 );
  Leg1->SetTextSize( 0.03 );
  Leg1->SetTextColor( EColor::kBlack );
  Leg1->DrawLatex(labAx, labAy-0.04, "#sqrt{s} = 13 TeV, "+lumi+" fb^{-1}");
  Leg1->AppendPad();

  TLatex *Leg2 = new TLatex();
  Leg2->SetNDC();
  Leg2->SetTextSize(0.03);
  Leg2->SetTextColor( 1 );
  Leg2->SetTextFont( 42 );
  Leg2->DrawLatex(lab2x, lab2y, prefix);
  if (prefix!=0) { Leg2->AppendPad(); }

  TLatex *Leg3 = new TLatex();
  Leg3->SetNDC();
  Leg3->SetTextFont( 42 );
  Leg3->SetTextSize( 0.03 );
  Leg3->SetTextColor( EColor::kBlack );
  Leg3->DrawLatex(lab3x, lab3y, "All limits at 95% CL.");
  Leg3->AppendPad();

  TLatex *Leg4 = new TLatex();
  Leg4->SetNDC();
  Leg4->SetTextFont( 42 );
  Leg4->SetTextSize( 0.03 );
  Leg4->SetTextColor( EColor::kBlack );
  Leg4->DrawLatex(0.13, 0.65, SRLabel);
  Leg4->AppendPad();

  TLatex *Leg5 = new TLatex();
  Leg5->SetNDC();
  Leg5->SetTextFont( 42 );
  Leg5->SetTextSize( 0.03 );
  Leg5->SetTextColor( EColor::kBlack );
  if( Grid.EndsWith("CC_mlsp60") ){
    Leg5->DrawLatex(0.13, 0.65, "m_{#tilde{#chi}_{1}^{0}} = 60 [GeV]");
    Leg5->AppendPad();
  }


  if (prefix!=0) { Leg2->AppendPad(); }

  if(showSR){
    std::cout << "--- printing best SRs" << std::endl;
    //if(Grid=="GG_onestepCC") Show_SR(fname0, c, frame, xmin, xmax, ymin, ymax-400, Grid, leg);
    //else                     Show_SR(fname0, c, frame, xmin, xmax, ymin, ymax    , Grid, leg);
    Show_SR(fname0, c, frame, xmin, xmax, ymin, ymax, Grid, leg);
    TLatex lat;
    //lat.SetTextAlign( 11 );
    lat.SetTextSize( 0.0265 );
    lat.SetTextColor( 12 );
    lat.SetTextFont( 42 ); 
  }

  frame->Draw("sameaxis"); 
  leg->Draw("same");

  c->Update();
  TString suffix=(showSR) ? "_showSR" : "";
  suffix+=(showAllSRExp) ? "_withAllSRExp" : "";
  TString plotname = outdir+"/atlasCLs_"+Grid+suffix+".eps";
  c->SaveAs(plotname);
  c->SaveAs(plotname.ReplaceAll(".eps",".pdf"));
  c->SaveAs(plotname.ReplaceAll(".pdf",".root"));


}


