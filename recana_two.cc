#include <iostream>
#include <string>
#include <TApplication.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TText.h>
#include <TLatex.h>

using namespace std;


double Peak(TH1D *h);

int main(int argc, char* argv[]){
  if(argc!=2){
    cout << "usage : ./recana_two (input filename)" << endl;
    return 1;
  }
  string filename = string(argv[1]);

  TApplication app("app", &argc, argv);

  TFile *ifile = new TFile(filename.c_str());
  TTree *photonsum_two_tree = (TTree*)ifile -> Get("photonsum_two_tree");
  int entry_num = photonsum_two_tree -> GetEntries();
  int photon_sum = 0;
  photonsum_two_tree -> SetBranchAddress("phsum",&photon_sum);
  
  TTree *photonsum_two_header = (TTree*)ifile -> Get("photonsum_two_header");
  double tau1 = 0., tau2 = 0., tau = 0.;
  double alpha = 0., beta = 0.;
  photonsum_two_header -> SetBranchAddress("tau1", &tau1);
  photonsum_two_header -> SetBranchAddress("tau2", &tau2);
  photonsum_two_header -> SetBranchAddress("alpha", &alpha);
  photonsum_two_header -> SetBranchAddress("beta", &beta);
  photonsum_two_header -> SetBranchAddress("tau", &tau);
  photonsum_two_header -> GetEntry(0);

  TH1D *h_phsum;
  h_phsum = new TH1D("phsum", "phsum", 8000,0,800e3);    

  for(int ev=0; ev<entry_num; ev++){
    photonsum_two_tree -> GetEntry(ev);
    h_phsum -> Fill(photon_sum);
  }

  TF1 *f_phsum;
  f_phsum = new TF1("phsum","gaus",0,700e3);
 
  double fit_mean = 0;
  double fit_sigma = 0;
  double fit_fwhm = 0;
  double fit_peak = 0;
  fit_peak = Peak(h_phsum);
  
  TCanvas *c_phsum = new TCanvas("phsum","phsum",100,100,500,250);
  c_phsum -> Divide(2,1,0.01,0.01);
  gStyle -> SetStatX(0.2);

  c_phsum -> cd(1) -> SetLogy(1);
  h_phsum -> Draw();
  h_phsum -> Fit("phsum","Q","",fit_peak-50e3,fit_peak+50e3);
  fit_mean = f_phsum -> GetParameter(1);
  fit_sigma = f_phsum -> GetParameter(2);
  fit_fwhm = fit_sigma/fit_mean*2.35*100;
  c_phsum -> cd(2);
  c_phsum -> DrawFrame(630e3, 0, 680e3, 300);
  h_phsum -> Draw("same");
  TText *tx_rectime_two = new TText(0.25,0.85,Form("tau1 = %.1f ns, tau2 = %.0f", tau1, tau2));
  TText *tx_albe = new TText(0.25,0.8,Form("(alpha = %.2f, beta = %.2f)", alpha, beta));
  TText *tx_rectime = new TText(0.25,0.75,Form("tau = %.1f ns", tau));
  TText *tx_mean = new TText(0.25,0.7,Form("mean = %.0f",fit_mean));
  TText *tx_sigma = new TText(0.25,0.65,Form("sigma = %.0f",fit_sigma));
  TText *tx_fwhm = new TText(0.25,0.6,Form("FWHM %.3f %%",fit_fwhm));
  tx_rectime_two -> SetNDC(1);
  tx_albe -> SetNDC(1);
  tx_rectime -> SetNDC(1);
  tx_mean -> SetNDC(1);
  tx_sigma -> SetNDC(1);
  tx_fwhm -> SetNDC(1);
  tx_rectime_two -> Draw("same");
  tx_albe -> Draw("same");
  tx_rectime -> Draw("same");
  tx_mean -> Draw("same");
  tx_sigma -> Draw("same");
  tx_fwhm -> Draw("same");
  
  app.Run();
  
  //photonsum_two_tree -> Draw("phsum_0ns");

  ifile -> Close();



  return 0;
}


double Peak(TH1D *h){
  int all_bin_num = h -> GetXaxis() -> GetNbins();
  int bin_content = 0;
  int max_bin_content = 0;
  int max_bin_num = 0;
  for(int i=0; i<all_bin_num; i++){
    bin_content = h->GetBinContent(i);
    if(bin_content > max_bin_content) {
      bin_content = max_bin_content;
      max_bin_num = i;
    }
  }
  int xmin = h -> GetXaxis() -> GetXmin();
  int xmax = h -> GetXaxis() -> GetXmax();
  return (double)xmin + (double)max_bin_num*(xmax - xmin)/(double)all_bin_num;
}
