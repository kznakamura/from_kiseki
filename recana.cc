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

const double rec_array[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 20.0, 40.0, 60.0};
const size_t rec_array_size = end(rec_array) - begin(rec_array);

double Peak(TH1D *h);

int main(int argc, char* argv[]){
  if(argc!=2){
    cout << "usage : ./recana (input filename)" << endl;
    return 1;
  }
  string filename = string(argv[1]);

  TApplication app("app", &argc, argv);

  TFile *ifile = new TFile(filename.c_str());
  TTree *photonsum_tree = (TTree*)ifile -> Get("photonsum_tree");

  int entry_num = photonsum_tree -> GetEntries();

  int photon_sum[rec_array_size] = {};
  for(size_t elem=0; elem<rec_array_size; elem++){
    photonsum_tree -> SetBranchAddress(Form("phsum_%dns",(int)rec_array[elem]),&photon_sum[elem]);
  }

  TH1D *h_phsum[rec_array_size];
  for(size_t elem=0; elem<rec_array_size; elem++){
    h_phsum[elem] = new TH1D(Form("phsum_%dns",(int)rec_array[elem]), Form("phsum_%dns",(int)rec_array[elem]), 8000,0,800e3);
    
  }

  for(int ev=0; ev<entry_num; ev++){
    photonsum_tree -> GetEntry(ev);
    for(size_t elem=0; elem<rec_array_size; elem++){
      h_phsum[elem] -> Fill(photon_sum[elem]);
      //   cout << photon_sum[elem] << endl;
    }
  }

  TF1 *f_phsum[rec_array_size];
  for(size_t elem=0; elem<rec_array_size; elem++){
    f_phsum[elem] = new TF1(Form("phsum_%dns",(int)rec_array[elem]),"gaus",0,700e3);
  }
  double fit_mean[rec_array_size] = {};
  double fit_sigma[rec_array_size] = {};
  double fit_fwhm[rec_array_size] = {};
  double fit_peak[rec_array_size] = {};
  for(size_t elem=0; elem<rec_array_size; elem++){
    fit_peak[elem] = Peak(h_phsum[elem]);
  }
  
  TCanvas *c_phsum = new TCanvas("phsum","phsum",100,100,1250,500);
  c_phsum -> Divide(5,2,0.01,0.01);
  gStyle -> SetStatX(0.2);
  TCanvas *c_phsum_q = new TCanvas("phsum_q","phsum_q",100,400,1250,500);
  c_phsum_q -> Divide(5,2,0.01,0.01);

  TText *tx_rectime[rec_array_size];
  TText *tx_mean[rec_array_size];
  TText *tx_sigma[rec_array_size];
  TText *tx_fwhm[rec_array_size];

  for(size_t elem=0; elem<rec_array_size; elem++){
    c_phsum -> cd(elem+1) -> SetLogy(1);
    h_phsum[elem] -> Draw();
    h_phsum[elem] -> Fit(Form("phsum_%dns",(int)rec_array[elem]),"Q","",fit_peak[elem]-50e3,fit_peak[elem]+50e3);
    fit_mean[elem] = f_phsum[elem] -> GetParameter(1);
    fit_sigma[elem] = f_phsum[elem] -> GetParameter(2);
    fit_fwhm[elem] = fit_sigma[elem]/fit_mean[elem]*2.35*100;
    c_phsum_q -> cd(elem+1) -> DrawFrame(630e3, 0, 680e3, 300);
    h_phsum[elem] -> Draw("same");
    tx_rectime[elem] = new TText(0.25,0.8,Form("rectime = %d ns",(int)rec_array[elem]));
    tx_mean[elem] = new TText(0.25,0.75,Form("mean = %.0f",fit_mean[elem]));
    tx_sigma[elem] = new TText(0.25,0.7,Form("sigma = %.0f",fit_sigma[elem]));
    tx_fwhm[elem] = new TText(0.25,0.65,Form("FWHM %.3f %%",fit_fwhm[elem]));
    tx_rectime[elem] -> SetNDC(1);
    tx_mean[elem] -> SetNDC(1);
    tx_sigma[elem] -> SetNDC(1);
    tx_fwhm[elem] -> SetNDC(1);
    tx_rectime[elem] -> Draw("same");
    tx_mean[elem] -> Draw("same");
    tx_sigma[elem] -> Draw("same");
    tx_fwhm[elem] -> Draw("same");
  }

  TGraph *g_fitmean = new TGraph(rec_array_size,rec_array,fit_mean);
  TGraph *g_fitsigma = new TGraph(rec_array_size,rec_array,fit_sigma);
  TGraph *g_fitfwhm = new TGraph(rec_array_size,rec_array,fit_fwhm);
  TCanvas *c_fit = new TCanvas("fit","fit",110,100,1200,300);
  c_fit -> Divide(3,1,0.01,0.01);
  c_fit -> cd(1);
  g_fitmean -> SetMarkerStyle(20);
  g_fitmean -> SetTitle("fit mean;recovery time error (ns);photon sum (1/event)");
  g_fitmean -> Draw("ap");
  c_fit -> cd(2);
  g_fitsigma -> SetMarkerStyle(20);
  g_fitsigma -> SetTitle("fit sigma;recovery time error (ns);photon Std Dev");
  g_fitsigma -> Draw("ap");
  c_fit -> cd(3);
  g_fitfwhm -> SetMarkerStyle(20);
  g_fitfwhm -> SetTitle("energy resolution (FWHM);recovery time error (ns);energy resolution (FWHM)");
  g_fitfwhm -> Draw("ap");
 

  app.Run();
  
  //photonsum_tree -> Draw("phsum_0ns");

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
