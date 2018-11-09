#include <iostream>
#include <string>
#include <TApplication.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLeaf.h>
#include <TRandom.h>
#include <vector>
#include <numeric>
#include <cmath>
#include <iomanip>
#include <TPaveText.h>
#include <TGraph.h>
#include <TStyle.h>

//const int MAXSAMPLING = 100000;
const int MAXAPERTURELEN = 20;

using namespace std;

class Analizer {
public:
  Analizer(string filename);
  ~Analizer();
private:
  double mem_aperture_array[MAXAPERTURELEN] = {};
  int mem_aperture_array_len = 0;
  double mem_tau1 = 0., mem_tau2 = 0.; //ns unit
  double mem_alpha = 0., mem_beta = 0.;
  int mem_Npix = 0;
  double mem_sampling_us = 0;
  //int mem_MAXSAMPLING = MAXSAMPLING;
  
  //  int mem_hitwav_clk_len = 0;
  // double mem_hitwav_satuph_num[MAXSAMPLING] = {};
  
  //  double mem_hitwav_ap_num[MAXSAMPLING] = {};
  double mem_hitwav_sum[MAXAPERTURELEN] = {};
 
  int mem_entry_num = 0;
  //double isK(double tau, double aperture){return tau/(mem_sampling_us*1000*mem_Npix*aperture);}
  /*double get2corFval(double photo_num, double aperture){return 
      ( -( (mem_alpha+mem_beta) - photo_num*(isK(mem_tau1, aperture)+isK(mem_tau2, aperture)) ) + sqrt( pow( (mem_alpha+mem_beta) - photo_num*( isK(mem_tau1, aperture)+isK(mem_tau2, aperture) ), 2.0) + 4*( ( mem_alpha*isK(mem_tau2, aperture) + mem_beta*isK(mem_tau1, aperture) ) - photo_num*isK(mem_tau1, aperture)*isK(mem_tau2, aperture) )*photo_num ) )
      /
      (2*( ( mem_alpha*isK(mem_tau2, aperture) + mem_beta*isK(mem_tau1, aperture) ) - photo_num*isK(mem_tau1, aperture)*isK(mem_tau2, aperture) ) );
      }*/
  double Peak(TH1D *h);
};

int main(int argc, char* argv[]){
  if(argc!=2){
    cout << "usage : ./photo_aperture (input *__phnum__2parsatu__aperture.root)" << endl;
    return 1;
  }
  
  string filename = string(argv[1]);
  Analizer *anal = new Analizer(filename);
  delete anal;
  return 0;
}

Analizer::Analizer(string filename){
  TApplication app("app",0,0,0,0);
  string ifilename1 = filename;
  cout << "Input file 1: " << ifilename1 << endl;
  TFile *ifile1 = new TFile(ifilename1.c_str());
  
  // string ofilename = filename + "__aperture.root";
  //TFile *ofile = new TFile(ofilename.c_str(),"recreate");
  
  //  TTree *particle_header = ((TTree*)ifile1 -> Get("particle_header")) -> CloneTree();
  //TTree *electron_header = ((TTree*)ifile1 -> Get("electron_header")) -> CloneTree();
  //TTree *waveform_header = ((TTree*)ifile1 -> Get("waveform_header")) -> CloneTree();
  //TTree *wav2ana_header  = ((TTree*)ifile1 -> Get("wav2ana_header"))  -> CloneTree();
  //TTree *photonnum_header= ((TTree*)ifile1 -> Get("photonnum_header")) -> CloneTree();
  TTree *aperture_header  = ((TTree*)ifile1 -> Get("aperture_header"));
  TTree *aperture_tree    = ((TTree*)ifile1 -> Get("aperture_tree"));
  //particle_header -> Write();
  //electron_header -> Write();
  //waveform_header -> Write();
  //wav2ana_header  -> Write();
  //photonnum_header-> Write();
  //ph2satu_header  -> Write();
   
  //  aperture_header -> SetBranchAddress("tau1", &mem_tau1);
  //aperture_header -> SetBranchAddress("tau2", &mem_tau2);
  //aperture_header -> SetBranchAddress("alpha", &mem_alpha);
  //aperture_header -> SetBranchAddress("beta", &mem_beta);
  //aperture_header -> SetBranchAddress("Npix", &mem_Npix);
  //aperture_header -> SetBranchAddress("sampling_us", &mem_sampling_us);
  aperture_header -> SetBranchAddress("aperture_array", mem_aperture_array);
  aperture_header -> SetBranchAddress("aperture_array_len", &mem_aperture_array_len);
  aperture_header -> SetBranchAddress("tau1", &mem_tau1);
  aperture_header -> SetBranchAddress("tau2", &mem_tau2);
  aperture_header -> SetBranchAddress("alpha", &mem_alpha);
  aperture_header -> SetBranchAddress("beta", &mem_beta);
  aperture_header -> SetBranchAddress("Npix", &mem_Npix);
  aperture_header -> SetBranchAddress("sampling_us", &mem_sampling_us);
  
  aperture_header -> GetEntry(0);
  for(int elem=0; elem<mem_aperture_array_len; elem++){
    aperture_tree   -> SetBranchAddress(Form("hitwav_sum_ap%d",(int)(mem_aperture_array[elem]*100)), &mem_hitwav_sum[elem]);
  }

  TH1D *h_phsum[mem_aperture_array_len];
  for(int elem=0; elem<mem_aperture_array_len; elem++){
    h_phsum[elem] = new TH1D(Form("phsum_ap%d",(int)(mem_aperture_array[elem]*100)),Form("phsum_ap%d",(int)(mem_aperture_array[elem]*100)), 8000, 0, 800e3);
    h_phsum[elem] -> SetTitle(Form("aperture ratio %d %%;sum of photons;# of entries",(int)(mem_aperture_array[elem]*100)));
  }

  mem_entry_num = aperture_tree -> GetEntries();
  for(int elem=0; elem<mem_aperture_array_len; elem++){
    for(int ev=0; ev<mem_entry_num; ev++){
      aperture_tree -> GetEntry(ev);
      h_phsum[elem] -> Fill(mem_hitwav_sum[elem]);
    }
  }

  TF1 *f_phsum[mem_aperture_array_len];
  double fit_mean[mem_aperture_array_len] = {};
  double fit_sigma[mem_aperture_array_len] = {};
  double fit_fwhm[mem_aperture_array_len] = {};
  double fit_peak[mem_aperture_array_len] = {};
  TText *tx_aperture[mem_aperture_array_len];
  TText *tx_mean[mem_aperture_array_len];
  TText *tx_sigma[mem_aperture_array_len];
  TText *tx_fwhm[mem_aperture_array_len];

  for(int elem=0; elem<mem_aperture_array_len; elem++){
    f_phsum[elem] = new TF1(Form("f_phsum_ap%d",(int)(mem_aperture_array[elem]*100)), "gaus",0,800e3);
    fit_peak[elem] = Peak(h_phsum[elem]);
    h_phsum[elem] -> Fit(Form("f_phsum_ap%d",(int)(mem_aperture_array[elem]*100)),"Q0","",fit_peak[elem]-50e3, fit_peak[elem]+50e3);
    fit_mean[elem] = f_phsum[elem] -> GetParameter(1);
    fit_sigma[elem] = f_phsum[elem] -> GetParameter(2);
    fit_fwhm[elem] = fit_sigma[elem]/fit_mean[elem]*2.35*100;
    tx_aperture[elem] = new TText(0.25,0.8,Form("aperture ratio = %d %%",(int)(mem_aperture_array[elem]*100)));
    tx_mean[elem] = new TText(0.25,0.75,Form("mean = %.0f",fit_mean[elem]));
    tx_sigma[elem] = new TText(0.25,0.7,Form("sigma = %.0f",fit_sigma[elem]));
    tx_fwhm[elem] = new TText(0.25,0.65,Form("FWHM %.3f %%",fit_fwhm[elem]));
    tx_aperture[elem] -> SetNDC(1);
    tx_mean[elem] -> SetNDC(1);
    tx_sigma[elem] -> SetNDC(1);
    tx_fwhm[elem] -> SetNDC(1);
  }

  TPaveText *pt_info = new TPaveText(0.1,0.1,0.9,0.9);
  pt_info -> SetFillColor(5);
  pt_info -> AddText(Form("tau1=%.1f ns, tau2=%.0f ns",mem_tau1, mem_tau2));
  pt_info -> AddText(Form("(alpha=%.2f, beta=%.2f)",mem_alpha, mem_beta));
  pt_info -> AddText(Form("Npix: %d pixel",mem_Npix));
  pt_info -> AddText(Form("sampling: %d ns/clock", (int)(mem_sampling_us*1000)));

  TCanvas *c_phsum = new TCanvas("phsum","phsum",100,100,1250,750);
  c_phsum -> Divide(4,3,0.01,0.01);
  gStyle -> SetStatX(0.3);
  TCanvas *c_phsum_q = new TCanvas("phsum_q","phsum_q",100,400,1250,750);
  c_phsum_q -> Divide(4,3,0.01,0.01);
  TH1 *frame[mem_aperture_array_len];
  for(int elem=0; elem<mem_aperture_array_len; elem++){
    c_phsum -> cd(elem+1) -> SetLogy(1);
    h_phsum[elem] -> Draw();
    f_phsum[elem] -> Draw("same");
    c_phsum -> cd(12);
    pt_info -> Draw();
    frame[elem] = c_phsum_q -> cd(elem+1) -> DrawFrame(630e3, 0, 680e3, 300);
    h_phsum[elem] -> Draw("same");
    f_phsum[elem] -> Draw("same");
    frame[elem] -> GetXaxis() -> SetTitle("sum pf photons");
    frame[elem] -> GetYaxis() -> SetTitle("# of entries");

    tx_aperture[elem] -> Draw();
    tx_mean[elem] -> Draw();
    tx_sigma[elem] -> Draw();
    tx_fwhm[elem] -> Draw();
    c_phsum_q -> cd(12);
    pt_info -> Draw();
  }

  TGraph *g_fitmean = new TGraph(mem_aperture_array_len,mem_aperture_array,fit_mean);
  TGraph *g_fitsigma = new TGraph(mem_aperture_array_len,mem_aperture_array,fit_sigma);
  TGraph *g_fitfwhm = new TGraph(mem_aperture_array_len,mem_aperture_array,fit_fwhm);
  TCanvas *c_fit = new TCanvas("fit","fit",110,100,1200,300);
  c_fit -> Divide(3,1,0.01,0.01);
  c_fit -> cd(1);
  g_fitmean -> SetMarkerStyle(20);
  g_fitmean -> SetTitle("fit mean;aperture ratio;photon sum (1/event)");
  g_fitmean -> Draw("ap");
  c_fit -> cd(2);
  g_fitsigma -> SetMarkerStyle(20);
  g_fitsigma -> SetTitle("fit sigma;aperture ratio;photon Std Dev");
  g_fitsigma -> Draw("ap");
  c_fit -> cd(3);
  g_fitfwhm -> SetMarkerStyle(20);
  g_fitfwhm -> SetTitle("energy resolution (FWHM);aperture ratio;energy resolution (FWHM)");
  g_fitfwhm -> Draw("ap");

  
  app.Run();
  //aperture_tree   -> SetBranchAddress("hitwav_satuph_num", mem_hitwav_satuph_num);
  /*
  TTree *aperture_header = new TTree("aperture_header", "aperture_header");
  aperture_header -> Branch("aperture_array_len", &mem_aperture_array_len, "aperture_array_len/I"); 
  aperture_header -> Branch("aperture_array", mem_aperture_array, "mem_aperture_array[aperture_array_len]/D");
  aperture_header -> Branch("tau1", &mem_tau1);
  aperture_header -> Branch("tau2", &mem_tau2);
  aperture_header -> Branch("alpha", &mem_alpha);
  aperture_header -> Branch("beta", &mem_beta);
  aperture_header -> Branch("Npix", &mem_Npix);
  aperture_header -> Branch("sampling_us", &mem_sampling_us);
  aperture_header -> Branch("MAXSAMPLING", &mem_MAXSAMPLING);
  aperture_header -> Fill();
  aperture_header -> Write();
  */
  /*  TTree *aperture_tree = new TTree("aperture_tree", "aperture_tree");
  for(int elem=0; elem<mem_aperture_array_len; elem++){
    aperture_tree -> Branch(Form("hitwav_sum_ap%d",(int)(mem_aperture_array[elem]*100.)), &mem_hitwav_sum[elem]);
  }
  */

  //cout << "#### Analysis start!! (" << mem_entry_num << " events) ####" << endl;
  /* 
  for(int ev=0; ev<mem_entry_num; ev++){
    //if(ev != 0 && ev%10 == 0) cout << "#### EVENT " << ev << " / " << mem_entry_num << " ####" << endl;
    aperture_tree -> GetEntry(ev);
    for(int elem=0; elem<mem_aperture_array_len; elem++){
      for(int clk=0; clk<mem_hitwav_clk_len; clk++){
	mem_hitwav_ap_num[clk] = get2corFval(mem_hitwav_satuph_num[clk], mem_aperture_array[elem]);
      }
      mem_hitwav_sum[elem] = accumulate(mem_hitwav_ap_num, mem_hitwav_ap_num + mem_hitwav_clk_len, 0.0);
    }  
    aperture_tree -> Fill();
  }
  aperture_tree -> Write();
  */
  ifile1 -> Close();
  //ofile -> Close(); 

}

Analizer::~Analizer(){
  cout << "Analizer is deleted" << endl;
}

 double Analizer::Peak(TH1D *h){
   int all_bin_num = h -> GetXaxis() -> GetNbins();
   //   int bin_content = 0;
   //int max_bin_content = 0;
   int max_bin_num = 0;
   /*
   for(int i=0; i<all_bin_num; i++){
     bin_content = h->GetBinContent(i);
     if(bin_content > max_bin_content) {
       bin_content = max_bin_content;
       max_bin_num = i;
     }
     }*/
   max_bin_num = h -> GetMaximumBin();
   int xmin = h -> GetXaxis() -> GetXmin();
   int xmax = h -> GetXaxis() -> GetXmax();
   return (double)xmin + (double)max_bin_num*(xmax - xmin)/(double)all_bin_num;
 }



