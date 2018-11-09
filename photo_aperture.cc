#include <iostream>
#include <string>
#include <TApplication.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TLeaf.h>
#include <TRandom.h>
#include <vector>
#include <numeric>
#include <cmath>
#include <iomanip>

const int MAXSAMPLING = 100000;
const int MAXAPERTURELEN = 20;
double mem_aperture_array[] = {1.0, 0.99, 0.98, 0.97, 0.96, 0.95, 0.9, 0.85, 0.8, 0.75};

using namespace std;

class Analizer {
public:
  Analizer(string filename);
  ~Analizer();
private:
  int mem_aperture_array_len = (int)(end(mem_aperture_array) - begin(mem_aperture_array));;
  double mem_tau1 = 0., mem_tau2 = 0.; //ns unit
  double mem_alpha = 0., mem_beta = 0.;
  int mem_Npix = 0;
  double mem_sampling_us = 0;
  int mem_MAXSAMPLING = MAXSAMPLING;
  int mem_MAXAPERTURELEN = MAXAPERTURELEN;
  
  int mem_hitwav_clk_len = 0;
  double mem_hitwav_satuph_num[MAXSAMPLING] = {};
  
  double mem_hitwav_ap_num[MAXSAMPLING] = {};
  double mem_hitwav_sum[MAXAPERTURELEN] = {};
 
  int mem_entry_num = 0;
  double isK(double tau, double aperture){return tau/(mem_sampling_us*1000*mem_Npix*aperture);}
  double get2corFval(double photo_num, double aperture){return 
      ( -( (mem_alpha+mem_beta) - photo_num*(isK(mem_tau1, aperture)+isK(mem_tau2, aperture)) ) + sqrt( pow( (mem_alpha+mem_beta) - photo_num*( isK(mem_tau1, aperture)+isK(mem_tau2, aperture) ), 2.0) + 4*( ( mem_alpha*isK(mem_tau2, aperture) + mem_beta*isK(mem_tau1, aperture) ) - photo_num*isK(mem_tau1, aperture)*isK(mem_tau2, aperture) )*photo_num ) )
      /
      (2*( ( mem_alpha*isK(mem_tau2, aperture) + mem_beta*isK(mem_tau1, aperture) ) - photo_num*isK(mem_tau1, aperture)*isK(mem_tau2, aperture) ) );
  }
};

int main(int argc, char* argv[]){
  if(argc!=2){
    cout << "usage : ./photo_aperture (input *__phnum__2parsatu)" << endl;
    return 1;
  }
  string filename = string(argv[1]);
  Analizer *anal = new Analizer(filename);
  delete anal;
  return 0;
}

Analizer::Analizer(string filename){
  string ifilename1 = filename + ".root";
  cout << "Input file 1: " << ifilename1 << endl;
  TFile *ifile1 = new TFile(ifilename1.c_str());
  
  string ofilename = filename + "__aperture.root";
  TFile *ofile = new TFile(ofilename.c_str(),"recreate");
  
  TTree *particle_header = ((TTree*)ifile1 -> Get("particle_header")) -> CloneTree();
  TTree *electron_header = ((TTree*)ifile1 -> Get("electron_header")) -> CloneTree();
  TTree *waveform_header = ((TTree*)ifile1 -> Get("waveform_header")) -> CloneTree();
  TTree *wav2ana_header  = ((TTree*)ifile1 -> Get("wav2ana_header"))  -> CloneTree();
  TTree *photonnum_header= ((TTree*)ifile1 -> Get("photonnum_header")) -> CloneTree();
  TTree *ph2satu_header  = ((TTree*)ifile1 -> Get("ph2satu_header")) -> CloneTree();
  TTree *ph2satu_tree    = ((TTree*)ifile1 -> Get("ph2satu_tree"));
  particle_header -> Write();
  electron_header -> Write();
  waveform_header -> Write();
  wav2ana_header  -> Write();
  photonnum_header-> Write();
  ph2satu_header  -> Write();
   
  ph2satu_header -> SetBranchAddress("tau1", &mem_tau1);
  ph2satu_header -> SetBranchAddress("tau2", &mem_tau2);
  ph2satu_header -> SetBranchAddress("alpha", &mem_alpha);
  ph2satu_header -> SetBranchAddress("beta", &mem_beta);
  ph2satu_header -> SetBranchAddress("Npix", &mem_Npix);
  ph2satu_header -> SetBranchAddress("sampling_us", &mem_sampling_us);
  ph2satu_header -> GetEntry(0);
  ph2satu_tree   -> SetBranchAddress("hitwav_clk_len", &mem_hitwav_clk_len);
  ph2satu_tree   -> SetBranchAddress("hitwav_satuph_num", mem_hitwav_satuph_num);
  
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
  aperture_header -> Branch("MAXAPERTURELEN", &mem_MAXAPERTURELEN);
  aperture_header -> Fill();
  aperture_header -> Write();
  
  TTree *aperture_tree = new TTree("aperture_tree", "aperture_tree");
  for(int elem=0; elem<mem_aperture_array_len; elem++){
    aperture_tree -> Branch(Form("hitwav_sum_ap%d",(int)(mem_aperture_array[elem]*100.)), &mem_hitwav_sum[elem]);
  }
  
  mem_entry_num = ph2satu_tree -> GetEntries();
  cout << "#### Analysis start!! (" << mem_entry_num << " events) ####" << endl;
  
  for(int ev=0; ev<mem_entry_num; ev++){
    if(ev != 0 && ev%10 == 0) cout << "#### EVENT " << ev << " / " << mem_entry_num << " ####" << endl;
    ph2satu_tree -> GetEntry(ev);
    for(int elem=0; elem<mem_aperture_array_len; elem++){
      for(int clk=0; clk<mem_hitwav_clk_len; clk++){
	mem_hitwav_ap_num[clk] = get2corFval(mem_hitwav_satuph_num[clk], mem_aperture_array[elem]);
      }
      mem_hitwav_sum[elem] = accumulate(mem_hitwav_ap_num, mem_hitwav_ap_num + mem_hitwav_clk_len, 0.0);
    }  
    aperture_tree -> Fill();
  }
  aperture_tree -> Write();
    
  ifile1 -> Close();
  ofile -> Close(); 
}

Analizer::~Analizer(){
  cout << "Analizer is deleted" << endl;
}
