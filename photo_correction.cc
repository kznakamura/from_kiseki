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

const int MAXHITCH = 3000;
const int MAXSAMPLING = 100000;
const int MAXPHOTONNUM = 100000;

using namespace std;

class Analizer {
public:
  Analizer(string filename);
  ~Analizer();
private:
  double mem_tau = 66.05; //ns unit
  double mem_tau1 = 0., mem_tau2 = 0.; //ns unit
  double mem_alpha = 0., mem_beta = 0.;
  int mem_Npix = 0;
  double mem_sampling_us = 0;
  int mem_highlow_bound_num = 10000; // [photons/us] less than it: low, larger than it: high
  int mem_MAXHITCH = MAXHITCH;
  int mem_MAXSAMPLING = MAXSAMPLING;
  int mem_MAXPHOTONNUM = MAXPHOTONNUM;
  
  int mem_wav_ch_num = 0;
  double mem_wav_ch_x[MAXHITCH] = {};
  double mem_wav_ch_y[MAXHITCH] = {};
  double mem_wav_satuph_num[MAXHITCH][MAXSAMPLING] = {};
  int mem_wav_start_clk = 0;
  int mem_wav_end_clk = 0;

  double mem_wav_low1cor_num[MAXHITCH][MAXSAMPLING] = {};
  double mem_wav_high1cor_num[MAXHITCH][MAXSAMPLING] = {};
  double mem_wav_all1cor_num[MAXHITCH][MAXSAMPLING] = {};
  double mem_wav_all2cor_num[MAXHITCH][MAXSAMPLING] = {};
  int mem_hitwav_clk_len = 0;
  double mem_hitwav_low1cor_num[500*MAXSAMPLING] = {};
  double mem_hitwav_high1cor_num[500*MAXSAMPLING] = {};
  double mem_hitwav_all1cor_num[500*MAXSAMPLING] = {};
  double mem_hitwav_all2cor_num[500*MAXSAMPLING] = {};
  double mem_hitwav_low1cor_sum = 0;
  double mem_hitwav_high1cor_sum = 0;
  double mem_hitwav_all1cor_sum = 0;
  double mem_hitwav_all2cor_sum = 0;

  int mem_entry_num = 0;
  double isK(double tau){return tau/(mem_sampling_us*1000*mem_Npix);}
  //string make2corForm();
  double get1corFval(double photo_num){return photo_num/(1-isK(mem_tau)*photo_num);}
  double get2corFval(double photo_num){return 
      ( -( (mem_alpha+mem_beta) - photo_num*(isK(mem_tau1)+isK(mem_tau2)) ) + sqrt( pow( (mem_alpha+mem_beta) - photo_num*( isK(mem_tau1)+isK(mem_tau2) ), 2.0) + 4*( ( mem_alpha*isK(mem_tau2) + mem_beta*isK(mem_tau1) ) - photo_num*isK(mem_tau1)*isK(mem_tau2) )*photo_num ) )
      /
      (2*( ( mem_alpha*isK(mem_tau2) + mem_beta*isK(mem_tau1) ) - photo_num*isK(mem_tau1)*isK(mem_tau2) ) );
  }
};

int main(int argc, char* argv[]){
  if(argc!=2){
    cout << "usage : ./ (input *__phnum__2parsatu)" << endl;
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
  
  string ofilename = filename + "__cor.root";
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
  ph2satu_tree -> SetBranchAddress("wav_ch_num", &mem_wav_ch_num);
  ph2satu_tree -> SetBranchAddress("wav_ch_x", mem_wav_ch_x);
  ph2satu_tree -> SetBranchAddress("wav_ch_y", mem_wav_ch_y);
  ph2satu_tree -> SetBranchAddress("wav_satuph_num", mem_wav_satuph_num);
  ph2satu_tree -> SetBranchAddress("wav_start_clk", &mem_wav_start_clk);
  ph2satu_tree -> SetBranchAddress("wav_end_clk", &mem_wav_end_clk);

  TTree *phcor_header = new TTree("phcor_header", "phcor_header");
  phcor_header -> Branch("tau", &mem_tau);
  phcor_header -> Branch("tau1", &mem_tau1);
  phcor_header -> Branch("tau2", &mem_tau2);
  phcor_header -> Branch("alpha", &mem_alpha);
  phcor_header -> Branch("beta", &mem_beta);
  phcor_header -> Branch("Npix", &mem_Npix);
  phcor_header -> Branch("sampling_us", &mem_sampling_us);
  phcor_header -> Branch("highlow_bound_num", &mem_highlow_bound_num);
  phcor_header -> Branch("MAXHITCH", &mem_MAXHITCH);
  phcor_header -> Branch("MAXSAMPLING", &mem_MAXSAMPLING);
  phcor_header -> Branch("MAXPHOTONNUM", &mem_MAXPHOTONNUM);
  phcor_header -> Fill();
  phcor_header   -> Write();
  
  TTree *phcor_tree = new TTree("phcor_tree", "phcor_tree");
  phcor_tree -> Branch("wav_ch_num", &mem_wav_ch_num, "wav_ch_num/I");
  phcor_tree -> Branch("wav_ch_x", mem_wav_ch_x, "mem_wav_ch_x[wav_ch_num]/D");
  phcor_tree -> Branch("wav_ch_y", mem_wav_ch_y, "mem_wav_ch_y[wav_ch_num]/D");
  phcor_tree -> Branch("wav_low1cor_num", mem_wav_low1cor_num, Form("mem_wav_low1cor_num[wav_ch_num][%d]/D",MAXSAMPLING));
  phcor_tree -> Branch("wav_high1cor_num", mem_wav_high1cor_num, Form("mem_wav_high1cor_num[wav_ch_num][%d]/D",MAXSAMPLING));
  phcor_tree -> Branch("wav_all1cor_num", mem_wav_all1cor_num, Form("mem_wav_all1cor_num[wav_ch_num][%d]/D",MAXSAMPLING));
  phcor_tree -> Branch("wav_all2cor_num", mem_wav_all2cor_num, Form("mem_wav_all2cor_num[wav_ch_num][%d]/D",MAXSAMPLING));
  phcor_tree -> Branch("wav_start_clk", &mem_wav_start_clk);
  phcor_tree -> Branch("wav_end_clk", &mem_wav_end_clk);
  phcor_tree -> Branch("hitwav_clk_len", &mem_hitwav_clk_len, "hitwav_clk_len/I");
  phcor_tree -> Branch("hitwav_low1cor_num", mem_hitwav_low1cor_num, "mem_hitwav_low1cor_num[hitwav_clk_len]/D");
  phcor_tree -> Branch("hitwav_high1cor_num", mem_hitwav_high1cor_num, "mem_hitwav_high1cor_num[hitwav_clk_len]/D");
  phcor_tree -> Branch("hitwav_all1cor_num", mem_hitwav_all1cor_num, "mem_hitwav_all1cor_num[hitwav_clk_len]/D");
  phcor_tree -> Branch("hitwav_all2cor_num", mem_hitwav_all2cor_num, "mem_hitwav_all2cor_num[hitwav_clk_len]/D");
  phcor_tree -> Branch("hitwav_low1cor_sum", &mem_hitwav_low1cor_sum);
  phcor_tree -> Branch("hitwav_high1cor_sum", &mem_hitwav_high1cor_sum);
  phcor_tree -> Branch("hitwav_all1cor_sum", &mem_hitwav_all1cor_sum);
  phcor_tree -> Branch("hitwav_all2cor_sum", &mem_hitwav_all2cor_sum);
  
  mem_entry_num = ph2satu_tree -> GetEntries();
  cout << "#### Analysis start!! (" << mem_entry_num << " events) ####" << endl;
    
  for(int ev=0; ev<mem_entry_num; ev++){
    if(ev != 0 && ev%10 == 0) cout << "#### EVENT " << ev << " / " << mem_entry_num << " ####" << endl;
    ph2satu_tree -> GetEntry(ev);
    int wav_ch_num = ph2satu_tree -> GetLeaf("wav_ch_num") -> GetValue(0);
    mem_hitwav_clk_len = 0;
    mem_hitwav_low1cor_sum = 0;
    mem_hitwav_high1cor_sum = 0;
    mem_hitwav_all1cor_sum = 0;
    for(int ch=0; ch<wav_ch_num; ch++){
      for(int smp=0; smp<MAXSAMPLING; smp++){
	if(mem_wav_satuph_num[ch][smp] >= 0 && mem_wav_satuph_num[ch][smp] < mem_highlow_bound_num*mem_sampling_us){
	  mem_wav_low1cor_num[ch][smp] = get1corFval(mem_wav_satuph_num[ch][smp]);
	  mem_wav_high1cor_num[ch][smp] = get2corFval(mem_wav_satuph_num[ch][smp]);
	}else if(mem_wav_satuph_num[ch][smp] >= mem_highlow_bound_num*mem_sampling_us && mem_wav_satuph_num[ch][smp] < MAXPHOTONNUM){
	  mem_wav_low1cor_num[ch][smp] = get2corFval(mem_wav_satuph_num[ch][smp]);
	  mem_wav_high1cor_num[ch][smp] = get1corFval(mem_wav_satuph_num[ch][smp]);
	}else{
	  cout << "### input photon num is out of range ###" << endl;
	}
	mem_wav_all1cor_num[ch][smp] = get1corFval(mem_wav_satuph_num[ch][smp]);
	mem_wav_all2cor_num[ch][smp] = get2corFval(mem_wav_satuph_num[ch][smp]);
	if(mem_wav_satuph_num[ch][smp] > 0){
	  mem_hitwav_low1cor_num[mem_hitwav_clk_len] = mem_wav_low1cor_num[ch][smp];
	  mem_hitwav_high1cor_num[mem_hitwav_clk_len] = mem_wav_high1cor_num[ch][smp];
	  mem_hitwav_all1cor_num[mem_hitwav_clk_len] = mem_wav_all1cor_num[ch][smp];
	  mem_hitwav_all2cor_num[mem_hitwav_clk_len] = mem_wav_all2cor_num[ch][smp];
	  mem_hitwav_clk_len++;
	}
      }
    }
    mem_hitwav_low1cor_sum = accumulate(mem_hitwav_low1cor_num, mem_hitwav_low1cor_num + mem_hitwav_clk_len, 0.0);
    mem_hitwav_high1cor_sum = accumulate(mem_hitwav_high1cor_num, mem_hitwav_high1cor_num + mem_hitwav_clk_len, 0.0);
    mem_hitwav_all1cor_sum = accumulate(mem_hitwav_all1cor_num, mem_hitwav_all1cor_num + mem_hitwav_clk_len, 0.0);
    mem_hitwav_all2cor_sum = accumulate(mem_hitwav_all2cor_num, mem_hitwav_all2cor_num + mem_hitwav_clk_len, 0.0);
    phcor_tree -> Fill();
  }
  phcor_tree -> Write();
    
  ifile1 -> Close();

  ofile -> Close(); 
}

Analizer::~Analizer(){
  cout << "Analizer is deleted" << endl;
}
