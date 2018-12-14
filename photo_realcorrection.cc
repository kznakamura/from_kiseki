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
#include <TGraph.h>

const int MAXHITCH = 3000;
const int MAXSAMPLING = 100000;
const int MAXPHOTONNUM = 100000;
const int MAXCYCLE = 25;

using namespace std;

class Analizer {
public:
  Analizer(TString photonfile, string fitfilename);
  ~Analizer();
private:
  double mem_sampling_us = 0;
  int mem_MAXHITCH = MAXHITCH;
  int mem_MAXSAMPLING = MAXSAMPLING;
  int mem_MAXPHOTONNUM = MAXPHOTONNUM;
  int mem_MAXCYCLE = MAXCYCLE;

  int mem_wav_ch_num = 0;
  double mem_wav_ch_x[MAXHITCH] = {};
  double mem_wav_ch_y[MAXHITCH] = {};
  int mem_wav_photo_num[MAXHITCH][MAXSAMPLING] = {};
  int mem_wav_start_clk = 0;
  int mem_wav_end_clk = 0;

  int mem_cycle = 0;
  double mem_onepar_nref_mean[MAXCYCLE] = {};
  double mem_onepar_residual[MAXCYCLE] = {};
  double mem_twopar_nref_mean[MAXCYCLE] = {};
  double mem_twopar_residual[MAXCYCLE] = {};

  double mem_wav_oneparcor_num[MAXHITCH][MAXSAMPLING] = {};
  double mem_wav_twoparcor_num[MAXHITCH][MAXSAMPLING] = {};
  double mem_hitwav_oneparcor_num[500*MAXSAMPLING] = {};
  double mem_hitwav_twoparcor_num[500*MAXSAMPLING] = {};
  int mem_hitwav_clk_len = 0;
  double mem_hitwav_oneparcor_sum = 0.;
  double mem_hitwav_twoparcor_sum = 0.;

  int mem_entry_num = 0;
};

int main(int argc, char* argv[]){
  if(argc!=3){
    cout << "usage : ./photo_realcorrection *****__phnum ../fitdata/****.root" << endl;
    return 1;
  }
  string photonfile = string(argv[1]);
  string fitfilepath = string(argv[2]);
  Analizer *anal = new Analizer(photonfile, fitfilepath);
  delete anal;
  return 0;
}

Analizer::Analizer(TString photonfile, string fitfilepath){
  TString ifilename = photonfile + ".root";
  TString fitfilename = fitfilepath.substr(11);
  cout << "Input photon file: " << ifilename << ", Input fit file: " << fitfilename << endl;
  TFile *ifile1 = new TFile(ifilename);
  TFile *ifile_fit = new TFile(fitfilepath.c_str());
  
  TString ofilename = photonfile + "__realcor.root";
  TFile *ofile = new TFile(ofilename,"recreate");
  
  TTree *particle_header = ((TTree*)ifile1 -> Get("particle_header")) -> CloneTree();
  TTree *electron_header = ((TTree*)ifile1 -> Get("electron_header")) -> CloneTree();
  TTree *waveform_header = ((TTree*)ifile1 -> Get("waveform_header")) -> CloneTree();
  TTree *wav2ana_header  = ((TTree*)ifile1 -> Get("wav2ana_header"))  -> CloneTree();
  TTree *photonnum_header= ((TTree*)ifile1 -> Get("photonnum_header")) -> CloneTree();
  TTree *photonnum_tree  = ((TTree*)ifile1 -> Get("photonnum_tree"));
  TTree *onepar_tree = ((TTree*)ifile_fit -> Get("onepar_tree"));
  TTree *twopar_tree = ((TTree*)ifile_fit -> Get("twopar_tree"));

  particle_header -> Write();
  electron_header -> Write();
  waveform_header -> Write();
  wav2ana_header  -> Write();
  photonnum_header-> Write();
    
  photonnum_header -> SetBranchAddress("sampling_us", &mem_sampling_us);
  photonnum_header -> GetEntry(0);
  photonnum_tree   -> SetBranchAddress("wav_ch_num", &mem_wav_ch_num);
  photonnum_tree   -> SetBranchAddress("wav_ch_x", mem_wav_ch_x);
  photonnum_tree   -> SetBranchAddress("wav_ch_y", mem_wav_ch_y);
  photonnum_tree   -> SetBranchAddress("wav_photo_num", mem_wav_photo_num);
  photonnum_tree   -> SetBranchAddress("wav_start_clk", &mem_wav_start_clk);
  photonnum_tree   -> SetBranchAddress("wav_end_clk", &mem_wav_end_clk);

  onepar_tree -> SetBranchAddress("cycle", &mem_cycle);
  onepar_tree -> SetBranchAddress("nref_mean", &mem_onepar_nref_mean);
  onepar_tree -> SetBranchAddress("residual", &mem_onepar_residual);
  twopar_tree -> SetBranchAddress("nref_mean", &mem_twopar_nref_mean);
  twopar_tree -> SetBranchAddress("residual", &mem_twopar_residual);
  onepar_tree -> GetEntry(0);
  twopar_tree -> GetEntry(0);

  TGraph *g_onepar_resi = new TGraph(mem_cycle, mem_onepar_nref_mean, mem_onepar_residual);
  TGraph *g_twopar_resi = new TGraph(mem_cycle, mem_twopar_nref_mean, mem_twopar_residual);

  TTree *realcor_header = new TTree("realcor_header", "realcor_header");
  realcor_header -> Branch("fitfilename", &fitfilename);
  realcor_header -> Branch("sampling_us", &mem_sampling_us);
  realcor_header -> Branch("MAXHITCH", &mem_MAXHITCH);
  realcor_header -> Branch("MAXSAMPLING", &mem_MAXSAMPLING);
  realcor_header -> Branch("MAXPHOTONNUMM", &mem_MAXPHOTONNUM);
  realcor_header -> Branch("MAXCYCLE", &mem_MAXCYCLE);
  realcor_header -> Fill();
  realcor_header  -> Write();
  
  TTree *realcor_tree = new TTree("realcor_tree", "realcor_tree");
  realcor_tree -> Branch("wav_ch_num", &mem_wav_ch_num, "wav_ch_num/I");
  realcor_tree -> Branch("wav_ch_x", mem_wav_ch_x, "wav_ch_x[wav_ch_num]/D");
  realcor_tree -> Branch("wav_ch_y", mem_wav_ch_y, "wav_ch_y[wav_ch_num]/D");
  realcor_tree -> Branch("wav_oneparcor_num", mem_wav_oneparcor_num, Form("wav_oneparcor_num[wav_ch_num][%d]/D",MAXSAMPLING));
  realcor_tree -> Branch("wav_twoparcor_num", mem_wav_twoparcor_num, Form("wav_twoparcor_num[wav_ch_num][%d]/D",MAXSAMPLING));
  realcor_tree -> Branch("wav_start_clk", &mem_wav_start_clk);
  realcor_tree -> Branch("wav_end_clk", &mem_wav_end_clk);
  realcor_tree -> Branch("hitwav_clk_len", &mem_hitwav_clk_len, "hitwav_clk_len/I");
  realcor_tree -> Branch("hitwav_oneparcor_num", mem_hitwav_oneparcor_num, "hitwav_oneparcor_num[hitwav_clk_len]/D");
  realcor_tree -> Branch("hitwav_twoparcor_num", mem_hitwav_twoparcor_num, "hitwav_twoparcor_num[hitwav_clk_len]/D");
  realcor_tree -> Branch("hitwav_oneparcor_sum", &mem_hitwav_oneparcor_sum);
  realcor_tree -> Branch("hitwav_twoparcor_sum", &mem_hitwav_twoparcor_sum);
  
  mem_entry_num = photonnum_tree -> GetEntries();
  cout << "#### Analysis start!! (" << mem_entry_num << " events) ####" << endl;
    
  for(int ev=0; ev<mem_entry_num; ev++){
    if(ev != 0 && ev%10 == 0) cout << "#### EVENT " << ev << " / " << mem_entry_num << " ####" << endl;
    photonnum_tree -> GetEntry(ev);
    int wav_ch_num = photonnum_tree -> GetLeaf("wav_ch_num") -> GetValue(0);
    mem_hitwav_clk_len = 0;
    mem_hitwav_oneparcor_sum = 0;
    mem_hitwav_twoparcor_sum = 0;
    for(int ch=0; ch<wav_ch_num; ch++){
      for(int smp=0; smp<MAXSAMPLING; smp++){
	mem_wav_oneparcor_num[ch][smp] = mem_wav_photo_num[ch][smp] * g_onepar_resi->Eval(mem_wav_photo_num[ch][smp]);
	mem_wav_twoparcor_num[ch][smp] = mem_wav_photo_num[ch][smp] * g_twopar_resi->Eval(mem_wav_photo_num[ch][smp]);
	if(mem_wav_photo_num[ch][smp] > MAXPHOTONNUM) cout << "MAXPHOTONNUM is too small" << endl;
	if(mem_wav_photo_num[ch][smp] > 0){
          mem_hitwav_oneparcor_num[mem_hitwav_clk_len] = mem_wav_oneparcor_num[ch][smp];
	  mem_hitwav_twoparcor_num[mem_hitwav_clk_len] = mem_wav_twoparcor_num[ch][smp];
          mem_hitwav_clk_len++;
        }
      }
    }
    mem_hitwav_oneparcor_sum = accumulate(mem_hitwav_oneparcor_num, mem_hitwav_oneparcor_num + mem_hitwav_clk_len, 0.0);
    mem_hitwav_twoparcor_sum = accumulate(mem_hitwav_twoparcor_num, mem_hitwav_twoparcor_num + mem_hitwav_clk_len, 0.0);
    realcor_tree -> Fill();
  }
  realcor_tree    -> Write();
    
  ifile1 -> Close();
  ifile_fit -> Close();
  ofile -> Close(); 
}

Analizer::~Analizer(){
  cout << "Analizer is deleted" << endl;
}

