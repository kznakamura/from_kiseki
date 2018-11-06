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

const int MAXHITCH = 3000;
const int MAXSAMPLING = 100000;
const int MAXPHOTONNUM = 100000;

using namespace std;

class Analizer {
public:
  Analizer(string filename);
  ~Analizer();
private:
  double mem_tau1 = 53.05, mem_tau2 = 527.02; //ns unit
  double mem_alpha = 0.8895, mem_beta = 0.1424;
  double mem_sampling_us = 0;
  int mem_Npix = 3600;
  int mem_MAXHITCH = MAXHITCH;
  int mem_MAXSAMPLING = MAXSAMPLING;
  int mem_MAXPHOTONNUM = MAXPHOTONNUM;

  int mem_wav_ch_num = 0;
  double mem_wav_ch_x[MAXHITCH] = {};
  double mem_wav_ch_y[MAXHITCH] = {};
  int mem_wav_photo_num[MAXHITCH][MAXSAMPLING] = {};
  int mem_wav_start_clk = 0;
  int mem_wav_end_clk = 0;

  double mem_wav_satuph_num[MAXHITCH][MAXSAMPLING] = {};
  double mem_hitwav_satuph_num[500*MAXSAMPLING] = {};
  int mem_hitwav_clk_len = 0;
  double mem_hitwav_satuph_sum = 0.;

  int mem_entry_num = 0;
  double isK(double tau){return tau/(mem_sampling_us*1000*mem_Npix);}
  double get2satuFval(double photo_num){return mem_alpha*photo_num/(1+isK(mem_tau1)*photo_num) + mem_beta*photo_num/(1+isK(mem_tau2)*photo_num);}
};

int main(int argc, char* argv[]){
  if(argc!=2){
    cout << "usage : ./ (input *__phnum)" << endl;
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
  
  string ofilename = filename + "__2parsatu.root";
  TFile *ofile = new TFile(ofilename.c_str(),"recreate");
  
  TTree *particle_header = ((TTree*)ifile1 -> Get("particle_header")) -> CloneTree();
  TTree *electron_header = ((TTree*)ifile1 -> Get("electron_header")) -> CloneTree();
  TTree *waveform_header = ((TTree*)ifile1 -> Get("waveform_header")) -> CloneTree();
  TTree *wav2ana_header  = ((TTree*)ifile1 -> Get("wav2ana_header"))  -> CloneTree();
  TTree *photonnum_header= ((TTree*)ifile1 -> Get("photonnum_header")) -> CloneTree();
  TTree *photonnum_tree  = ((TTree*)ifile1 -> Get("photonnum_tree"));
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

  TTree *ph2satu_header = new TTree("ph2satu_header", "ph2satu_header");
  ph2satu_header -> Branch("tau1", &mem_tau1);
  ph2satu_header -> Branch("tau2", &mem_tau2);
  ph2satu_header -> Branch("alpha", &mem_alpha);
  ph2satu_header -> Branch("beta", &mem_beta);
  ph2satu_header -> Branch("Npix", &mem_Npix);
  ph2satu_header -> Branch("sampling_us", &mem_sampling_us);
  ph2satu_header -> Branch("MAXHITCH", &mem_MAXHITCH);
  ph2satu_header -> Branch("MAXSAMPLING", &mem_MAXSAMPLING);
  ph2satu_header -> Branch("MAXPHOTONNUMM", &mem_MAXPHOTONNUM);
  ph2satu_header -> Fill();
  ph2satu_header  -> Write();
  
  TTree *ph2satu_tree = new TTree("ph2satu_tree", "ph2satu_tree");
  ph2satu_tree -> Branch("wav_ch_num", &mem_wav_ch_num, "wav_ch_num/I");
  ph2satu_tree -> Branch("wav_ch_x", mem_wav_ch_x, "mem_wav_ch_x[wav_ch_num]/D");
  ph2satu_tree -> Branch("wav_ch_y", mem_wav_ch_y, "mem_wav_ch_y[wav_ch_num]/D");
  ph2satu_tree -> Branch("wav_satuph_num", mem_wav_satuph_num, Form("mem_wav_satuph_num[wav_ch_num][%d]/D",MAXSAMPLING));
  ph2satu_tree -> Branch("wav_start_clk", &mem_wav_start_clk);
  ph2satu_tree -> Branch("wav_end_clk", &mem_wav_end_clk);
  ph2satu_tree -> Branch("hitwav_clk_len", &mem_hitwav_clk_len, "hitwav_clk_len/I");
  ph2satu_tree -> Branch("hitwav_satuph_num", mem_hitwav_satuph_num, "mem_hitwav_satuph_num[hitwav_clk_len]/D");
  ph2satu_tree -> Branch("hitwav_satuph_sum", &mem_hitwav_satuph_sum);
  
  mem_entry_num = photonnum_tree -> GetEntries();
  cout << "#### Analysis start!! (" << mem_entry_num << " events) ####" << endl;
    
  for(int ev=0; ev<mem_entry_num; ev++){
    if(ev != 0 && ev%10 == 0) cout << "#### EVENT " << ev << " / " << mem_entry_num << " ####" << endl;
    photonnum_tree -> GetEntry(ev);
    int wav_ch_num = photonnum_tree -> GetLeaf("wav_ch_num") -> GetValue(0);
    mem_hitwav_clk_len = 0;
    mem_hitwav_satuph_sum = 0;
    for(int ch=0; ch<wav_ch_num; ch++){
      for(int smp=0; smp<MAXSAMPLING; smp++){
	mem_wav_satuph_num[ch][smp] =  get2satuFval(mem_wav_photo_num[ch][smp]);
	if(mem_wav_photo_num[ch][smp] > MAXPHOTONNUM) cout << "MAXPHOTONNUM is too small" << endl;
	if(mem_wav_photo_num[ch][smp] > 0){
          mem_hitwav_satuph_num[mem_hitwav_clk_len] = mem_wav_satuph_num[ch][smp];
          mem_hitwav_clk_len++;
        }
      }
    }
    mem_hitwav_satuph_sum = accumulate(mem_hitwav_satuph_num, mem_hitwav_satuph_num + mem_hitwav_clk_len, 0.0);
    ph2satu_tree -> Fill();
  }
  ph2satu_tree    -> Write();
    
  ifile1 -> Close();
  ofile -> Close(); 
}

Analizer::~Analizer(){
  cout << "Analizer is deleted" << endl;
}

