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

using namespace std;

class Analizer {
public:
  Analizer(string filename);
  ~Analizer();
private:
  int mem_MAXHITCH = MAXHITCH;
  int mem_MAXSAMPLING = MAXSAMPLING;
  double mem_EL_gain = 6.0; 
  double mem_sampling_us = 0.0;
  int mem_fidcut_condition = 0;

  int mem_wav_ch_num = 0;
  double mem_wav_ch_x[MAXHITCH] = {};
  double mem_wav_ch_y[MAXHITCH] = {};
  int mem_wav_ele_num[MAXHITCH][MAXSAMPLING] = {};
  int mem_wav_start_clk = 0;
  int mem_wav_end_clk = 0;

  int mem_wav_photo_num[MAXHITCH][MAXSAMPLING] = {};
  int mem_hitwav_photo_num[500*MAXSAMPLING] = {};
  int mem_hitwav_clk_len = 0;
  int mem_hitwav_photo_sum = 0;

  int mem_fid[10] = {};
  int mem_entry_num = 0;
  int mem_skipped_entry_num = 0;
};

int main(int argc, char* argv[]){
  if(argc!=2){
    cout << "usage : ./ele_to_photon (inputfile name except for extension)" << endl;
    return 1;
  }
  string filename = string(argv[1]);
  Analizer *anal = new Analizer(filename);
  delete anal;
  return 0;
}

Analizer::Analizer(string filename){
  string ifilename1 = filename + ".root";
  string ifilename2 = filename + "__05.1wav2ana.root";
  cout << "Input file 1: " << ifilename1 << endl;
  cout << "Input file 2: " << ifilename2 << endl;
 
  TFile *ifile1 = new TFile(ifilename1.c_str());
  TFile *ifile2 = new TFile(ifilename2.c_str());
  
  string ofilename = filename + "__phnum.root";
  TFile *ofile = new TFile(ofilename.c_str(),"recreate");
  
  TTree *particle_header = ((TTree*)ifile1 -> Get("particle_header")) -> CloneTree();
  TTree *electron_header = ((TTree*)ifile1 -> Get("electron_header")) -> CloneTree();
  TTree *waveform_header = ((TTree*)ifile1 -> Get("waveform_header")) -> CloneTree();
  TTree *waveform_tree   = ((TTree*)ifile1 -> Get("waveform_tree"));
  TTree *wav2ana_header  = ((TTree*)ifile2 -> Get("wav2ana_header"))  -> CloneTree();
  TTree *wav2ana_tree    = ((TTree*)ifile2 -> Get("wav2ana_tree"));
  particle_header -> Write();
  electron_header -> Write();
  waveform_header -> Write();
  wav2ana_header  -> Write();
  
  waveform_header -> SetBranchAddress("sampling_us", &mem_sampling_us);
  waveform_header -> GetEntry(0);
  waveform_tree   -> SetBranchAddress("wav_ch_num", &mem_wav_ch_num);
  waveform_tree   -> SetBranchAddress("wav_ch_x", mem_wav_ch_x);
  waveform_tree   -> SetBranchAddress("wav_ch_y", mem_wav_ch_y);
  waveform_tree   -> SetBranchAddress("wav_ele_num", mem_wav_ele_num);
  waveform_tree   -> SetBranchAddress("wav_start_clk", &mem_wav_start_clk);
  waveform_tree   -> SetBranchAddress("wav_end_clk", &mem_wav_end_clk);
  wav2ana_tree    -> SetBranchAddress("fid", mem_fid);

  TTree *photonnum_header = new TTree("photonnum_header", "photonnum_header");
  photonnum_header -> Branch("EL_gain", &mem_EL_gain);
  photonnum_header -> Branch("sampling_us", &mem_sampling_us);
  photonnum_header -> Branch("fidcut_condition", &mem_fidcut_condition);
  photonnum_header -> Branch("MAXHITCH", &mem_MAXHITCH);
  photonnum_header -> Branch("MAXSAMPLING", &mem_MAXSAMPLING);
  photonnum_header -> Fill();
  photonnum_header-> Write();
  
  TTree *photonnum_tree = new TTree("photonnum_tree", "photonnum_tree");
  photonnum_tree -> Branch("wav_ch_num", &mem_wav_ch_num, "wav_ch_num/I");
  photonnum_tree -> Branch("wav_ch_x", mem_wav_ch_x, "mem_wav_ch_x[wav_ch_num]/D");
  photonnum_tree -> Branch("wav_ch_y", mem_wav_ch_y, "mem_wav_ch_y[wav_ch_num]/D");
  photonnum_tree -> Branch("wav_photo_num", mem_wav_photo_num, Form("mem_wav_photo_num[wav_ch_num][%d]/I",MAXSAMPLING));
  photonnum_tree -> Branch("wav_start_clk", &mem_wav_start_clk);
  photonnum_tree -> Branch("wav_end_clk", &mem_wav_end_clk);
  photonnum_tree -> Branch("hitwav_clk_len", &mem_hitwav_clk_len, "hitwav_clk_len/I");
  photonnum_tree -> Branch("hitwav_photo_num", mem_hitwav_photo_num, "mem_hitwav_photo_num[hitwav_clk_len]/I");  
  photonnum_tree -> Branch("hitwav_photo_sum", &mem_hitwav_photo_sum);  

  mem_entry_num = waveform_tree -> GetEntries();
  cout << "#### Analysis start!! (" << mem_entry_num << " events) ####" << endl;
  gRandom -> SetSeed(time(NULL));
  
  for(int ev=0; ev<mem_entry_num; ev++){
    if(ev != 0 && ev%10 == 0) cout << "#### EVENT " << ev << " / " << mem_entry_num << " ####" << endl;
    wav2ana_tree  -> GetEntry(ev);
    if(mem_fid[mem_fidcut_condition] != 0){
      //cout << "#### Skip event " << ev << " ####" << endl;
      mem_skipped_entry_num+=1;
      continue;
    }
    waveform_tree -> GetEntry(ev);
    int wav_ch_num = waveform_tree -> GetLeaf("wav_ch_num") -> GetValue(0);
    mem_hitwav_clk_len = 0;
    mem_hitwav_photo_sum = 0;
    for(int ch=0; ch<wav_ch_num; ch++){
      for(int smp=0; smp<MAXSAMPLING; smp++){
	mem_wav_photo_num[ch][smp] =  gRandom->Poisson(mem_wav_ele_num[ch][smp]*mem_EL_gain);
	if(mem_wav_photo_num[ch][smp] > 0){
	  mem_hitwav_photo_num[mem_hitwav_clk_len] = mem_wav_photo_num[ch][smp];
	  mem_hitwav_clk_len++;
	}
      }
    }
    mem_hitwav_photo_sum = accumulate(mem_hitwav_photo_num, mem_hitwav_photo_num + mem_hitwav_clk_len, 0);
    photonnum_tree -> Fill();
  }

  cout << "#### Analysis finished (" << mem_skipped_entry_num << " / "<< mem_entry_num << " was skipped) ####" << endl;

  TTree *photonnum_tailer = new TTree("photonnum_tailer", "photonnum_tailer");
  photonnum_tailer -> Branch("analyzed_entries",&mem_entry_num);
  photonnum_tailer -> Branch("skipped_entries", &mem_skipped_entry_num);  
  photonnum_tailer -> Fill();
  
  photonnum_tree  -> Write();
  photonnum_tailer-> Write();
    
  ifile1 -> Close();
  ifile2 -> Close();

  ofile -> Close(); 
}

Analizer::~Analizer(){
  cout << "Analizer is deleted" << endl;
}
