#include <iostream>
#include <string>
#include <TApplication.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TLeaf.h>
#include <TRandom.h>
#include <vector>


const int MAXHITCH = 3000;
const int MAXSAMPLING = 100000;
const int MAXPHOTONNUM = 100000;
double EL_gain = 6.0; 
int Npix = 3600;
int fidcut_condition = 0;

using namespace std;

class Analizer {
public:
  Analizer(string filename);
  ~Analizer();


private:
  vector<double> rec_array{0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 20.0, 40.0, 60.0}; //ns unit
  double mem_sampling_us = 0;
  int mem_wav_ele_num[MAXHITCH][MAXSAMPLING]={};
  int PhotonSum(TF1 *f_rectime, int wav_ch_num);
  int mem_fid[10] = {};
  int mem_entry_num = 0;
  int mem_skipped_entry_num = 0;
};


int main(int argc, char* argv[]){
  if(argc!=2){
    cout << "usage : ./photonsum.cc (inputfile name except for extension)" << endl;
    return 1;
  }
  string filename = string(argv[1]);
  //  TApplication app("app", &argc, argv );

  Analizer *anal = new Analizer(filename);
  delete anal;
  //  app.Run();
  return 0;
}

Analizer::Analizer(string filename){
  string ifilename1 = filename + ".root";
  string ifilename2 = filename + "__05.1wav2ana.root";
  cout << "Input file 1: " << ifilename1 << endl;
  cout << "Input file 2: " << ifilename2 << endl;
 
  TFile *ifile1 = new TFile(ifilename1.c_str());
  TFile *ifile2 = new TFile(ifilename2.c_str());
  
  string ofilename = filename + "__phsum.root";
  TFile *ofile = new TFile(ofilename.c_str(),"recreate");
  
  TTree *particle_header = ((TTree*)ifile1 -> Get("particle_header")) -> CloneTree();
  TTree *electron_header = ((TTree*)ifile1 -> Get("electron_header")) -> CloneTree();
  TTree *waveform_header = ((TTree*)ifile1 -> Get("waveform_header")) -> CloneTree();
  TTree *waveform_tree   = ((TTree*)ifile1 -> Get("waveform_tree"));
  TTree *wav2ana_header  = ((TTree*)ifile2 -> Get("wav2ana_header"))  -> CloneTree();
  TTree *wav2ana_tree    = ((TTree*)ifile2 -> Get("wav2ana_tree"));
  
  waveform_header -> SetBranchAddress("sampling_us", &mem_sampling_us);
  waveform_header -> GetEntry(0);
  waveform_tree   -> SetBranchAddress("wav_ele_num", mem_wav_ele_num);
  wav2ana_tree    -> SetBranchAddress("fid", mem_fid);

  const size_t rec_array_size = rec_array.size();

  TTree *photonsum_header = new TTree("photonsum_header", "photonsum_header");
  TTree *photonsum_tree = new TTree("photonsum_tree", "photonsum_tree");
  //  photonsum_header -> Branch("rec_array_size", &mem_rec_array_size, "mem_rec_array_size/I");
  photonsum_header -> Branch("rectime", &rec_array);
  // photonsum_header -> Branch("rectime", rec_array, "rec_array[3]/D");
  photonsum_header -> Branch("EL_gain", &EL_gain);
  photonsum_header -> Branch("Npix", &Npix);
  photonsum_header -> Branch("sampling_us", &mem_sampling_us);
  photonsum_header -> Branch("fidcut_condition", &fidcut_condition);
  photonsum_header -> Fill();
    
  int photon_sum[rec_array_size] = {};
  for(size_t elem=0; elem<rec_array_size; elem++){
    photonsum_tree -> Branch(Form("phsum_%dns",(int)rec_array[elem]), &photon_sum[elem]);
  }
  
  TF1 *f_rec[rec_array_size];
  for(size_t elem=0; elem<rec_array_size; elem++){
    f_rec[elem] = new TF1(Form("%dns",(int)rec_array[elem]),"x/(1+[0]*x)", 0, MAXPHOTONNUM);
    f_rec[elem] -> SetParameter(0,rec_array[elem]/(mem_sampling_us*1000*Npix));
  }
  
  mem_entry_num = waveform_tree -> GetEntries();
  cout << "#### Analysis start!! (" << mem_entry_num << " events) ####" << endl;
  gRandom -> SetSeed(time(NULL));
  
  for(int ev=0; ev<mem_entry_num; ev++){
    if(ev != 0 && ev%10 == 0) cout << "#### EVENT " << ev << " / " << mem_entry_num << " ####" << endl;
    wav2ana_tree  -> GetEntry(ev);
    if(mem_fid[fidcut_condition] != 0){
      //cout << "#### Skip event " << ev << " ####" << endl;
      mem_skipped_entry_num+=1;
      continue;
    }
    waveform_tree -> GetEntry(ev);
    int wav_ch_num = waveform_tree -> GetLeaf("wav_ch_num") -> GetValue(0);

    for(size_t elem=0; elem<rec_array_size; elem++){
      //      cout << rec_array[elem] << " " << PhotonSum(f_rec[elem], wav_ch_num) << endl;      
      photon_sum[elem] = PhotonSum(f_rec[elem], wav_ch_num);
    }
    photonsum_tree -> Fill();
  }

  cout << "#### Analysis finished (" << mem_skipped_entry_num << " / "<< mem_entry_num << " was skipped) ####" << endl;

  TTree *photonsum_tailer = new TTree("photonsum_tailer", "photonsum_tailer");
  photonsum_tailer -> Branch("analyzed_entries",&mem_entry_num);
  photonsum_tailer -> Branch("skipped_entries", &mem_skipped_entry_num);
  photonsum_tailer -> Fill();
  
  particle_header -> Write();
  electron_header -> Write();
  waveform_header -> Write();
  wav2ana_header  -> Write();
  photonsum_header-> Write();
  photonsum_tree  -> Write();
  photonsum_tailer-> Write();
    
  ifile1 -> Close();
  ifile2 -> Close();

  ofile -> Close(); 
}

Analizer::~Analizer(){
  cout << "Analizer is deleted" << endl;
}

int Analizer::PhotonSum(TF1 *f_rectime, int wav_ch_num){
  double photon_sum=0.0;
  for(int ch=0; ch<wav_ch_num; ch++){
    for(int smp=0; smp<MAXSAMPLING; smp++){
      photon_sum += f_rectime->Eval(gRandom->Poisson(mem_wav_ele_num[ch][smp]*EL_gain));
      /*      if(mem_wav_ele_num[ch][smp]*random_num >300){
	cout<< mem_wav_ele_num[ch][smp]*random_num << " " << f_rectime->Eval(mem_wav_ele_num[ch][smp]*random_num) << endl;
	}*/
      if(mem_wav_ele_num[ch][smp]>MAXPHOTONNUM){
	cout << "MAXPHOTONNUM is too small" << endl;
	return 1;
      }
    }
  }
  return (int)photon_sum;
  }
