void mcr_draw_photonum_resolution(string draw_file="dummy"){
  if(draw_file=="dummy"){
    cout << "usage: root 'mcr_draw_photonum_resolution.cc(\"input_filename\")'" << endl;
    return 1;
  }
  
  TFile *ifile = new TFile(draw_file.c_str());
  TTree *photonnum_tree = ((TTree*)ifile -> Get("photonnum_tree"));
  
  TCanvas *c_photonsum = new TCanvas("photonsum", "photonsum", 100, 100, 800, 600);

  //**** analysis parameters ***//
  int bin_min = 600e3, bin_max = 700e3, bin_num = 1000;
  int fit_width = 5000;
    
  c_photonsum -> cd();
  photonnum_tree -> Draw(Form("hitwav_photo_sum>>h1(%d, %d, %d)",bin_num, bin_min, bin_max),"","goff");
  TH1D *h_photosum = (TH1D*)gROOT -> FindObject("h1");
  double photosum_peak = bin_min + h_photosum->GetMaximumBin()*(bin_max - bin_min)/bin_num;
  TF1 *f_photosum = new TF1("f_photosum","gaus",bin_min, bin_max);
  h_photosum -> Fit("f_photosum","","",photosum_peak-fit_width, photosum_peak+fit_width);
  double photosum_mean = f_photosum -> GetParameter(1);
  double photosum_sigma = f_photosum -> GetParameter(2);
  double photosum_fwhm = photosum_sigma*2.35/photosum_mean*100;
  TText *tx_photosum_mean = new TText(0.25,0.75,Form("mean = %.0f",photosum_mean));
  TText *tx_photosum_sigma = new TText(0.25,0.7,Form("sigma = %.0f",photosum_sigma));
  TText *tx_photosum_fwhm = new TText(0.25,0.65,Form("FWHM %.3f %%",photosum_fwhm));
  tx_photosum_mean -> SetNDC(1);
  tx_photosum_sigma -> SetNDC(1);
  tx_photosum_fwhm -> SetNDC(1);  
  tx_photosum_mean -> Draw();
  tx_photosum_sigma -> Draw();
  tx_photosum_fwhm -> Draw();
}
