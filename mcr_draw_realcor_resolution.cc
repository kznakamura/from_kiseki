void mcr_draw_realcor_resolution(string draw_file="dummy"){
  if(draw_file=="dummy"){
    cout << "usage: root 'mcr_draw_realcor_resolution.cc(\"input_filename\")'" << endl;
    return 1;
  }
  
  TFile *ifile = new TFile(draw_file.c_str());
  TTree *realcor_tree = ((TTree*)ifile -> Get("realcor_tree"));
  
  TCanvas *c_onepar = new TCanvas("onepar", "onepar", 100, 100, 800, 600);
  TCanvas *c_twopar = new TCanvas("twopar", "twopar", 100, 400, 800, 600);

  //**** analysis parameters ***//
  int bin_min = 600e3, bin_max = 700e3, bin_num = 1000;
  int fit_width = 5000;
    
  c_onepar -> cd();
  realcor_tree -> Draw(Form("hitwav_oneparcor_sum>>h1(%d, %d, %d)",bin_num, bin_min, bin_max),"","goff");
  TH1D *h_onepar = (TH1D*)gROOT -> FindObject("h1");
  h_onepar -> SetTitle("one parameter;clock sum;# of entries");
  double onepar_peak = bin_min + h_onepar->GetMaximumBin()*(bin_max - bin_min)/bin_num;
  TF1 *f_onepar = new TF1("f_onepar","gaus",bin_min, bin_max);
  h_onepar -> Fit("f_onepar","","",onepar_peak-fit_width, onepar_peak+fit_width);
  double onepar_mean = f_onepar -> GetParameter(1);
  double onepar_sigma = f_onepar -> GetParameter(2);
  double onepar_fwhm = onepar_sigma*2.35/onepar_mean*100;
  TText *tx_onepar_mean = new TText(0.25,0.75,Form("mean = %.0f",onepar_mean));
  TText *tx_onepar_sigma = new TText(0.25,0.7,Form("sigma = %.0f",onepar_sigma));
  TText *tx_onepar_fwhm = new TText(0.25,0.65,Form("FWHM %.3f %%",onepar_fwhm));
  tx_onepar_mean -> SetNDC(1);
  tx_onepar_sigma -> SetNDC(1);
  tx_onepar_fwhm -> SetNDC(1);  
  tx_onepar_mean -> Draw();
  tx_onepar_sigma -> Draw();
  tx_onepar_fwhm -> Draw();
 
  c_twopar -> cd();
  realcor_tree -> Draw(Form("hitwav_twoparcor_sum>>h2(%d, %d, %d)",bin_num, bin_min, bin_max),"","goff");
  TH1D *h_twopar = (TH1D*)gROOT -> FindObject("h2");
  h_twopar -> SetTitle("two parameter;clock sum;# of entries");
  double twopar_peak = bin_min + h_twopar->GetMaximumBin()*(bin_max - bin_min)/bin_num;
  TF1 *f_twopar = new TF1("f_twopar","gaus",bin_min, bin_max);
  h_twopar -> Fit("f_twopar","","",twopar_peak-fit_width, twopar_peak+fit_width);
  double twopar_mean = f_twopar -> GetParameter(1);
  double twopar_sigma = f_twopar -> GetParameter(2);
  double twopar_fwhm = twopar_sigma*2.35/twopar_mean*100;
  TText *tx_twopar_mean = new TText(0.25,0.75,Form("mean = %.0f",twopar_mean));
  TText *tx_twopar_sigma = new TText(0.25,0.7,Form("sigma = %.0f",twopar_sigma));
  TText *tx_twopar_fwhm = new TText(0.25,0.65,Form("FWHM %.3f %%",twopar_fwhm));
  tx_twopar_mean -> SetNDC(1);
  tx_twopar_sigma -> SetNDC(1);
  tx_twopar_fwhm -> SetNDC(1);  
  tx_twopar_mean -> Draw();
  tx_twopar_sigma -> Draw();
  tx_twopar_fwhm -> Draw();
}
