//double coeff_a = 0.05, coeff_b = 2800., coeff_c = 10000., coeff_d = 1000.;
double coeff_a = 0.05, coeff_b = 1000., coeff_c = 10000., coeff_d = 1500.;

void correction(){

  TF1 *f_resi = new TF1("residual","-[0]*exp(-x/[1])*cos([2]*(x+[3]))+1", 0, 40000);
  f_resi -> SetParameters(coeff_a, coeff_b, 2.*TMath::Pi()/coeff_c, coeff_d);
  TLatex *tl = new TLatex();
  TCanvas *c_resi = new TCanvas("residual","residual",800,600);
  f_resi -> SetTitle("real residual");
  f_resi -> SetMaximum(1.2);
  f_resi -> SetMinimum(0.96);
  f_resi -> SetMaximum(1.02);
  f_resi -> Draw();
  tl -> SetTextSize(0.05);
  tl -> DrawLatex(5000,0.98,"Data/Fit = -Ae^{-x/B}cos#left[C(x+D)#right]+1");
  tl -> DrawLatex(5000, 0.975,Form("(A=%0.2f, B=%d, C=2#pi/%d, D=%d)",coeff_a, (int)coeff_b, (int)coeff_c, (int)coeff_d));

}
