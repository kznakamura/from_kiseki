const int delta = 1000; //ns unit
const int Npix = 3600;
const int Nref_max = 40000;

double isK(double tau){return tau/(delta*Npix);}

void residual(){
  const double rec_array[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 20.0, 40.0, 60.0};
  const size_t rec_array_size = end(rec_array) - begin(rec_array);
  
  const double tau1 = 53.05, tau2 = 527.02, tau = 66.053;
  const double alpha = 0.889575, beta = 0.14244;
  
  //  const double coeff_a = 0.05, coeff_b = 2800., coeff_c = 2.*TMath::Pi()/10000., coeff_d = 1000.;
  
  TF1 *f_resi[rec_array_size];
  TLine *tl_base = new TLine(0,1,Nref_max,1);
  tl_base -> SetLineStyle(2);
  for(size_t elem=0; elem<rec_array_size; elem++){
    f_resi[elem] = new TF1(Form("resi_%dns",(int)rec_array[elem]),"1/(1+[0]*x)",0,Nref_max);
    f_resi[elem] -> SetTitle(Form("residual %dns",(int)rec_array[elem]));
    f_resi[elem] -> SetParameter(0,isK(rec_array[elem]));
    f_resi[elem] -> SetMaximum(1.02);
    f_resi[elem] -> SetMinimum(0.98); 
  }
  TCanvas *c_resi = new TCanvas("residual","residual",100,100,1250,500);
  c_resi -> Divide(5,2,0.01,0.01);
  for(size_t elem=0; elem<rec_array_size; elem++){
    c_resi -> cd(elem+1);
    f_resi[elem] -> Draw();
    tl_base -> Draw();
  }

  TF1 *f_resi_two = new TF1("resi_two", "([0]/(1+[2]*x)+[1]/(1+[3]*x))*(1+[4]*x)",0,Nref_max);
  f_resi_two -> SetParameters(alpha,beta,isK(tau1),isK(tau2),isK(tau));
  f_resi_two -> SetTitle("two residual");
  f_resi_two -> SetMaximum(1.02);
  f_resi_two -> SetMinimum(0.98);
  TText *tx_tau_two = new TText(0.25,0.8,Form("tau1 = %.1f ns, tau2 = %.0f ns", tau1, tau2));
  TText *tx_albe_two = new TText(0.25,0.75,Form("(alpha = %.2f, beta = %.2f)", alpha, beta));
  tx_tau_two -> SetNDC(1);
  tx_albe_two -> SetNDC(1);				  
  TCanvas *c_resi_two = new TCanvas("residual_two","residual_two",1350,100,250,250);
  c_resi_two -> cd();
  f_resi_two -> Draw();
  tx_tau_two -> Draw();
  tx_albe_two -> Draw();
  tl_base -> Draw();


  TLine *tl_base_cor = new TLine(0,1,Nref_max/4,1);
  tl_base_cor -> SetLineStyle(2);

  TF1 *f_cor[rec_array_size];
  for(size_t elem=0; elem<rec_array_size; elem++){
    f_cor[elem] = new TF1(Form("cor_%dns",(int)rec_array[elem]),"1/(1+[0]*x)",0,Nref_max/4);
    f_cor[elem] -> SetTitle(Form("residual %dns",(int)rec_array[elem]));
    f_cor[elem] -> SetParameter(0,isK(rec_array[elem]));
    f_cor[elem] -> SetMaximum(1.04);
    f_cor[elem] -> SetMinimum(0.92); 
  }
  TCanvas *c_cor = new TCanvas("correction","correction",100,700,1250,500);
  c_cor -> Divide(5,2,0.01,0.01);
  for(size_t elem=0; elem<rec_array_size; elem++){
    c_cor -> cd(elem+1);
    f_cor[elem] -> Draw();
    tl_base_cor -> Draw();
  }

  TF1 *f_cor_two = new TF1("cor_two", "([0]+[1]*x)/((1+[2]*x)*(1+[3]*x)-[4]*x*([0]+[1]*x))",0,Nref_max/4);
  f_cor_two -> SetParameters(alpha+beta, alpha*isK(tau2)+beta*isK(tau1),isK(tau1),isK(tau2),isK(tau));
  f_cor_two -> SetTitle("two correction");
  f_cor_two -> SetMaximum(1.04);
  f_cor_two -> SetMinimum(0.92);
  TCanvas *c_cor_two = new TCanvas("cor_two","cor_two",1350,700,250,250);
  c_cor_two -> cd();
  f_cor_two -> Draw();
  tx_tau_two -> Draw();
  tx_albe_two -> Draw();
  tl_base_cor -> Draw();

  TF1 *f_resp_linear = new TF1("resp_linear", "x", 0, Nref_max);
  f_resp_linear -> SetTitle("responce");
  TF1 *f_resp_1par = new TF1("resp_1par", "x/(1+[0]*x)", 0, Nref_max);
  f_resp_1par -> SetParameter(0,isK(tau));
  TF1 *f_resp_2par = new TF1("resp_2par", "[0]*x/(1+[2]*x)+[1]*x/(1+[3]*x)", 0, Nref_max);
  f_resp_2par -> SetParameters(alpha, beta, isK(tau1), isK(tau2));
  f_resp_linear -> SetLineStyle(2);
  f_resp_linear -> SetLineColor(1);
  f_resp_2par -> SetLineColor(3);
  TLegend *leg = new TLegend(0.2,0.7,0.4,0.85);
  leg -> AddEntry(f_resp_linear,"linear","l");
  leg -> AddEntry(f_resp_1par,"1par","l");
  leg -> AddEntry(f_resp_2par,"2par","l");

  TF1 *f_dif_2par = new TF1("dif_2par", "[0]*x/(1+[2]*x)+[1]*x/(1+[3]*x) - x", 0, Nref_max);
  f_dif_2par -> SetTitle("difference");
  f_dif_2par -> SetParameters(alpha, beta, isK(tau1), isK(tau2));
  TF1 *f_base = new TF1("base","0", 0, Nref_max);
  f_base -> SetLineStyle(2);
  f_base -> SetLineColor(1);
  TLegend *leg_dif = new TLegend(0.5,0.75,0.8,0.85);
  leg_dif -> AddEntry(f_dif_2par, "2par - linear","l"); 

  TCanvas *c_resp = new TCanvas("resp", "resp", 1600, 500, 800, 800);
  c_resp -> Divide(2,2,0.01,0.01);
  c_resp -> cd(1);
  f_resp_linear -> Draw();
  f_resp_1par -> Draw("same");
  f_resp_2par -> Draw("same");
  leg -> Draw();
  c_resp -> cd(2) -> DrawFrame(0,0,2000,2000);
  f_resp_linear -> Draw("same");
  f_resp_1par -> Draw("same");
  f_resp_2par -> Draw("same");
  leg -> Draw("same");
  c_resp -> cd(3);
  f_dif_2par -> Draw();
  f_base -> Draw("same");
  leg_dif -> Draw();
  c_resp -> cd(4) -> DrawFrame(0,-60,2000,20);
  f_dif_2par -> Draw("same");
  f_base -> Draw("same");
  leg_dif -> Draw("same");
}


