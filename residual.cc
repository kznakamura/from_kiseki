const int delta = 1000; //ns unit
const int Npix = 3600;
const int Nref_max = 40000;

double isK(double tau){return tau/(delta*Npix);}

void residual(){
  const double rec_array[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 20.0, 40.0, 60.0};
  const size_t rec_array_size = end(rec_array) - begin(rec_array);
  
  const double tau1 = 53.05, tau2 = 527.02, tau = 66.053;
  const double alpha = 0.889575, beta = 0.14244;
  
  const double coeff_a = 0.05, coeff_b = 2800., coeff_c = 2.*TMath::Pi()/10000., coeff_d = 1000.;
  
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
  
}


