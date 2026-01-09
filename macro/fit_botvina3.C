TH2* cutTH2(const std::unique_ptr<TH2D>& hFull, TCutG* cut) {
  std::string name = hFull->GetName();
  name += "_cut";
  TH2* hCut = (TH2*) hFull->Clone(name.c_str());
  int nBinsX = hCut->GetNbinsX();
  int nBinsY = hCut->GetNbinsY();
  double x, y;

  for (int i = 1; i <= nBinsX; i++) {
    for (int j = 1; j <= nBinsY; j++) {
      x = hCut->GetXaxis()->GetBinCenter(i);
      y = hCut->GetYaxis()->GetBinCenter(j);
      if (!cut->IsInside(x, y)) {
        hCut->SetBinContent(i, j, 0.);
        hCut->SetBinError(i, j, 0.);
      }
    }
    //    hCut -> Draw ("colz");
    //    gPad -> SetLogz ();
    //    cut -> SetLineColor (kRed);
    //    cut -> Draw ("same");
  }
  return hCut;
}

double asymmGauss(double* xx, double* p) {
  double x = xx[0];
  double amp = p[0];
  double mean = p[1];
  double sig = p[2];
  double delta = p[3];
  if (x > mean) sig *= 1 + delta;
  else
    sig *= 1 - delta;
  return amp * exp(-0.5 * (x - mean) * (x - mean) / sig / sig);
}

void fit_botvina3() {
  auto start = std::chrono::system_clock::now();
  gROOT->SetBatch(true);
  gROOT->SetEscape(true);
  ROOT::EnableImplicitMT(2);

  TString inputFileName = "../input/botv_3.3agev_default.root";
  TString cutsFileName = "../input/cuts.root";
  std::unique_ptr<TFile> fIn{TFile::Open(inputFileName, "read")};
  std::unique_ptr<TFile> fCuts{TFile::Open(cutsFileName, "read")};

  Pid::Fitter tof;
  Pid::GetterTof getter;
  float xmin, xmax, ymin, ymax;

  cout << "\n\npionpos\n";
  xmin = 0.25, xmax = 12., ymin = -2., ymax = 2.;
  tof.SetChi2Max(10000);
  std::unique_ptr<TH2D> hpionpos{(TH2D*) fIn->Get("reco_vs_sim_info/h2TofM2_piplus")};
  std::unique_ptr<TH2D> hpionpos_cut{(TH2D*) cutTH2(hpionpos, (TCutG*) fCuts->Get("piplus"))};
  hpionpos_cut->Rebin2D(10, 1);
  TF1 fit_pionpos("fit_pionpos", "gaus", ymin, ymax);
  fit_pionpos.SetParNames("p0", "p1", "p2");
  fit_pionpos.SetParLimits(0, 0., 2.e6);
  fit_pionpos.SetParLimits(1, -.1, 0.03);
  fit_pionpos.SetParLimits(2, 0., 0.35);
  TF1 pionpos_0("pionpos_0", "0", xmin, xmax);
  TF1 pionpos_1("pionpos_1", "pol2", xmin, xmax);
  TF1 pionpos_2("pionpos_2", "pol5", xmin, xmax);

  Pid::ParticleFit pionpos(PidTofParticles::kPionPos);
  pionpos.SetParametrization({pionpos_0, pionpos_1, pionpos_2});
  pionpos.SetFitFunction(fit_pionpos);
  pionpos.SetRange(xmin, xmax);
  pionpos.SetIsFitted();

  tof.AddParticle(pionpos, PidTofParticles::kPionPos);
  tof.SetHisto2D(std::move(hpionpos_cut));
  tof.SetRangeX(xmin, 8.);
  tof.SetRangeY(ymin, ymax);
  tof.SetOutputFileName("pionpos.root");
  tof.Fit();
  pionpos = tof.GetParticle(0);
  tof.Clear();

  cout << "\n\nkaonpos\n";
  xmin = 0.5, xmax = 4.5, ymin = -2., ymax = 2.;
  tof.SetChi2Max(1000);
  std::unique_ptr<TH2D> hkaonpos{(TH2D*) fIn->Get("reco_vs_sim_info/h2TofM2_kplus")};
  std::unique_ptr<TH2D> hkaonpos_cut{(TH2D*) cutTH2(hkaonpos, (TCutG*) fCuts->Get("kplus"))};
  hkaonpos_cut->Rebin2D(10, 1);
  TF1 fit_kaonpos("fit_kaonpos", "gaus", ymin, ymax);
  fit_kaonpos.SetParNames("p3", "p4", "p5");
  fit_kaonpos.SetParLimits(0, 0., 1.e5);
  fit_kaonpos.SetParLimits(1, -0.05, 0.3);
  fit_kaonpos.SetParLimits(2, 0., 3.);
  TF1 kaonpos_0("kaonpos_0", "0", xmin, xmax);
  TF1 kaonpos_1("kaonpos_1", "pol2", xmin, xmax);
  TF1 kaonpos_2("kaonpos_2", "pol6", xmin, xmax);

  Pid::ParticleFit kaonpos(PidTofParticles::kKaonPos);
  kaonpos.SetParametrization({kaonpos_0, kaonpos_1, kaonpos_2});
  kaonpos.SetFitFunction(fit_kaonpos);
  kaonpos.SetRange(xmin, xmax);
  kaonpos.SetIsFitted();

  tof.AddParticle(kaonpos, PidTofParticles::kKaonPos);
  tof.SetHisto2D(std::move(hkaonpos_cut));
  tof.SetRangeX(xmin, xmax);
  tof.SetRangeY(ymin, ymax);
  tof.SetOutputFileName("kaonpos.root");
  tof.Fit();
  kaonpos = tof.GetParticle(0);
  tof.Clear();

  cout << "\n\nproton\n";
  xmin = 0.3, xmax = 12., ymin = -.5, ymax = 3.;
  tof.SetChi2Max(10000);
  std::unique_ptr<TH2D> hproton{(TH2D*) fIn->Get("reco_vs_sim_info/h2TofM2_p")};
  std::unique_ptr<TH2D> hproton_cut{(TH2D*) cutTH2(hproton, (TCutG*) fCuts->Get("proton"))};
  hproton_cut->Rebin2D(10, 5);
  TF1 fit_proton("fit_proton", "gaus", ymin, ymax);
  fit_proton.SetParNames("p6", "p7", "p8");
  fit_proton.SetParLimits(0, 0., 2.e6);
  fit_proton.SetParLimits(1, 0., 1.1);
  fit_proton.SetParLimits(2, 0., 1.2);
  TF1 proton_0("proton_0", "0", xmin, xmax);
  TF1 proton_1("proton_1", "pol11", .3, 12.);
  TF1 proton_2("proton_2", "pol8", .3, 12.);

  Pid::ParticleFit proton(PidTofParticles::kProton);
  proton.SetParametrization({proton_0, proton_1, proton_2});
  proton.SetFitFunction(fit_proton);
  proton.SetRange(xmin, xmax);
  proton.SetIsFitted();

  tof.AddParticle(proton, PidTofParticles::kPionPos);
  tof.SetHisto2D(std::move(hproton_cut));
  tof.SetRangeX(xmin, 8.);
  tof.SetRangeY(ymin, ymax);
  tof.SetOutputFileName("proton.root");
  tof.Fit();
  proton = tof.GetParticle(0);
  tof.Clear();

  cout << "\n\nhe3\n";
  std::unique_ptr<TH2D> hhe3{(TH2D*) fIn->Get("reco_vs_sim_info/h2TofM2_He3")};
  std::unique_ptr<TH2D> hhe3_cut{(TH2D*) cutTH2(hhe3, (TCutG*) fCuts->Get("he3"))};
  hhe3_cut->Rebin2D(10, 10);
  xmin = 1.5, xmax = 7., ymin = 1., ymax = 3.;
  tof.SetChi2Max(10);
  Pid::ParticleFit he3(PidTofParticles::kHe3);
  TF1 he3_fit("fit_he3", "gaus", ymin, ymax);
  he3_fit.SetParNames("p9", "p10", "p11");
  he3_fit.SetParLimits(0, 0., 1.e6);
  he3_fit.SetParLimits(1, 1.6, 2.2);
  he3_fit.SetParLimits(2, 0.08, 0.6);
  TF1 he3_0("he3_0", "0", xmin, xmax);
  TF1 he3_1("he3_1", "pol9", xmin, xmax);
  TF1 he3_2("he3_2", "pol5", xmin, xmax);

  he3.SetParametrization({he3_0, he3_1, he3_2});
  he3.SetFitFunction(he3_fit);
  he3.SetRange(xmin, xmax);
  he3.SetIsFitted();
  tof.AddParticle(he3, PidTofParticles::kHe3);
  tof.SetHisto2D(std::move(hhe3_cut));
  tof.SetRangeX(xmin, xmax);
  tof.SetRangeY(ymin, ymax);
  tof.SetOutputFileName("he3.root");
  tof.Fit();
  he3 = tof.GetParticle(0);
  tof.Clear();

  cout << "\n\ndeutron\n";
  xmin = 1.3, xmax = 11.5, ymin = 1., ymax = 6.;
  tof.SetChi2Max(4000);
  std::unique_ptr<TH2D> hdeutron{(TH2D*) fIn->Get("reco_vs_sim_info/h2TofM2_d")};
  std::unique_ptr<TH2D> hdeutron_cut{(TH2D*) cutTH2(hdeutron, (TCutG*) fCuts->Get("deutron"))};
  hdeutron_cut->Rebin2D(10, 25);
  TF1 fit_deutron("fit_deutron", "gaus", ymin, ymax);
  fit_deutron.SetParNames("p12", "p13", "p14");
  fit_deutron.SetParLimits(0, 0., 6.e4);
  fit_deutron.SetParLimits(1, 2., 4);
  fit_deutron.SetParLimits(2, 0., 1.);
  TF1 deutron_0("deutron_0", "0", xmin, xmax);
  TF1 deutron_1("deutron_1", "pol10", xmin, xmax);
  TF1 deutron_2("deutron_2", "pol3", xmin, xmax);

  Pid::ParticleFit deutron(PidTofParticles::kProton);
  deutron.SetParametrization({deutron_0, deutron_1, deutron_2});
  deutron.SetFitFunction(fit_deutron);
  deutron.SetRange(xmin, xmax);
  deutron.SetIsFitted();

  tof.AddParticle(deutron, PidTofParticles::kDeutron);
  tof.SetHisto2D(std::move(hdeutron_cut));
  tof.SetRangeX(xmin, xmax);
  tof.SetRangeY(ymin, ymax);
  tof.SetOutputFileName("deutron.root");
  tof.Fit();
  deutron = tof.GetParticle(0);
  tof.Clear();

  cout << "\n\nbgpos\n";
  xmin = 0.3, xmax = 20., ymin = -3., ymax = 3.;
  Pid::ParticleFit bgpos(PidTofParticles::kBgPos);
  TF1 bgpos_fit("fit_bgpos", "pol2", ymin, ymax);
  bgpos_fit.SetParNames("p15", "p16", "p17");
  TF1 bgpos_0("bgpos_0", "pol3", xmin, xmax);
  TF1 bgpos_1("bgpos_1", "pol3", xmin, xmax);
  TF1 bgpos_2("bgpos_2", "pol3", xmin, xmax);

  bgpos.SetParametrization({bgpos_0, bgpos_1, bgpos_2});
  bgpos.SetFitFunction(bgpos_fit);
  bgpos.SetRange(xmin, xmax);
  bgpos.SetIsFitted();

  cout << "\n\nallpos\n";
  xmin = 0.3, xmax = 12., ymin = -6., ymax = 6.;
  std::unique_ptr<TH2D> hpos{(TH2D*) fIn->Get("reco_info/h2TofM2")};
  hpos->Rebin2D(10, 1);
  tof.SetChi2Max(1e6);
  pionpos.SetIsFixed({false, true, true});
  kaonpos.SetIsFixed({false, true, true});
  proton.SetIsFixed({false, true, true});
  deutron.SetIsFixed({false, true, true});
  he3.SetIsFixed({false, true, true});
  tof.AddParticle(pionpos, PidTofParticles::kPionPos);
  tof.AddParticle(kaonpos, PidTofParticles::kKaonPos);
  tof.AddParticle(proton, PidTofParticles::kProton);
  tof.AddParticle(deutron, PidTofParticles::kDeutron);
  tof.AddParticle(he3, PidTofParticles::kHe3);
  tof.AddParticle(bgpos, PidTofParticles::kBgPos);

  tof.SetHisto2D(std::move(hpos));
  tof.SetRangeX(xmin, xmax);
  tof.SetRangeY(ymin, ymax);
  tof.SetOutputFileName("allpos.root");
  tof.Fit();

  getter.AddParticle(tof.GetParticleSpecie(PidTofParticles::kPionPos), PidTofParticles::kPionPos);
  getter.AddParticle(tof.GetParticleSpecie(PidTofParticles::kKaonPos), PidTofParticles::kKaonPos);
  getter.AddParticle(tof.GetParticleSpecie(PidTofParticles::kProton), PidTofParticles::kProton);
  getter.AddParticle(tof.GetParticleSpecie(PidTofParticles::kHe3), PidTofParticles::kHe3);
  getter.AddParticle(tof.GetParticleSpecie(PidTofParticles::kDeutron), PidTofParticles::kDeutron);
  getter.AddParticle(tof.GetParticleSpecie(PidTofParticles::kBgPos), PidTofParticles::kBgPos);

  tof.Clear();

  cout << "\n\npionneg\n";
  xmin = -10., xmax = -0.25, ymin = -1., ymax = 1.;
  tof.SetChi2Max(10000);
  std::unique_ptr<TH2D> hpionneg{(TH2D*) fIn->Get("reco_vs_sim_info/h2TofM2_piminus")};
  std::unique_ptr<TH2D> hpionneg_cut{(TH2D*) cutTH2(hpionneg, (TCutG*) fCuts->Get("piminus"))};
  hpionneg_cut->Rebin2D(10, 1);
  TF1 fit_pionneg("fit_pionneg", "gaus", ymin, ymax);
  fit_pionneg.SetParNames("p0", "p1", "p2");
  fit_pionneg.SetParLimits(0, 0., 5.e6);
  fit_pionneg.SetParLimits(1, -.075, 0.25);
  fit_pionneg.SetParLimits(2, 0., .33);
  TF1 pionneg_0("pionneg_0", "0", xmin, xmax);
  TF1 pionneg_1("pionneg_1", "pol2", xmin, xmax);
  TF1 pionneg_2("pionneg_2", "pol5", xmin, xmax);

  Pid::ParticleFit pionneg(PidTofParticles::kPionNeg);
  pionneg.SetParametrization({pionneg_0, pionneg_1, pionneg_2});
  pionneg.SetFitFunction(fit_pionneg);
  pionneg.SetRange(xmin, xmax);
  pionneg.SetIsFitted();

  tof.AddParticle(pionneg, PidTofParticles::kPionNeg);
  tof.SetHisto2D(std::move(hpionneg_cut));
  tof.SetRangeX(xmin, xmax);
  tof.SetRangeY(ymin, ymax);
  tof.SetOutputFileName("pionneg.root");
  tof.Fit();
  pionneg = tof.GetParticle(0);
  tof.Clear();

  cout << "\n\nkaonneg\n";
  xmin = -4., xmax = -0.25, ymin = -1., ymax = 1.5;
  tof.SetChi2Max(100);
  std::unique_ptr<TH2D> hkaonneg{(TH2D*) fIn->Get("reco_vs_sim_info/h2TofM2_kminus")};
  std::unique_ptr<TH2D> hkaonneg_cut{(TH2D*) cutTH2(hkaonneg, (TCutG*) fCuts->Get("kminus"))};
  hkaonneg_cut->Rebin2D(10, 2);
  TF1 fit_kaonneg("fit_kaonneg", "gaus", ymin, ymax);
  fit_kaonneg.SetParNames("p3", "p4", "p5");
  fit_kaonneg.SetParLimits(0, 0., 4.e4);
  fit_kaonneg.SetParLimits(1, 0.21, 0.25);
  fit_kaonneg.SetParLimits(2, 0.02, .16);
  TF1 kaonneg_0("kaonneg_0", "0", xmin, xmax);
  TF1 kaonneg_1("kaonneg_1", "pol2", xmin, xmax);
  TF1 kaonneg_2("kaonneg_2", "pol2", xmin, xmax);

  Pid::ParticleFit kaonneg(PidTofParticles::kKaonNeg);
  kaonneg.SetParametrization({kaonneg_0, kaonneg_1, kaonneg_2});
  kaonneg.SetFitFunction(fit_kaonneg);
  kaonneg.SetRange(xmin, xmax);
  kaonneg.SetIsFitted();

  tof.AddParticle(kaonneg, PidTofParticles::kKaonNeg);
  tof.SetHisto2D(std::move(hkaonneg_cut));
  tof.SetRangeX(xmin, xmax);
  tof.SetRangeY(ymin, ymax);
  tof.SetOutputFileName("kaonneg.root");
  tof.Fit();
  kaonneg = tof.GetParticle(0);
  tof.Clear();

  cout << "\n\nbgneg\n";
  xmin = -10., xmax = -0.25, ymin = -1., ymax = 2.;
  Pid::ParticleFit bgneg(PidTofParticles::kBgNeg);
  TF1 bgneg_0("bgneg_0", "pol3", xmin, xmax);//bgneg_0.SetParameters(100, 0, 0);
  TF1 bgneg_1("bgneg_1", "pol5", xmin, xmax);//bgneg_1.SetParameters(0, 0, 0);
  TF1 bgneg_2("bgneg_2", "pol5", xmin, xmax);//bgneg_2.SetParameters(0.0, 0.0, 0);

  TF1 bgneg_fit("fit_bgneg", "pol2", ymin, ymax);
  //TF1 bgneg_fit ("fit_bgneg", "[0]+[2]*(x-[1])*(x-[1])", ymin, ymax);
  //bgneg_fit.SetParLimits(0, 80., 150.);
  //bgneg_fit.SetParLimits(1, -.1, .1);
  //bgneg_fit.SetParLimits(2, -2000., 0.);
  bgneg_fit.SetParNames("p6", "p7", "p8");
  bgneg.SetParametrization({bgneg_0, bgneg_1, bgneg_2});
  bgneg.SetFitFunction(bgneg_fit);
  bgneg.SetRange(xmin, xmax);
  bgneg.SetIsFitted();

  cout << "\n\nallneg\n";
  xmin = -12, xmax = -0.25, ymin = -2., ymax = 6.;
  std::unique_ptr<TH2D> hneg{(TH2D*) fIn->Get("reco_info/h2TofM2")};
  hneg->Rebin2D(10, 1);
  tof.SetChi2Max(1e6);
  pionneg.SetIsFixed({false, true, true});
  kaonneg.SetIsFixed({false, true, true});
  tof.AddParticle(pionneg, PidTofParticles::kPionNeg);
  tof.AddParticle(kaonneg, PidTofParticles::kKaonNeg);
  tof.AddParticle(bgneg, PidTofParticles::kBgNeg);

  tof.SetHisto2D(std::move(hneg));
  tof.SetRangeX(xmin, xmax);
  tof.SetRangeY(ymin, ymax);
  tof.SetOutputFileName("allneg.root");
  tof.Fit();

  getter.AddParticle(tof.GetParticleSpecie(PidTofParticles::kPionNeg), PidTofParticles::kPionNeg);
  getter.AddParticle(tof.GetParticleSpecie(PidTofParticles::kKaonNeg), PidTofParticles::kKaonNeg);
  getter.AddParticle(tof.GetParticleSpecie(PidTofParticles::kBgNeg), PidTofParticles::kBgNeg);

  tof.Clear();

  std::unique_ptr<TFile> outfile{TFile::Open("pid_getter.root", "recreate")};
  getter.Write("pid_getter");
  outfile->Close();
}
