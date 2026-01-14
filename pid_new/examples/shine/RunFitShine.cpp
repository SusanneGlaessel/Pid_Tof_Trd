//
// Created by eugene on 26/11/2019.
//

#include <TFile.h>
#include <TH2D.h>

#include <RooDataHist.h>
#include <RooPlot.h>
#include <RooWorkspace.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TRatioPlot.h>
#include <pid_new/core/FitterHelper.h>

#include <TMath.h>
#include <TPaveText.h>
#include <pid_new/core/FitParameter.h>

#include <pid_new/core/ParticleFitModel.h>
#include <pid_new/model/AsymmetricGaussianPDF.h>

double poly(const std::vector<double>& pv, double x) {
  double result = 0.;
  double powX = 1.;

  for (auto p : pv) {
    result += p * powX;
    powX *= x;
  }

  return result;
}

FitParameter::ConstraintFct_t polyF(std::initializer_list<double> params) {
  std::vector<double> pv(params.begin(), params.end());
  return [pv](double x) { return poly(pv, x); };
}

FitParameter::ConstraintFct_t wrapToX(FitParameter::ConstraintFct_t f, int charge) {
  return [=](double x) { return f(TMath::Exp(charge * x) / 20); };
}

int main(int argc, char** argv) {

  /* Building RooFit model */

  /* observable */
  RooRealVar v_dEdx("x", "dEdx", 0, 2.5, "MIP");
  RooDataHist dh_dEdx("dh_dEdx", "", v_dEdx);

  /* common for all particle species asymmetry parameter */
  RooRealVar v_d("d", "#delta", -0.5, 0.5, "");

  /* negative pions model */
  RooRealVar v_pion_neg_bb("pion_neg_bb", "BB(#pi^{-})", 0, "MIP");
  RooRealVar v_pion_neg_sigma("pion_neg_sigma", "#sigma(#pi^{-})", 0, "MIP");
  RooRealVar v_pion_neg_n("pion_neg_n", "N(#pi^{-})", 0, "");
  AsymmetricGaussianPDF pdf_pion_neg("pion_neg", "", v_dEdx, v_pion_neg_bb, v_pion_neg_sigma, v_d);
  RooExtendPdf epdf_pion_neg("ext_pion_neg", "", pdf_pion_neg, v_pion_neg_n);

  /* negative kaons model */
  RooRealVar v_kaon_neg_bb("kaon_neg_bb", "BB(K^{-})", 0, "MIP");
  RooFormulaVar v_kaon_neg_sigma("kaon_neg_sigma", "#sigma(K^{-})", "@0*TMath::Power(@1/@2, 0.65)",
                                 RooArgSet(v_pion_neg_sigma, v_kaon_neg_bb, v_pion_neg_bb));
  RooRealVar v_kaon_neg_n("kaon_neg_n", "N(K^{-})", 0, "");
  AsymmetricGaussianPDF pdf_kaon_neg("kaon_neg", "", v_dEdx, v_kaon_neg_bb, v_kaon_neg_sigma, v_d);
  RooExtendPdf epdf_kaon_pdf("ext_kaon_neg", "", pdf_kaon_neg, v_kaon_neg_n);

  /* electrons model */
  RooRealVar v_electron_bb("electron_bb", "BB(e^{-})", 0, "MIP");
  RooFormulaVar v_electron_sigma("electron_sigma", "#sigma(e^{-})", "@0*TMath::Power(@1/@2, 0.65)",
                                 RooArgSet(v_pion_neg_sigma, v_electron_bb, v_pion_neg_bb));
  RooRealVar v_electron_n("electron_n", "N(e^{-})", 0, "");
  AsymmetricGaussianPDF pdf_electron("electron", "", v_dEdx, v_electron_bb, v_electron_sigma, v_d);
  RooExtendPdf epdf_electron("ext_electron", "", pdf_electron, v_electron_n);

  /* positive pions */
  RooRealVar v_pion_pos_bb("pion_pos_bb", "BB(#pi^{+})", 0, "MIP");
  RooRealVar v_pion_pos_sigma("pion_pos_sigma", "#sigma(#pi^{+})", 0, "MIP");
  RooRealVar v_pion_pos_n("pion_pos_n", "N(#pi^{+})", 0, "");
  AsymmetricGaussianPDF pdf_pion_pos("pion_pos", "", v_dEdx, v_pion_pos_bb, v_pion_pos_sigma, v_d);
  RooExtendPdf epdf_pion_pos("ext_pion_pos", "", pdf_pion_pos, v_pion_pos_n);

  /* protons model */
  RooRealVar v_proton_bb("proton_bb", "BB(p)", 0, "MIP");
  RooFormulaVar v_proton_sigma("proton_sigma", "#sigma(p)", "@0*TMath::Power(@1/@2, 0.65)",
                               RooArgSet(v_pion_pos_sigma, v_proton_bb, v_pion_pos_bb));
  RooRealVar v_proton_n("proton_n", "N(p)", 0, "");
  AsymmetricGaussianPDF pdf_proton("proton", "", v_dEdx, v_proton_bb, v_proton_sigma, v_d);
  RooExtendPdf epdf_proton("ext_proton", "", pdf_proton, v_proton_n);

  /* positrons model */
  RooRealVar v_positron_bb("positron_bb", "BB(e^{-})", 0, "MIP");
  RooFormulaVar v_positron_sigma("positron_sigma", "#sigma(e^{-})", "@0*TMath::Power(@1/@2, 0.65)",
                                 RooArgSet(v_pion_pos_sigma, v_positron_bb, v_pion_pos_bb));
  RooRealVar v_positron_n("positron_n", "N(e^{-})", 0, "");
  AsymmetricGaussianPDF pdf_positron("positron", "", v_dEdx, v_positron_bb, v_positron_sigma, v_d);
  RooExtendPdf epdf_positron("ext_positron", "", pdf_positron, v_positron_n);

  /* kaons model */
  RooRealVar v_kaon_pos_bb("kaon_pos_bb", "BB(K^{+})", 0, "MIP");
  RooFormulaVar v_kaon_pos_sigma("kaon_pos_sigma", "#sigma(K^{+})", "@0*TMath::Power(@1/@2, 0.65)",
                                 RooArgSet(v_pion_pos_sigma, v_kaon_pos_bb, v_pion_pos_bb));
  RooRealVar v_kaon_pos_n("kaon_pos_n", "N(K^{+})", 0, "");
  AsymmetricGaussianPDF pdf_kaon_pos("kaon_pos", "", v_dEdx, v_kaon_pos_bb, v_kaon_pos_sigma, v_d);
  RooExtendPdf epdf_kaon_pos_pdf("ext_kaon_pos", "", pdf_kaon_pos, v_kaon_pos_n);

  std::list<ParticleFitModel> fit_models;

  /* negative */
  fit_models.emplace_back("pion_neg", epdf_pion_neg, dh_dEdx, "pion_neg_");
  fit_models.back().setRange(-5.5, -0.7);
  fit_models.back().parameter("bb").fix(wrapToX(BetheBlochHelper::makeBBForPdg(-211), -1));
  fit_models.back().parameter("sigma").fix(polyF({1.17732e-01, 8.11273e-03}), -5.5, -3.0);
  fit_models.back().parameter("sigma").fix(polyF({1.28975e-01, 1.14362e-02}), -3.0, -2.6);
  fit_models.back().parameter("sigma").fix(polyF({4.81931e-01, 7.07309e-01, 5.30289e-01, 1.81315e-01, 2.31737e-02}), -2.6, -0.4);
  fit_models.back().parameter("n").fix(polyF({5.98777e+05, 2.17637e+05, 2.35418e+04, 1.51676e+03, 1.53068e+02}), -5.5, -4.5);
  fit_models.back().parameter("n").fix(polyF({-5.56917e+04, -1.64847e+05, -1.78915e+05, -8.46502e+04, -1.34296e+04}), -2.2, -0.4);

  fit_models.emplace_back("kaon_neg", epdf_kaon_pdf, dh_dEdx, "kaon_neg_");
  fit_models.back().setRange(-5.3, -2.5);
  fit_models.back().parameter("bb").fix(wrapToX(BetheBlochHelper::makeBBForPdg(-321), -1));
  fit_models.back().parameter("n").range(0, 2000);

  fit_models.emplace_back("electron", epdf_electron, dh_dEdx, "electron_");
  fit_models.back().setRange(-5., -0.4);
  fit_models.back().parameter("bb").fix(wrapToX(BetheBlochHelper::makeBBForPdg(-11), -1));

  /* positive */
  fit_models.emplace_back("pion_pos", epdf_pion_pos, dh_dEdx, "pion_pos_");
  fit_models.back().setRange(0.4, 5.5);
  fit_models.back().parameter("bb").fix(wrapToX(BetheBlochHelper::makeBBForPdg(211), 1));
  fit_models.back().parameter("sigma").fix(polyF({4.81931e-01, -7.07309e-01, 5.30289e-01, -1.81315e-01, 2.31737e-02}), 0.4, 2.6);
  fit_models.back().parameter("sigma").fix(polyF({1.28975e-01, -1.14362e-02}), 2.6, 3.0);
  fit_models.back().parameter("sigma").fix(polyF({1.17732e-01, -8.11273e-03}), 3.0, 5.5);

  fit_models.emplace_back("proton", epdf_proton, dh_dEdx, "proton_");
  fit_models.back().setRange(2.2, 5.5);
  fit_models.back().parameter("bb").fix(wrapToX(BetheBlochHelper::makeBBForPdg(2212), 1));

  fit_models.emplace_back("positron", epdf_positron, dh_dEdx, "positron_");
  fit_models.back().setRange(0.4, 5.0);
  fit_models.back().parameter("bb").fix(wrapToX(BetheBlochHelper::makeBBForPdg(11), 1));

  fit_models.emplace_back("kaon_pos", epdf_kaon_pos_pdf, dh_dEdx, "kaon_pos_");
  fit_models.back().setRange(2.5, 5.5);
  fit_models.back().parameter("bb").fix(wrapToX(BetheBlochHelper::makeBBForPdg(321), 1));
  fit_models.back().parameter("n").range(0, 2000);

  FitParameterRegistry::instance().par_find("d")->fix(0.05);

  /* IO */
  auto inputFile = TFile::Open(argv[1]);
  TFile outputFile("output.root", "RECREATE");

  TH2D* inputHistogram = nullptr;
  inputFile->GetObject("reco_info/hTrackdEdx_allTPC", inputHistogram);
  inputHistogram->SetDirectory(nullptr);
  inputHistogram->Print();

  inputHistogram->RebinX(1);
  inputHistogram->RebinY(1);

  auto c = new TCanvas;
  c->Print("output.pdf(", "pdf");

  inputHistogram->Draw("colz");
  gPad->SetLogz();
  inputHistogram->GetYaxis()->SetRangeUser(0, 4.);

  {
    TF1 bbPionNeg(
        "bbPionNeg", [](const Double_t* x, const Double_t* p) {
          return wrapToX(BetheBlochHelper::makeBBForPdg(-211), -1).operator()(x[0]);
        },
        -6., -0.5, 0);
    bbPionNeg.SetLineColor(kRed);
    bbPionNeg.Draw("same");

    TF1 bbPionPos(
        "bbPionPos", [](const Double_t* x, const Double_t* p) {
          return wrapToX(BetheBlochHelper::makeBBForPdg(211), 1).operator()(x[0]);
        },
        0.5, 6, 0);
    bbPionPos.SetLineColor(kMagenta + 1);
    bbPionPos.Draw("same");

    TF1 bbProton(
        "bbProton", [](const Double_t* x, const Double_t* p) {
          return wrapToX(BetheBlochHelper::makeBBForPdg(2212), 1).operator()(x[0]);
        },
        2.5, 6, 0);
    bbProton.SetLineColor(kBlack);
    bbProton.Draw("same");

    TF1 bbPositron(
        "bbPositron", [](const Double_t* x, const Double_t* p) {
          return wrapToX(BetheBlochHelper::makeBBForPdg(11), 1).operator()(x[0]);
        },
        0., 5., 0);
    bbPositron.SetLineColor(kBlue + 1);
    bbPositron.Draw("same");

    TF1 bbElectron(
        "bbElectron", [](const Double_t* x, const Double_t* p) {
          return wrapToX(BetheBlochHelper::makeBBForPdg(-11), -1).operator()(x[0]);
        },
        -5., 0., 0);
    bbElectron.SetLineColor(kBlue + 1);
    bbElectron.Draw("same");

    TF1 bbDeuteron(
        "bbDeuteron", [](const Double_t* x, const Double_t* p) {
          return wrapToX(BetheBlochHelper::makeBBForPdg(1000010020), 1).operator()(x[0]);
        },
        3., 6., 0);
    bbDeuteron.SetLineColor(kOrange + 1);
    bbDeuteron.Draw("same");

    TF1 bbKaonPos(
        "bbKaonPos", [](const Double_t* x, const Double_t* p) {
          return wrapToX(BetheBlochHelper::makeBBForPdg(321), 1).operator()(x[0]);
        },
        2., 5., 0);
    bbKaonPos.SetLineColor(kCyan + 1);
    bbKaonPos.Draw("same");

    TF1 bbKaonNeg(
        "bbKaonNeg", [](const Double_t* x, const Double_t* p) {
          return wrapToX(BetheBlochHelper::makeBBForPdg(-321), -1).operator()(x[0]);
        },
        -5., -2., 0);
    bbKaonNeg.SetLineColor(kCyan + 1);
    bbKaonNeg.Draw("same");

    c->Print("output.pdf(", "pdf");
    c->Clear();
  }

  for (int i = 1; i < inputHistogram->GetXaxis()->GetNbins(); ++i) {
    double x = inputHistogram->GetXaxis()->GetBinCenter(i);
    double xlo = inputHistogram->GetXaxis()->GetBinLowEdge(i);
    double xhi = inputHistogram->GetXaxis()->GetBinUpEdge(i);

    std::cout << "Bin: " << i << "; x = " << x << std::endl;

    /* building composite model */
    RooArgList composite_model_components;
    std::for_each(fit_models.begin(), fit_models.end(), [x, &composite_model_components](ParticleFitModel& m) {
      if (m.isDefinedAt(x)) {
        composite_model_components.add(m.pdf());
      }
    });

    if (composite_model_components.empty()) {
      std::cout << "Nothing to fit in this momentum bin" << std::endl;
      continue;
    }

    RooAddPdf pdf_composite(Form("composite_%d", i), "", composite_model_components);
    ParticleFitModel m_composite("composite", pdf_composite, dh_dEdx);

    /* preparing input */
    auto py = inputHistogram->ProjectionY("tmp", i, i);
    py->SetDirectory(nullptr);

    dh_dEdx.reset();
    RooDataHist slice_data("slice_dh", "", v_dEdx, py);
    dh_dEdx.add(slice_data);

    auto qa_dir = outputFile.mkdir(Form("bin_%d", i));

    auto slice_integral = dh_dEdx.sumEntries();

    /* set range for the component integrals */
    auto integral_params = m_composite.parameters_matching(".*_n$");
    std::for_each(integral_params.begin(), integral_params.end(), [=](FitParameter* p) {
      auto& constr_lo = p->findConstraint(xlo);
      auto& constr_hi = p->findConstraint(xhi);

      if (
          constr_lo.type_ == FitParameter::EConstraintType::kNone && constr_hi.type_ == FitParameter::EConstraintType::kNone) {
        p->range(0., slice_integral, xlo, xhi);
      }
    });
    /* constraints */
    auto model_params = m_composite.parameters();
    std::for_each(model_params.begin(), model_params.end(), [=](FitParameter* p) { p->applyConstraint(x); });

    m_composite.pdf().fitTo(dh_dEdx, RooFit::Extended());
    std::for_each(model_params.begin(), model_params.end(), [=](FitParameter* p) { p->pickFitResultAt(x); });

    auto frame_data_fit = v_dEdx.frame();
    dh_dEdx.plotOn(frame_data_fit, RooFit::Name("data"));
    m_composite.pdf().plotOn(frame_data_fit, RooFit::Name("fit"));

    for (auto& m : fit_models) {
      if (m.isDefinedAt(x)) {
        m.pdf().plotOn(frame_data_fit,
                       RooFit::Normalization(1., RooAbsReal::RelativeExpected),
                       RooFit::LineColor(m.getColor()));
      }
    }

    qa_dir->WriteObject(frame_data_fit, "data_fit");

    TPaveText pave_text(0.1, 0.9, 0.4, 0.98);
    pave_text.AddText(Form("x = %f", x));

    c->cd();
    c->Divide(2, 2);

    c->cd(1);
    frame_data_fit->Draw();
    pave_text.DrawClone();

    c->cd(2);
    frame_data_fit->Draw();
    gPad->SetLogy(1);

    c->cd(3);
    auto hist_data = dh_dEdx.createHistogram("tmp_data", v_dEdx, RooFit::IntrinsicBinning());
    auto hist_fit = m_composite.pdf().createHistogram("tmp_fit", v_dEdx, RooFit::IntrinsicBinning(),
                                                      RooFit::Extended(kTRUE), RooFit::Scaling(kTRUE));
    hist_fit->Scale(v_dEdx.getBinWidth(1));

    TRatioPlot rp(hist_fit, hist_data, "divsym");
    rp.Draw("nohide,grid");
    rp.GetLowerRefGraph()->SetMinimum(0.5);
    rp.GetLowerRefGraph()->SetMaximum(1.5);

    c->cd(4);

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    TLegend species_legend(0.3, 0.85, 0.7, 0.95);
    species_legend.SetNColumns(2);

    bool same_flag = false;
    bool legend_flag = false;
    /* calculate purity */
    for (auto& m : fit_models) {
      if (m.isDefinedAt(x)) {
        auto model_color = m.getColor();

        auto particle_pdf_hist = (TH1F*) m.pdf().createHistogram("tmp_component", v_dEdx,
                                                                 RooFit::IntrinsicBinning());
        auto composite_pdf_hist = (TH1F*) m_composite.pdf().createHistogram("tmp_composite", v_dEdx,
                                                                            RooFit::IntrinsicBinning(),
                                                                            RooFit::Extended());
        auto data_hist = (TH1F*) dh_dEdx.createHistogram("tmp_data", v_dEdx, RooFit::IntrinsicBinning());
        qa_dir->WriteObject(particle_pdf_hist, Form("pdf_%s", m.getName().c_str()));

        auto purity_hist = (*particle_pdf_hist) / (*composite_pdf_hist);
        purity_hist.SetName(Form("purity_%s", m.getName().c_str()));
        purity_hist.SetLineColor(model_color);
        purity_hist.SetMarkerColor(model_color);
        purity_hist.SetLineWidth(2);
        qa_dir->WriteObject(&purity_hist, purity_hist.GetName());

        particle_pdf_hist->Scale(v_dEdx.getBinWidth(1));
        auto purity_hist_data = (*particle_pdf_hist) / (*data_hist);
        purity_hist_data.SetName(Form("purity_to_data_%s", m.getName().c_str()));
        purity_hist_data.SetLineColor(model_color);
        purity_hist_data.SetMarkerColor(model_color);
        qa_dir->WriteObject(&purity_hist_data, purity_hist_data.GetName());

        purity_hist_data.GetYaxis()->SetRangeUser(0., 2.);
        purity_hist_data.DrawClone(same_flag ? "same,LE1" : "LE1");
        same_flag = true;
        purity_hist.DrawClone("same");

        species_legend.AddEntry(Form("purity_%s", m.getName().c_str()), m.getName().c_str(), "l");

        delete particle_pdf_hist;
        delete composite_pdf_hist;
        delete data_hist;
      }
    }

    species_legend.Draw();

    c->Print("output.pdf", "pdf");
    c->Clear();

    delete hist_fit;
    delete hist_data;
    delete py;
  }

  std::for_each(FitParameterRegistry::instance().par_begin(), FitParameterRegistry::instance().par_end(),
                [&outputFile, c](FitParameter& p) {
                  p.dumpResult(outputFile, c);
                  c->Print("output.pdf", "pdf");
                  c->Clear();
                });

  c->Print("output.pdf)", "pdf");
  return 0;
}
