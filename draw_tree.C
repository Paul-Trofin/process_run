#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TLatex.h>

void draw_tree() {
    // OPEN ROOT FILE
    TFile *file = TFile::Open("run.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file 'run.root'" << std::endl;
        return;
    }

    // GET TREES
    TTree *e_tree = (TTree*)file->Get("e_tree");
    TTree *mu_tree = (TTree*)file->Get("mu_tree");
    TTree *tau_tree = (TTree*)file->Get("tau_tree");
    TTree *l_tree = (TTree*)file->Get("l_tree");
    TTree *l_k_tree = (TTree*)file->Get("l_k_tree");

    if (!e_tree || !mu_tree || !tau_tree || !l_tree || !l_k_tree) {
        std::cerr << "Error: Could not find one of the trees in 'run.root'" << std::endl;
        file->Close();
        return;
    }

    // Variables for branches
    double m_ee, m_mm, m_ll, m_tt;
    double l_px, l_py, l_pz, l_e, l_pT, l_eta, l_phi;
    
    gStyle->SetOptStat(0);

    // GET BRANCHES FROM TREES
    e_tree->SetBranchAddress("m_ee", &m_ee);
    mu_tree->SetBranchAddress("m_mm", &m_mm);
    tau_tree->SetBranchAddress("m_tt", &m_tt);
    l_tree->SetBranchAddress("m_ll", &m_ll);
    l_k_tree->SetBranchAddress("l_px", &l_px);
    l_k_tree->SetBranchAddress("l_py", &l_py);
    l_k_tree->SetBranchAddress("l_pz", &l_pz);
    l_k_tree->SetBranchAddress("l_e", &l_e);
    l_k_tree->SetBranchAddress("l_pT", &l_pT);
    l_k_tree->SetBranchAddress("l_eta", &l_eta);
    l_k_tree->SetBranchAddress("l_phi", &l_phi);

    // INIT HISTOGRAMS
    TH1D *h_ee = new TH1D("h_ee", "Invariant Mass e- e+;m_{ee} (GeV);(Events/GeV)", 240, 10, 120);
    TH1D *h_mm = new TH1D("h_mm", "Invariant Mass mu- mu+;m_{#mu#mu} (GeV);(Events/GeV)", 240, 10, 120);
    TH1D *h_tt = new TH1D("h_tt", "Invariant Mass ta- ta+;m_{#ta#ta} (GeV);(Events/GeV)", 240, 10, 120);
    TH1D *h_ll = new TH1D("h_ll", "Invariant Mass l- l+;m_{ll} (GeV);(Events/GeV)", 240, 10, 120);
    TH1D *h_l_px = new TH1D("h_px", "Lepton p_{x};p_{x} (GeV);Entries", 100, -100, 100);
    TH1D *h_l_py = new TH1D("h_py", "Lepton p_{y};p_{y} (GeV);Entries", 100, -100, 100);
    TH1D *h_l_pz = new TH1D("h_pz", "Lepton p_{z};p_{z} (GeV);Entries", 100, -500, 500);
    TH1D *h_l_e = new TH1D("h_e", "Lepton Energy;E (GeV);Entries", 100, 0, 300);
    TH1D *h_l_pT = new TH1D("h_pT", "Lepton p_{T};p_{T} (GeV);Entries", 100, 0, 100);
    TH1D *h_l_eta = new TH1D("h_eta", "Lepton #eta;#eta;Entries", 100, -10, 10);
    TH1D *h_l_phi = new TH1D("h_phi", "Lepton #phi;#phi;Entries", 100, -4, 4);

    // FILL HISTOGRAMS
    Long64_t e_entries = e_tree->GetEntries();
    for (Long64_t i = 0; i < e_entries; i++) {
        e_tree->GetEntry(i);
        h_ee->Fill(m_ee);
    }

    Long64_t mu_entries = mu_tree->GetEntries();
    for (Long64_t i = 0; i < mu_entries; i++) {
        mu_tree->GetEntry(i);
        h_mm->Fill(m_mm);
    }
    
    Long64_t tau_entries = tau_tree->GetEntries();
    for (Long64_t i = 0; i < tau_entries; i++) {
        tau_tree->GetEntry(i);
        h_tt->Fill(m_tt);
    }
    
    Long64_t l_entries = l_tree->GetEntries();
    for (Long64_t i = 0; i < l_entries; i++) {
        l_tree->GetEntry(i);
        h_ll->Fill(m_ll);
    }
    
    Long64_t l_k_entries = l_k_tree->GetEntries();
    for (Long64_t i = 0; i < l_k_entries; i++) {
        l_k_tree->GetEntry(i);
        h_l_px->Fill(l_px);
        h_l_py->Fill(l_py);
        h_l_pz->Fill(l_pz);
        h_l_e->Fill(l_e);
        h_l_pT->Fill(l_pT);
        h_l_eta->Fill(l_eta);
        h_l_phi->Fill(l_phi);
    }

    // CUSTOMIZE HISTOGRAMS
    h_ee->SetFillColor(kBlue-8);
    h_mm->SetFillColor(kRed-8);   
    h_tt->SetFillColor(kGreen-8);
    h_ll->SetFillColor(kMagenta-8);
    

    // CANVAS 1
    TCanvas *c1 = new TCanvas("c1", "m_ee Plot", 1980, 1080);
    h_ee->Draw();
    // Fit with a Breit-Wigner
    TF1 *bw_fit_ee = new TF1("bw_fit_ee", "[0] * [2] / ( (x - [1])^2 + [2]^2 ) + [3]/x^2", 60, 120); 
    bw_fit_ee->SetParameters(1000, 91.1876, 2.495);
    bw_fit_ee->SetNpx(1000);
    bw_fit_ee->SetLineWidth(4);
    bw_fit_ee->SetLineColor(kBlack);
    h_ee->Fit(bw_fit_ee);
	// SAVE CANVAS 1
    c1->SaveAs("m_ee_plot.png");
	
	// CANVAS 2
    TCanvas *c2 = new TCanvas("c2", "m_mm Plot", 1980, 1080);
    h_mm->Draw();
    // Fit with a Breit-Wigner
    TF1 *bw_fit_mm = new TF1("bw_fit_mm", "[0] * [2] / ( (x - [1])^2 + [2]^2 ) + [3]/x^2", 60, 120);
    bw_fit_mm->SetParameters(1000, 91.1876, 2.495);
    bw_fit_mm->SetNpx(1000);
    bw_fit_mm->SetLineWidth(4);
    bw_fit_mm->SetLineColor(kBlack);
    h_mm->Fit(bw_fit_mm);
    // SAVE CANVAS 2
    c2->SaveAs("m_mm_plot.png");
	
	// CANVAS 3
    TCanvas *c3 = new TCanvas("c3", "m_ll Plot", 1980, 1080);
    h_ll->Draw();
    // Fit with a Breit-Wigner
    TF1 *bw_fit_ll = new TF1("bw_fit_ll", "[0] * [2] / ( (x - [1])^2 + [2]^2 ) + [3]/x^2", 60, 120);
    bw_fit_ll->SetParameters(1000, 91.1876, 2.495);
    bw_fit_ll->SetNpx(1000);
    bw_fit_ll->SetLineWidth(4);
    bw_fit_ll->SetLineColor(kBlack);
    h_ll->Fit(bw_fit_ll);
	// SAVE CANVAS 4
    c3->SaveAs("m_ll_plot.png");
    
    // CANVAS 4
    TCanvas *c4 = new TCanvas("c4", "m_tt Plot", 1980, 1080);
    h_tt->Draw();
    // Fit with a Breit-Wigner
    TF1 *bw_fit_tt = new TF1("bw_fit_tt", "[0] * [2] / ( (x - [1])^2 + [2]^2 ) + [3]/x^2", 60, 120);
    bw_fit_tt->SetParameters(1000, 91.1876, 2.495);
    bw_fit_tt->SetNpx(1000);
    bw_fit_tt->SetLineWidth(4);
    bw_fit_tt->SetLineColor(kBlack);
    h_tt->Fit(bw_fit_ll);
	// SAVE CANVAS 4
    c4->SaveAs("m_tt_plot.png");
    
    // CANVAS 5
	TCanvas *c5 = new TCanvas("c5", "l_px Plot", 1980, 1080);
	h_l_px->Draw();
	c5->SaveAs("l_px_plot.png");

	// CANVAS 6
	TCanvas *c6 = new TCanvas("c6", "l_py Plot", 1980, 1080);
	h_l_py->Draw();
	c6->SaveAs("l_py_plot.png");

	// CANVAS 7
	TCanvas *c7 = new TCanvas("c7", "l_pz Plot", 1980, 1080);
	h_l_pz->Draw();
	c7->SaveAs("l_pz_plot.png");

	// CANVAS 8
	TCanvas *c8 = new TCanvas("c8", "l_e Plot", 1980, 1080);
	h_l_e->Draw();
	c8->SaveAs("l_e_plot.png");

	// CANVAS 9
	TCanvas *c9 = new TCanvas("c9", "l_pT Plot", 1980, 1080);
	h_l_pT->Draw();
	c9->SaveAs("l_pT_plot.png");

	// CANVAS 10
	TCanvas *c10 = new TCanvas("c10", "l_eta Plot", 1980, 1080);
	h_l_eta->Draw();
	c10->SaveAs("l_eta_plot.png");

	// CANVAS 11
	TCanvas *c11 = new TCanvas("c11", "l_phi Plot", 1980, 1080);
	h_l_phi->Draw();
	c11->SaveAs("l_phi_plot.png");

    // Clean up
    file->Close();
}
