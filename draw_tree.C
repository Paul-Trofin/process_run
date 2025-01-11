#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TLatex.h>

void draw_tree() {
    // Open the ROOT file
    TFile *file = TFile::Open("run.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file 'run.root'" << std::endl;
        return;
    }

    // Get the trees for electron, muon, and lepton masses
    TTree *e_tree = (TTree*)file->Get("e_tree");
    TTree *mu_tree = (TTree*)file->Get("mu_tree");
    TTree *l_tree = (TTree*)file->Get("l_tree");

    if (!e_tree || !mu_tree || !l_tree) {
        std::cerr << "Error: Could not find one of the trees in 'run.root'" << std::endl;
        file->Close();
        return;
    }

    // Variables for branches
    double m_ee, m_mm, m_ll;
    
    gStyle->SetOptStat(0);

    // Set branch addresses for electron, muon, and lepton trees
    e_tree->SetBranchAddress("m_ee", &m_ee);
    mu_tree->SetBranchAddress("m_mm", &m_mm);
    l_tree->SetBranchAddress("m_ll", &m_ll);

    // Create histograms for each mass in the range 60-120 GeV
    TH1D *h_ee = new TH1D("h_ee", "Invariant Mass e- e+;m_{ee} (GeV);(Events/GeV)", 240, 60, 120);
    TH1D *h_mm = new TH1D("h_mm", "Invariant Mass mu- mu+;m_{#mu#mu} (GeV);(Events/GeV)", 240, 60, 120);
    TH1D *h_ll = new TH1D("h_ll", "Invariant Mass l- l+;m_{ll} (GeV);(Events/GeV)", 240, 60, 120);

    // Fill histograms from e_tree, mu_tree, and l_tree
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
    
    Long64_t l_entries = l_tree->GetEntries();
    for (Long64_t i = 0; i < l_entries; i++) {
        l_tree->GetEntry(i);
        h_ll->Fill(m_ll);
    }

    // Color and fill histograms
    h_ee->SetFillColor(kBlue-8); // Set fill color to blue for m_ee
    h_ee->SetFillStyle(1001); // Set solid fill style

    h_mm->SetFillColor(kRed-8); // Set fill color to red for m_mm
    h_mm->SetFillStyle(1001); // Set solid fill style

    h_ll->SetFillColor(kMagenta-8); // Set fill color to a mix of blue and red (purple) for m_ll
    h_ll->SetFillStyle(1001); // Set solid fill style

    // Create canvases
    TCanvas *c1 = new TCanvas("c1", "m_ee Plot", 1980, 1080);
    h_ee->Draw();

    // Fit the histogram with a Breit-Wigner function
    TF1 *bw_fit_ee = new TF1("bw_fit_ee", "[0] * [2] / ( (x - [1])^2 + [2]^2 ) + [3]/x^2", 60, 120); 
    bw_fit_ee->SetParameters(1000, 91.1876, 2.495); // Parameters: normalization, peak mass (Z mass), and width (Z width)
    h_ee->Fit(bw_fit_ee);

    // Add text for Z mass, sqrt(s), and integrated luminosity
    TLatex tex1;
    tex1.SetTextSize(0.04);
    tex1.SetTextAlign(12); // Left-aligned
    tex1.SetNDC();  // Ensure normalized coordinates
    tex1.DrawLatex(0.65, 0.88, "#sqrt{s} = 7 TeV");
    tex1.DrawLatex(0.7, 0.86, "Integrated Luminosity = 0.04 fb^{-1}");
    tex1.DrawLatex(0.7, 0.82, Form("Z Mass = %.2f GeV", bw_fit_ee->GetParameter(1)));

    c1->SaveAs("m_ee_plot.png");

    TCanvas *c2 = new TCanvas("c2", "m_mm Plot", 1980, 1080);
    h_mm->Draw();

    // Fit the histogram with a Breit-Wigner function
    TF1 *bw_fit_mm = new TF1("bw_fit_mm", "[0] * [2] / ( (x - [1])^2 + [2]^2 ) + [3]/x^2", 60, 120);
    bw_fit_mm->SetParameters(1000, 91.1876, 2.495); // Parameters: normalization, peak mass (Z mass), and width (Z width)
    h_mm->Fit(bw_fit_mm);

    // Add text for Z mass, sqrt(s), and integrated luminosity
    tex1.DrawLatex(0.65, 0.88, "#sqrt{s} = 7 TeV");
    tex1.DrawLatex(0.7, 0.86, "Integrated Luminosity = 0.04 fb^{-1}");
    tex1.DrawLatex(0.7, 0.82, Form("Z Mass = %.2f GeV", bw_fit_mm->GetParameter(1)));

    c2->SaveAs("m_mm_plot.png");

    TCanvas *c3 = new TCanvas("c3", "m_ll Plot", 1980, 1080);
    h_ll->Draw();

    // Fit the histogram with a Breit-Wigner function
    TF1 *bw_fit_ll = new TF1("bw_fit_ll", "[0] * [2] / ( (x - [1])^2 + [2]^2 ) + [3]/x^2", 60, 120);
    bw_fit_ll->SetParameters(1000, 91.1876, 2.495); // Parameters: normalization, peak mass (Z mass), and width (Z width)
    h_ll->Fit(bw_fit_ll);

    // Add text for Z mass, sqrt(s), and integrated luminosity
    tex1.DrawLatex(0.65, 0.88, "#sqrt{s} = 7 TeV");
    tex1.DrawLatex(0.7, 0.86, "Integrated Luminosity = 0.04 fb^{-1}");
    tex1.DrawLatex(0.7, 0.82, Form("Z Mass = %.2f GeV", bw_fit_ll->GetParameter(1)));

    c3->SaveAs("m_ll_plot.png");

    // Clean up
    file->Close();
}

