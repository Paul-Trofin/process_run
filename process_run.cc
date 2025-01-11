#include <iostream>
#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TTree.h"
#include <cmath>
using namespace Pythia8;

int main() {
    Pythia pythia;

    //----- Suppress Pythia Output -------//
    pythia.readString("Next:numberShowEvent = 0");
    pythia.readString("Init:showAllSettings = off");
    pythia.readString("Init:showChangedSettings = off");
    pythia.readString("Print:quiet = on");
    
    //----- MY OUTPUT -----//
    std::cout << "//----- MY OUTPUT -----//";

    //----- READ LHE -------//
    pythia.readString("Beams:frameType = 4");            
    pythia.readString("Beams:LHEF = unweighted_events.lhe.gz");

    //----- INITIALIZE -----//
    if (!pythia.init()) {
        std::cout << "Pythia initialization failed." << std::endl;
        return 1;
    }

    // ROOT FILE
    TFile *output = new TFile("run.root", "RECREATE");

    // ELECTRON TREE
    TTree *e_tree = new TTree("e_tree", "e_tree");
    double m_ee;
    e_tree->Branch("m_ee", &m_ee, "m_ee/D");
    
    // MUON TREE
    TTree *mu_tree = new TTree("mu_tree", "mu_tree");
    double m_mm;
    mu_tree->Branch("m_mm", &m_mm, "m_mm/D");
    
    // LEPTON TREE
    TTree *l_tree = new TTree("l_tree", "l_tree");
    double m_ll;
    l_tree->Branch("m_ll", &m_ll, "m_ll/D");
    

    //----- FILL TREE -----//
    while (pythia.next()) {

        std::vector<Particle> electrons;
        std::vector<Particle> positrons;
        std::vector<Particle> muons;
        std::vector<Particle> antimuons;

        // EACH PARTICLE IN CURRENT EVENT
        for (int i = 0; i < pythia.event.size(); ++i) {
            if (pythia.event[i].isFinal()) {
                if (pythia.event[i].id() == 11) { // Electron
                    electrons.push_back(pythia.event[i]);
                }
                if (pythia.event[i].id() == -11) { // Positron
                    positrons.push_back(pythia.event[i]);
                }
                if (pythia.event[i].id() == 13) { // Muon
                    muons.push_back(pythia.event[i]);
                }
                if (pythia.event[i].id() == -13) { // Antimuon
                    antimuons.push_back(pythia.event[i]);
                }
            }
        }

        // e- e+ invariant mass m_ee
        for (auto& e : electrons) {
            for (auto& e_bar : positrons) {
                double px_tot = e.px() + e_bar.px();
                double py_tot = e.py() + e_bar.py();
                double pz_tot = e.pz() + e_bar.pz();
                double e_tot = e.e() + e_bar.e();
                m_ee = sqrt(e_tot * e_tot - px_tot * px_tot - py_tot * py_tot - pz_tot * pz_tot);
                m_ll = m_ee;
                e_tree->Fill();
                l_tree->Fill();
                
            }
        }

        // mu- mu+ invariant mass m_mm
        for (auto& mu : muons) {
            for (auto& mu_bar : antimuons) {
                double px_tot = mu.px() + mu_bar.px();
                double py_tot = mu.py() + mu_bar.py();
                double pz_tot = mu.pz() + mu_bar.pz();
                double e_tot = mu.e() + mu_bar.e();
                m_mm = sqrt(e_tot * e_tot - px_tot * px_tot - py_tot * py_tot - pz_tot * pz_tot);
                m_ll = m_mm;
                mu_tree->Fill();
                l_tree->Fill();
                
            }
        }
    }   
	
    output->Write();
    output->Close();
    std::cout << "The file run.root has been created.\n";
    return 0;
}

