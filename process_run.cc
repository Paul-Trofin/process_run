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
    std::cout << "//----- MY OUTPUT -----//\n";

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

    // ELECTRON TREES
    // Electron Invariant Mass
    TTree *e_tree = new TTree("e_tree", "e_tree"); // e inv mass
    double m_ee;
    e_tree->Branch("m_ee", &m_ee, "m_ee/D");
    // Electron kinematics
    TTree *e_k_tree = new TTree("e_k_tree", "e_k_tree");
    double e_px, e_py, e_pz, e_e, e_eta, e_pT, e_phi;
    e_k_tree->Branch("e_px", &e_px, "e_px/D");
    e_k_tree->Branch("e_py", &e_py, "e_py/D");
    e_k_tree->Branch("e_pz", &e_pz, "e_pz/D");
    e_k_tree->Branch("e_e", &e_e, "e_e/D");
    e_k_tree->Branch("e_pT", &e_pT, "e_pT/D");
    e_k_tree->Branch("e_eta", &e_eta, "e_eta/D");
    e_k_tree->Branch("e_phi", &e_phi, "e_phi/D");
    
    // MUON TREE
    TTree *mu_tree = new TTree("mu_tree", "mu_tree");
    double m_mm;
    mu_tree->Branch("m_mm", &m_mm, "m_mm/D");
    // Muon kinematics
    TTree *mu_k_tree = new TTree("mu_k_tree", "mu_k_tree");
	double mu_px, mu_py, mu_pz, mu_e, mu_eta, mu_pT, mu_phi;
	mu_k_tree->Branch("mu_px", &mu_px, "mu_px/D");
	mu_k_tree->Branch("mu_py", &mu_py, "mu_py/D");
	mu_k_tree->Branch("mu_pz", &mu_pz, "mu_pz/D");
	mu_k_tree->Branch("mu_e", &mu_e, "mu_e/D");
	mu_k_tree->Branch("mu_pT", &mu_pT, "mu_pT/D");
	mu_k_tree->Branch("mu_eta", &mu_eta, "mu_eta/D");
	mu_k_tree->Branch("mu_phi", &mu_phi, "mu_phi/D");

    
    // TAU TREE
    TTree *tau_tree = new TTree("tau_tree", "tau_tree");
    double m_tt;
    tau_tree->Branch("m_tt", &m_tt, "m_tt/D");
    // Tau kinematics
    TTree *tau_k_tree = new TTree("tau_k_tree", "tau_k_tree");
	double tau_px, tau_py, tau_pz, tau_e, tau_eta, tau_pT, tau_phi;
	tau_k_tree->Branch("tau_px", &tau_px, "tau_px/D");
	tau_k_tree->Branch("tau_py", &tau_py, "tau_py/D");
	tau_k_tree->Branch("tau_pz", &tau_pz, "tau_pz/D");
	tau_k_tree->Branch("tau_e", &tau_e, "tau_e/D");
	tau_k_tree->Branch("tau_pT", &tau_pT, "tau_pT/D");
	tau_k_tree->Branch("tau_eta", &tau_eta, "tau_eta/D");
	tau_k_tree->Branch("tau_phi", &tau_phi, "tau_phi/D");

    
    // LEPTON TREE
    TTree *l_tree = new TTree("l_tree", "l_tree");
    double m_ll;
    l_tree->Branch("m_ll", &m_ll, "m_ll/D");
    // Lepton kinematics
    TTree *l_k_tree = new TTree("l_k_tree", "l_k_tree");
	double l_px, l_py, l_pz, l_e, l_eta, l_pT, l_phi;
	l_k_tree->Branch("l_px", &l_px, "l_px/D");
	l_k_tree->Branch("l_py", &l_py, "l_py/D");
	l_k_tree->Branch("l_pz", &l_pz, "l_pz/D");
	l_k_tree->Branch("l_e", &l_e, "l_e/D");
	l_k_tree->Branch("l_pT", &l_pT, "l_pT/D");
	l_k_tree->Branch("l_eta", &l_eta, "l_eta/D");
	l_k_tree->Branch("l_phi", &l_phi, "l_phi/D");


    //----- FILL TREE -----//
    while (pythia.next()) {

        std::vector<Particle> electrons;
        std::vector<Particle> positrons;
        std::vector<Particle> muons;
        std::vector<Particle> antimuons;
        std::vector<Particle> taus;
        std::vector<Particle> antitaus;

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
                int motherIndex = pythia.event[i].mother1(); // Parent index
                if (motherIndex > 0) {
                    Particle &mother = pythia.event[motherIndex];

                    // Check if the parent is a tau lepton
                    if (mother.id() == 15) { // Tau
                        taus.push_back(mother);
                    } else if (mother.id() == -15) { // Anti-tau
                        antitaus.push_back(mother);
                    }
                }
            }
        }

        // e- e+ invariant mass m_ee
        for (auto &e : electrons) {
            for (auto &e_bar : positrons) {
            	// Calculate inv mass
                double px_tot = e.px() + e_bar.px();
                double py_tot = e.py() + e_bar.py();
                double pz_tot = e.pz() + e_bar.pz();
                double e_tot = e.e() + e_bar.e();
                m_ee = sqrt(e_tot * e_tot - px_tot * px_tot - py_tot * py_tot - pz_tot * pz_tot);
                m_ll = m_ee;
                e_tree->Fill();
                l_tree->Fill();
                // Calculate e kinematics
                e_px = e.px();
                e_py = e.py();
                e_pz = e.pz();
                e_e = e.e();
                e_pT = e.pT();
                e_phi = e.phi();
                e_eta = e.eta();
                e_k_tree->Fill();
                // Update lepton kinematics
                l_px = e_px;
				l_py = e_py;
				l_pz = e_pz;
				l_e = e_e;
				l_pT = e_pT;
				l_phi = e_phi;
				l_eta = e_eta;
				l_k_tree->Fill();
            }
        }

        // mu- mu+ invariant mass m_mm
        for (auto &mu : muons) {
            for (auto &mu_bar : antimuons) {
                double px_tot = mu.px() + mu_bar.px();
                double py_tot = mu.py() + mu_bar.py();
                double pz_tot = mu.pz() + mu_bar.pz();
                double e_tot = mu.e() + mu_bar.e();
                m_mm = sqrt(e_tot * e_tot - px_tot * px_tot - py_tot * py_tot - pz_tot * pz_tot);
                m_ll = m_mm;
                mu_tree->Fill();
                l_tree->Fill();
                // Calculate mu kinematics
				mu_px = mu.px();
				mu_py = mu.py();
				mu_pz = mu.pz();
				mu_e = mu.e();
				mu_pT = mu.pT();
				mu_phi = mu.phi();
				mu_eta = mu.eta();
				mu_k_tree->Fill();
				// Update lepton kinematics
				l_px = mu_px;
				l_py = mu_py;
				l_pz = mu_pz;
				l_e = mu_e;
				l_pT = mu_pT;
				l_phi = mu_phi;
				l_eta = mu_eta;
				l_k_tree->Fill();

            }
        }
        
        // tau- tau+ invariant mass m_tt
        for (auto &tau : taus) {
            for (auto &tau_bar : antitaus) {
                double px_tot = tau.px() + tau_bar.px();
                double py_tot = tau.py() + tau_bar.py();
                double pz_tot = tau.pz() + tau_bar.pz();
                double e_tot = tau.e() + tau_bar.e();
                m_tt = sqrt(e_tot * e_tot - px_tot * px_tot - py_tot * py_tot - pz_tot * pz_tot);
                m_ll = m_tt;
                tau_tree->Fill();
                l_tree->Fill();
                // Calculate tau kinematics
				tau_px = tau.px();
				tau_py = tau.py();
				tau_pz = tau.pz();
				tau_e = tau.e();
				tau_pT = tau.pT();
				tau_phi = tau.phi();
				tau_eta = tau.eta();
				tau_k_tree->Fill();
				// Update lepton kinematics
				l_px = tau_px;
				l_py = tau_py;
				l_pz = tau_pz;
				l_e = tau_e;
				l_pT = tau_pT;
				l_phi = tau_phi;
				l_eta = tau_eta;
				l_k_tree->Fill();

            }
        }
    }

    output->Write();
    output->Close();
    std::cout << "The file run.root has been created.\n";
    return 0;
}
