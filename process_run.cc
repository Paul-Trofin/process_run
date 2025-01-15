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
	
	// COUNT THE NUMBER OF PARTICLES
	int Ne = 0, Ne_bar = 0, Nmu = 0, Nmu_bar = 0, Nta = 0, Nta_bar = 0, Nl = 0, Nl_bar = 0, Nha = 0;
	int moNe = 0, moNe_bar = 0, moNmu = 0, moNmu_bar = 0, moNta = 0, moNta_bar = 0, moNl = 0, moNl_bar = 0, moNha = 0;
	int gmoNe = 0, gmoNe_bar = 0, gmoNmu = 0, gmoNmu_bar = 0, gmoNta = 0, gmoNta_bar = 0, gmoNl = 0, gmoNl_bar = 0, gmoNha = 0;
	


    //----- FILL TREE -----//
    while (pythia.next()) {

        std::vector<Particle> electrons;
        std::vector<Particle> positrons;
        std::vector<Particle> muons;
        std::vector<Particle> antimuons;
        std::vector<Particle> taus;
        std::vector<Particle> antitaus;
        std::vector<Particle> hadrons;
        std::vector<Particle> zDecays;

        // EACH PARTICLE IN CURRENT EVENT
        for (int i = 0; i < pythia.event.size(); ++i) {
            if (pythia.event[i].isFinal()) {
                if (pythia.event[i].id() == 11) { // Electron
                    electrons.push_back(pythia.event[i]);
                    Ne ++;
                }
                if (pythia.event[i].id() == -11) { // Positron
                    positrons.push_back(pythia.event[i]);
                    Ne_bar ++;
                }
                if (pythia.event[i].id() == 13) { // Muon
                    muons.push_back(pythia.event[i]);
                    Nmu ++;
                }
                if (pythia.event[i].id() == -13) { // Antimuon
                    antimuons.push_back(pythia.event[i]);
                    Nmu_bar ++;
                }
                if (pythia.event[i].id() == 15) { // Antimuon
                    taus.push_back(pythia.event[i]);
                    Nta ++;
                }
                if (pythia.event[i].id() == -15) { // Antimuon
                    antitaus.push_back(pythia.event[i]);
                    Nta_bar ++;
                }
                if (pythia.event[i].isHadron()) { // Hadrons
                	hadrons.push_back(pythia.event[i]);
                	Nha ++;
                }
                // MOTHER-1 ANALYSIS
                int motherIndex = pythia.event[i].mother1(); // MOTHER-1
                if (motherIndex > 0) {
                    Particle &mother = pythia.event[motherIndex];
                    
                    // Mother = Electron
                    if (mother.id() == 11) { // Electron
                        electrons.push_back(mother);
                        moNe ++;
                    } else if (mother.id() == -11) { // Positron
                        positrons.push_back(mother);
                        moNe_bar ++;
                    }
                    
                    // Mother = Muon
                    if (mother.id() == 13) { // Muon
                        muons.push_back(mother);
                        moNmu ++;
                    } else if (mother.id() == -13) { // Anti-muon
                        antimuons.push_back(mother);
                        moNmu_bar ++;
                    }
                    
                    // Mother = Tau
                    if (mother.id() == 15) { // Tau
                        taus.push_back(mother);
                        moNta ++;
                    } else if (mother.id() == -15) { // Anti-tau
                        antitaus.push_back(mother);
                        moNta_bar ++;
                    }
                    
                    if (mother.id() == 23) { // Z boson
				        zDecays.push_back(pythia.event[i]);
				    }
				}
				// MOTHER-2 ANALYSIS "GRANMA :)"
                int motherIndex2 = pythia.event[i].mother2(); // MOTHER-2
                if (motherIndex2 > 0) {
                	Particle &mother2 = pythia.event[motherIndex2];
                    if (mother2.id() == 23) { // Z boson
				        zDecays.push_back(pythia.event[i]);
				    }
				    if (mother2.id() == 11) { // electron
				        electrons.push_back(pythia.event[i]);
				        gmoNe ++;
				    }
				    if (mother2.id() == -11) { // positron
				        positrons.push_back(pythia.event[i]);
				        gmoNe_bar ++;
				    }
				    if (mother2.id() == 13) { // muon
				        muons.push_back(pythia.event[i]);
				        gmoNmu ++;
				    }
				    if (mother2.id() == -13) { // antimuon
				        antimuons.push_back(pythia.event[i]);
				        gmoNmu_bar ++;
				    }
				    if (mother2.id() == 15) { // tau
				        taus.push_back(pythia.event[i]);
				        gmoNta ++;
				    }
				    if (mother2.id() == -15) { // antitau
				        antitaus.push_back(pythia.event[i]);
				        gmoNta_bar ++;
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

    //----- NUMBER OF PARTICLES -----//
    Nl = Ne + Nmu + Nta;
    Nl_bar = Ne_bar + Nmu_bar + Nta_bar;
    moNl = moNe + moNmu + moNta;
    moNl_bar = moNe_bar + moNmu_bar + moNta_bar;
    gmoNl = gmoNe + gmoNmu + gmoNta;
    gmoNl_bar = gmoNe_bar + gmoNmu_bar + gmoNta_bar;
    
    std::cout << "//----- NUMBER OF PARTICLES -----//\n";
    std::cout << "//----- LEPTONS -----//\n";
    std::cout << "//----- Total leptons: " << Nl + Nl_bar + moNl + moNl_bar + gmoNl + gmoNl_bar << "\n";
    std::cout << "//----- Daughter -----//\n";
    std::cout << "//----- l-: " << Nl << "\n";
    std::cout << "//----- l+: " << Nl_bar << "\n";
    std::cout << "//----- e-: " << Ne << "\n";
    std::cout << "//----- e+: " << Ne_bar << "\n";
    std::cout << "//----- mu-: " << Nmu << "\n";
    std::cout << "//----- mu+: " << Nmu_bar << "\n";
    std::cout << "//----- ta-: " << Nta << "\n";
    std::cout << "//----- ta+: " << Nta_bar << "\n";
    std::cout << "//----- Mother -----//\n";
    std::cout << "//----- l-: " << moNl << "\n";
    std::cout << "//----- l+: " << moNl_bar << "\n";
    std::cout << "//----- e-: " << moNe << "\n";
    std::cout << "//----- e+: " << moNe_bar << "\n";
    std::cout << "//----- mu-: " << moNmu << "\n";
    std::cout << "//----- mu+: " << moNmu_bar << "\n";
    std::cout << "//----- ta-: " << moNta << "\n";
    std::cout << "//----- ta+: " << moNta_bar << "\n";
    std::cout << "//----- GrandMother -----//\n";
    std::cout << "//----- l-: " << gmoNl << "\n";
    std::cout << "//----- l+: " << gmoNl_bar << "\n";
    std::cout << "//----- e-: " << gmoNe << "\n";
    std::cout << "//----- e+: " << gmoNe_bar << "\n";
    std::cout << "//----- mu-: " << gmoNmu << "\n";
    std::cout << "//----- mu+: " << gmoNmu_bar << "\n";
    std::cout << "//----- ta-: " << gmoNta << "\n";
    std::cout << "//----- ta+: " << gmoNta_bar << "\n";
    std::cout << "//----- HADRONS -----//\n";
    std::cout << "//----- Total hadrons: " << Nha << "\n";
    std::cout << "//-----------------------------//";

    output->Write();
    output->Close();
    std::cout << "The file run.root has been created.\n";
    return 0;
}
