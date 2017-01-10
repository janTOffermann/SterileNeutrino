//
//  DetectorSim3.c
//
//
//  Created by Jan Offermann on 11/17/16.
//
//

#include <iostream>
#include <fstream>
#include <string>
#include <cctype>
#include <vector>
#include <algorithm>
#include <random>

#include "TMath.h"
#include "TStyle.h"
#include "TF1.h"
#include "TSystem.h"
#include "TString.h"
#include "TError.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TPad.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TFitResult.h"

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"

// Global variables

// distances from source to multi-channel plates for recoil ion (R), and Auger electron (e),
// and radius of apparatus
const Double_t mcp_R = 0.5; // [m]
const Double_t mcp_e = -3.; // [m]
const Double_t radius = 3.; // [m], can simply decrease radius and increase Bz if needed


// magnetic field [T]
const Double_t Bx = 0.;
const Double_t By = 0.;
const Double_t Bz = 1.E-4; // [T], ~1 Gauss, paper says a "few Gauss"

// electric field [V/m]
const Double_t Ex = 0.;
const Double_t Ey = 0.;
const Double_t Ez = 500.; // 200 [V/m] in simulation done in paper, see Fig 3.2

// Physical constants
const Double_t c   = 299792458; // speed of light [m/s]
const Double_t q   = 1.6022E-19; // elementary charge, 1.6022E-19 [C]
const Double_t m_e = 510.9989; // electron mass [keV]
const Double_t m_Cs = 1.219376542142855E8; // Cesium mass [keV]
const Double_t Q = 352.; //[keV], from CJM_final.docx
const Double_t Q1 = Q - 34.4 - 0.12; // - (N shell x-ray) - (N shell Auger)
const Double_t m_Xe = m_Cs - Q;

const Double_t m_proton = 938272; // [keV]
const Double_t m_neutron = 939565.4; // [keV]
const Double_t m_v_sterile = 10.; // sterile neutrino mass [keV]

const Double_t kgms_to_keV = 1.8711573624779732066468E24; // convert [kg m/s] to [keV]
const Double_t keV_to_kg   = 1.7826619E-33; // convert [keV] to [kg]



/* Tabulated gamma ray energies, for gammas emmited as excited Xe131
 * transitions towards the ground state (Auger electrons will also be released).
 *
 * As per CJM_final.docx,
 *
 * L-shell -> 29.4, 29.8: these will be rejected later on
 *
 * M-shell -> 33.56,33.62
 *
 * N-shell -> 33.4
 *
 * The values below are taken from
 * http://www.nndc.bnl.gov/chart/decaysearchdirect.jsp?nuc=131CS&unc=nds
 * (more precise than those in CJM_final.docx).
 */
std::vector<Double_t> E_g_vals = {29.461,29.782,33.562,33.624,34.419}; // [keV]


/*
 * The Event object stores detector information, which is used
 * to reconstruct the neutrino square-mass histogram. It only stores
 * measured/known quantities.
 *
 */
class Event : public TObject {
    
public:
    
    // Note: For ion and electron, x and y give position.
    // For the gamma, RLorentzVector gives direction. We don't actually know the energy,
    // only have  some idea of range (our detector isn't sensitive enough to distinguish
    // between all individual energies), so we use the direction and apply a couple different
    // energies based on the actual one.
    Int_t type; // 0 = gamma, 1 = electron, 2 = ion
    Double_t x;
    Double_t y;
    Double_t t;
    Double_t m;
    Double_t w;
    Double_t phi;
    TLorentzVector gamma;
    
    ClassDef(Event,1)
};

ClassImp(Event)

/*
 * The RealEvent object stores TLorentzVector objects
 * for all complete events - this includes measured
 * TLorenzVectors, but unlike the Event object, information
 * from this object will not be susceptible to issues arising
 * from pile-up, or unknown gamma ray energies.
 * For analysis purposes.
 */
class RealEvent : public TObject {
    
public:
    
    TLorentzVector lor_g;
    TLorentzVector lor_R_m_c;
    TLorentzVector lor_e_m_c;
    TLorentzVector lor_v_m_c;
    
    ClassDef(RealEvent,1)
};

ClassImp(RealEvent)


void InitStyle()
{
    gStyle->SetPaperSize(TStyle::kUSLetter);
    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42,"xy");
    gStyle->SetTitleFont(42,"xy");
    gStyle->SetLabelSize(0.05,"xy");
    gStyle->SetTextSize(1.5);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleOffset(1.75,"Y");
    gStyle->SetTitleOffset(1.25,"X");
    gStyle->SetTitleAlign(13);
    gStyle->SetTitleSize(0.05);
    gStyle->SetTitleSize(0.025,"t");
    gStyle->SetCanvasColor(kWhite);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetPadColor(kWhite);
    gStyle->SetPadTopMargin(0.125);
    gStyle->SetPadLeftMargin(0.125);
    gStyle->SetPadBottomMargin(0.125);
    gStyle->SetPadRightMargin(0.10);
    gStyle->SetFillColor(kWhite);
    gStyle->SetTitleFillColor(kWhite);
    gStyle->SetStatX(0.9);
}

TH1D *DefineHist1D(
                   const TString     &h1d_name,
                   const TString     &xAxis,
                   Int_t        nbins,
                   Double_t     xMin,
                   Double_t     xMax)
{
    TH1D *h1d = new TH1D(h1d_name,h1d_name,nbins,xMin,xMax);
    
    h1d->SetBit(TH1::kNoTitle);
    h1d->SetDirectory(0);
    h1d->SetMarkerSize(0.5);
    h1d->SetMarkerStyle(20);
    h1d->SetLineStyle(1);
    h1d->SetMarkerColor(kBlack);
    
    h1d->Draw("");
    h1d->GetXaxis()->SetNoExponent();
    h1d->GetXaxis()->SetTitle(xAxis);
    
    return h1d;
}

static bool sort_using_less_than(double u, double v) {
    return u < v;
}

// cross product function
std::vector<Double_t> crossProduct(Double_t v1, Double_t v2, Double_t v3, Double_t b1, Double_t b2, Double_t b3) {
    
    const Double_t comp1 = v2 * b3 - v3 * b2;
    const Double_t comp2 = v3 * b1 - v1 * b3;
    const Double_t comp3 = v1 * b2 - v2 * b1;
    
    std::vector<Double_t> vec = {comp1, comp2, comp3};
    
    return vec;
}

// simulate the event
void Sim(UInt_t limit, Double_t mixAngle = 0.1, Bool_t no_l_shell = kTRUE, Bool_t ERRORS = kFALSE) {
    
    std::random_device rd;     // only used once to initialise (seed) engine
    //    std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
    std::mt19937 rng(1);    // random-number engine used (Mersenne-Twister in this case)
    
    std::uniform_real_distribution<Double_t> uni(0.,1.); // guaranteed unbiased
    std::normal_distribution<Double_t> electron_gauss(0.,4.E-5);
    std::normal_distribution<Double_t> ion_gauss(0.,4.E-5);
    std::normal_distribution<Double_t> TOF_gauss(0.,2.E-10);
    std::uniform_int_distribution<Int_t> uni_int(0,1);
    
    Double_t m_v;
    
    // charge of ion, will depend on number of Auger electrons released
    Double_t q_ion = 0.;
    
    // used for cross products (to calculate accelerations if needed)
    std::vector<Double_t> cP;
    
    Double_t total_auger_energy; // changes per run
    
    UInt_t gamma_type; // 0 = L shell, 1 = M shell, 2 = N shell
    
    UInt_t Nelectrons; // number of Auger electrons, changes per run
    
    Double_t global_time = 0.; // "global time" (used to determine times of all events)
    
    Double_t Natoms = 1.E8;
    //    Double_t lambda = 8.28005E-7; // taken from half life of 9.689 days
    Double_t lambda = 8.28005; // increased to force overlaps
    Double_t rate = lambda * Natoms;
    Double_t time_window = (Double_t)Natoms / rate;
    
    // create the times of decays
    std::vector<Double_t> decay_times;
    
    // for now, randomly distribute over time window using flat distribution
    for (UInt_t i = 0; i < limit; i++) {
        
        decay_times.push_back(uni(rng) * time_window);
    }
    
    std::sort(decay_times.begin(), decay_times.end(), sort_using_less_than);
    
    TFile f("sim.root","recreate");
    
    // TTree used to store detector info (for reconstructing)
    TTree *t1 = new TTree("tree1", "tree1");
    Event *event = new Event();
    
    t1->Branch("event",&event);
    
    // TTree used to store TLorentzVectors (for more plotting)
    TTree *t2 = new TTree("tree2", "tree2");
    RealEvent *revent = new RealEvent();
    
    t2->Branch("revent",&revent);
    
    TH1D* histogram = new TH1D("", "m_{#nu}^{2};[keV^{2}];", 400, -2E2, 2E2);
    TH1D* histogram2 = new TH1D("ion TOF", "Recoil Ion Flight Times;[s];", 500, 1.5E-5, 7.5E-5);
    TH1D* histogram3 = new TH1D("Auger TOF", "Auger Electron Flight Times;[s];", 500, 1.5E-7, 4.E-7);
    
    //    const Double_t m_v_original = m_v;
    
    for (UInt_t counter = 0; counter < limit; counter ++) {
        
        if (uni(rng) > mixAngle) m_v = 0.;
        else m_v = m_v_sterile;
        
        // Cs131 Lorentz vector
        const TLorentzVector lor_Cs(0.,0.,0.,m_Cs);
        
        // Xe131: Setup the Lorentz vector which is formed by the decay of Cs
        // through neutrino emission
        TLorentzVector lor_R;
        {
            const Double_t m_R_ex = m_Cs - m_v - Q1;
            
            const Double_t E_v = (Q1 * Q1 + 2. * Q1 * m_R_ex + 2. * Q1 * m_v + 2. * m_R_ex * m_v)/(2. * (Q1 + m_R_ex + m_v));
            const Double_t p_v = sqrt(E_v * E_v - m_v * m_v);
            
            const Double_t E_R = sqrt(m_R_ex * m_R_ex + p_v * p_v);
            
            // now get direction by choosing a random direction on a sphere
            const Double_t theta_R = acos(2.*uni(rng) - 1.);
            const Double_t phi_R   = 2. * M_PI * uni(rng);
            
            const Double_t px_R = - p_v * cos(phi_R) * sin(theta_R);
            const Double_t py_R = - p_v * sin(phi_R) * sin(theta_R);
            const Double_t pz_R = - p_v * cos(theta_R);
            lor_R.SetPxPyPzE(px_R,py_R,pz_R,E_R);
        }
        
        TLorentzVector lor_vt;
        
        // neutrino Lorentz vector: same 3-momentum as Xe but opposite
        {
            const Double_t E_v = (Q1 * Q1 + 2. * Q1 * m_Xe + 2. * Q1 * m_v + 2. * m_Xe * m_v)/(2. * (Q1 + m_Xe + m_v));
            
            const Double_t px_v = -lor_R.Px();
            const Double_t py_v = -lor_R.Py();
            const Double_t pz_v = -lor_R.Pz();
            lor_vt.SetPxPyPzE(px_v,py_v,pz_v,E_v);
        }
        
        //        std::cout << std::endl;
        //        std::cout << "Check 0." <<  std::endl;
        //        std::cout << "m^2  = " << lor_vt.M2() << std::endl;
        //        std::cout << std::endl;
        //
        //        lor_vt = lor_Cs - lor_R;
        //        std::cout << std::endl;
        //        std::cout << "Check 1." <<  std::endl;
        //        std::cout << "m^2  = " << lor_vt.M2() << std::endl;
        //        std::cout << std::endl;
        
        // Next step: Gamma ray
        // We'll need to boost into recoil nucleus's rest frame, we know possible gamma energies there
        
        // gamma ray Lorentz vector in lab frame
        TLorentzVector lor_g;
        {
            // pick gamma energy in recoil nucleus frame from tabulated values, determine # of Augers accordingly
            // using x-ray probabilities given on page 12-13 of CJM_final.docx
            
            Double_t rand = 100. * uni(rng);
            Int_t rand_int = uni_int(rng);
            
            Double_t E_g;
            
            if (rand <= 81.76) {
                
                E_g = E_g_vals.at(0 + rand_int);
                gamma_type = 0;
            }
            else if (rand <= 97.05) {
                
                E_g = E_g_vals.at(2 + rand_int);
                gamma_type = 1;
            }
            else {
                
                E_g = E_g_vals.at(4);
                gamma_type = 2;
            }
            
            // in case of debug, don't get x-rays from L-shell since we exclude those events later on anyway
            if (no_l_shell) {
                
                if (rand <= 83.9) {
                    
                    E_g = E_g_vals.at(2 + rand_int);
                    gamma_type = 1;
                }
                else {
                    
                    E_g = E_g_vals.at(4);
                    gamma_type = 2;
                }
            }
            
            if (E_g < 30.) Nelectrons = 3; // actual number not given in CJM_final.docx, but presumably >= 3 Auger electrons
            
            else if (E_g < 34.) {
                
                // 2-3 Auger electrons
                if (uni(rng) >= 0.5) Nelectrons = 2;
                else Nelectrons = 3;
            }
            
            else {
                
                // 1-2 Auger electrons
                if (uni(rng) >= 0.5) Nelectrons = 1;
                else Nelectrons = 2;
            }
            
            // Since gamma mass is zero, momentum magnitude == energy
            const Double_t theta_g = acos(2.*uni(rng) - 1.);
            const Double_t phi_g   = 2. * M_PI * uni(rng);
            
            const Double_t px_g = E_g * cos(phi_g) * sin(theta_g);
            const Double_t py_g = E_g * sin(phi_g) * sin(theta_g);
            const Double_t pz_g = E_g * cos(theta_g);
            lor_g.SetPxPyPzE(px_g,py_g,pz_g,E_g);
            
            // now apply lorentz transform to lab frame - we say the frame is moving (opposite direction of motion of recoil nucleus in the lab frame)
            const Double_t vx = -lor_R.Px() / lor_R.M();
            const Double_t vy = -lor_R.Py() / lor_R.M();
            const Double_t vz = -lor_R.Pz() / lor_R.M(); // please note a factor of 1/c built-in since we use [keV] for p and m
            
            const TLorentzRotation boost_g(vx,vy,vz);
            lor_g.Transform(boost_g);
            
            // set total Auger energy, and divide it up amongst Augers (to be applied later)
            total_auger_energy = Q - Q1 - E_g;
        }
        
        // Apply conservation of momentum and energy - adjust momentum and energy of recoil nucleus
        // energy of gamma resulting in lowered mass for recoil nucleus, should follow from lor_R.M()
        lor_R -= lor_g;
        
        lor_vt = lor_Cs - lor_R - lor_g;
        //        std::cout << std::endl;
        //        std::cout << "Second neutrino check." <<  std::endl;
        //        std::cout << "m^2  = " << lor_vt.M2() << std::endl;
        //        std::cout << std::endl;
        
        // NOTE: I am not sure what precisely to set for the Auger electron energy. I do know the total
        // Auger energy, as mentioned on page 13 of CJM_final.docx, however the bands are not given.
        // I will try to use ranges given on page 13, but the below code may need modification
        // to be more physically accurate (esp. if bands are known).
        
        std::vector<Double_t> auger_energies;
        
        for (UInt_t i = 0; i < Nelectrons; i++) {
            
            Double_t rand = uni(rng);
            
            Double_t individual_energy;
            
            // assign individual energy based on which x-ray was emitted
            
            // L - shell
            if (gamma_type == 0) {
                
                if (uni(rng) <= 0.14) individual_energy = 0.05 + rand * 0.35;
                else individual_energy = 3. + rand;
            }
            else if (gamma_type == 1) individual_energy = 0.3 + rand * 0.1;
            
            else individual_energy = 0.05 + rand * 0.07;
            
            //            Double_t individual_energy = 0.05 + rand * 0.25;
            
            //            std::cout << "i = " << i << std::endl;
            //            std::cout << "\t" << individual_energy * 1000. << std::endl;
            
            if (individual_energy > total_auger_energy) {
                
                // if it falls within allowed range, just assign all remaining energy to Auger
                if (gamma_type != 1 && total_auger_energy >= 0.05) {
                    
                    individual_energy = total_auger_energy;
                }
                
                else if (gamma_type == 1 && total_auger_energy >= 0.3) {
                    
                    individual_energy = total_auger_energy;
                }
                
                // if remaining energy is below (supposed) lower bounds, add it to previous Auger
                
                else {
                    auger_energies.back() += total_auger_energy;
                    Nelectrons = i;
                    break;
                }
            }
            
            auger_energies.push_back(individual_energy);
            
            total_auger_energy = total_auger_energy - individual_energy;
        }
        
        // set charge of ion now that we certainly know Nelectrons
        q_ion = (Double_t)(Nelectrons) * q;
        
        // Next, release  Auger electron(s). Again, we have to boost frames.
        std::vector<TLorentzVector> lor_e;
        
        for (UInt_t i = 0; i < Nelectrons; i++) {
            
            // in boosted frame
            const Double_t E_e = auger_energies.at(i) + m_e; // need to fix
            const Double_t p_e = sqrt(E_e * E_e - m_e * m_e);
            
            const Double_t theta_e = acos(2.*uni(rng) - 1.);
            const Double_t phi_e   = 2. * M_PI * uni(rng);
            
            const Double_t px_e = p_e * cos(phi_e) * sin(theta_e);
            const Double_t py_e = p_e * sin(phi_e) * sin(theta_e);
            const Double_t pz_e = p_e * cos(theta_e);
            TLorentzVector lor_e_single;
            
            lor_e_single.SetPxPyPzE(px_e,py_e,pz_e,E_e);
            
            // now apply lorentz transform to lab frame - we say the frame is moving (opposite direction of motion of recoil nucleus in the lab frame)
            const Double_t vx = -lor_R.Px() / lor_R.M();
            const Double_t vy = -lor_R.Py() / lor_R.M();
            const Double_t vz = -lor_R.Pz() / lor_R.M(); // please note a factor of 1/c built-in since we use keV for p and m
            TLorentzRotation boost_e(vx,vy,vz);
            lor_e_single.Transform(boost_e);
            
            lor_e.push_back(lor_e_single);
            
            // Apply conservation of momentum and energy - adjust momentum and energy of recoil nucleus
            lor_R -= lor_e_single;
        }
        
        
        // We now have a recoil ion, Auger electron(s) and gamma (we will not simulate where neutrino goes)
        // and we can calculate the Lorentz vector for the neutrino
        TLorentzVector lor_v = lor_Cs-lor_R-lor_g;
        for (UInt_t i = 0; i < lor_e.size(); i++) lor_v = lor_v - lor_e.at(i);
        
        //        std::cout << std::endl;
        //        std::cout << "Third neutrino check." <<  std::endl;
        //        std::cout << "m^2  = " << lor_v.M2() << std::endl;
        //        std::cout << std::endl;
        
        // At this point, K-capture is complete.
        
        /*
         * At this point, K-capture is complete.
         * Now, we need to propagate the recoil nucleus and the Auger electron(s) through the apparatus, into the detectors.
         * We will use discrete time steps - some code below currently takes advantage of the constant fields to get better results
         * but I have left the stepping architecture so that this can be easily modified for more complicated field maps.
         *
         * As an important note, some lines are commented out for the recoil ion and electron because I use constant fields,
         * so I can easily get analytical expressions for their motions. However this is used more extensively for the electron,
         * as seemingly due to some issues of numerical precision, this analytical approach gives incorrect results for the ion -
         * the errors are relatively small but the m^2 calculation is very sensitive to these. Decreasing B-field to 1.E-6 [T]
         * seems to remove these errors.
         *
         */
        
        Long64_t timeCount = 0; // counts time in nanoseconds
        
        const Double_t m_R_kg = lor_R.M() * keV_to_kg;
        const Double_t m_e_kg = m_e * keV_to_kg;
        
        // cyclotron frequency
        const Double_t w_R = q_ion * Bz /(m_R_kg);
        const Double_t w_e = -q * Bz /(m_e_kg);
        
        const Long64_t maxTime = 1E8; // 1/10 second (if event takes longer than this, we reject it)
        
        const Double_t timeStep_R = 1.E-11; // time step of 0.01 nanoseconds
        const Double_t timeStep_e = 1.E-11; // time step of 0.01 nanoseconds (previously was 0.001 nanoseconds)
        
        // actual initial ejection angles (for simulation, not for reconstruction since we don't directly observe these)
        const Double_t phi_R_th = atan2(lor_R.Py(),lor_R.Px());
        std::vector<Double_t> phi_e_th;
        
        for (UInt_t i =0; i < Nelectrons; i++) phi_e_th.push_back(atan2(lor_e.at(i).Py(),lor_e.at(i).Px()));
        
        // set initial positions for recoil ion, get initial velocities [m/s] from the momentum
        Double_t x_R = 0.;
        Double_t y_R = 0.;
        Double_t z_R = 0.;
        Double_t angle_R = phi_R_th;
        
        Double_t vx_R = lor_R.Px()/lor_R.M() * c;
        Double_t vy_R = lor_R.Py()/lor_R.M() * c;
        Double_t vz_R = lor_R.Pz()/lor_R.M() * c;
        
        // radius of path when projected onto detector plane (circular since fields are constant)
        const Double_t rad_R = abs(sqrt(vx_R * vx_R + vy_R * vy_R) / w_R);
        
        // first, propagate the recoil nucleus
        Bool_t impact = kFALSE;
        
        // take advantage of constant field here
        const Double_t az_R = Ez * q_ion / m_R_kg;
        
        Double_t TOFR_real = (-vz_R + sqrt(vz_R * vz_R + 2. * az_R * mcp_R))/az_R;
        
        while(!impact) {
            
            timeCount ++;
            
            x_R += vx_R * timeStep_R;
            y_R += vy_R * timeStep_R;
            z_R += vz_R * timeStep_R;
            
            //            z_R += 0.5 * az_R * timeStep_R * timeStep_R + vz_R * timeStep_R;
            
            angle_R += w_R * timeStep_R;
            
            // get new velocities via acceleration
            const std::vector<Double_t> cP = crossProduct(vx_R, vy_R, vz_R, Bx, By, Bz);
            
            const Double_t ax_R = (Ex * q_ion + q_ion * cP[0])/m_R_kg;
            const Double_t ay_R = (Ey * q_ion + q_ion * cP[1])/m_R_kg;
            //            const Double_t az_R = (Ez * q_ion + q_ion * cP[2])/m_R_kg;
            
            vx_R += ax_R * timeStep_R;
            vy_R += ay_R * timeStep_R;
            vz_R += az_R * timeStep_R;
            
            // collide with wrong detector -> not captured
            if (z_R < mcp_e) {
                timeCount = -1;
                break;
            }
            
            // collide with sides of chamber -> not captured
            if (x_R * x_R + y_R * y_R >= radius * radius) {
                timeCount = -1;
                break;
            }
            
            if (z_R < mcp_R && timeCount <= maxTime) continue;
            
            impact = kTRUE;
            
            if (timeCount > maxTime) {
                timeCount = -1;
            }
            
            break;
        }
        
        const Double_t TOFR = (Double_t)(timeCount) * timeStep_R;
        
        // as we have possibly multiple Augers, we use vectors to keep their info
        //        std::vector<Double_t> x_e_vec;
        //        std::vector<Double_t> y_e_vec;
        std::vector<Double_t> z_e_vec;
        std::vector<Double_t> vz_e_vec;
        std::vector<Double_t> TOFE_vec;
        std::vector<Double_t> rad_e_vec;
        
        // again, take advantage of constant field
        const Double_t az_e = (- Ez * q)/m_e_kg;
        
        for (UInt_t i = 0; i < Nelectrons; i++) {
            
            // reset to run for electron
            impact = kFALSE;
            timeCount = 0;
            
            Double_t x_e = 0.;
            Double_t y_e = 0.;
            Double_t z_e = 0.;
            
            Double_t vx_e = lor_e.at(i).Px()/lor_e.at(i).M() * c;
            Double_t vy_e = lor_e.at(i).Py()/lor_e.at(i).M() * c;
            Double_t vz_e = lor_e.at(i).Pz()/lor_e.at(i).M() * c;
            
            const Double_t rad_e = abs(sqrt(vx_e * vx_e + vy_e * vy_e) / w_e);
            
            Double_t TOFE_real = (-vz_e + sqrt(vz_e * vz_e + 2. * az_e * mcp_e))/az_e;
            
            while(!impact) {
                
                timeCount ++;
                
                //            x_e += vx_e * timeStep_e;
                //            y_e += vy_e * timeStep_e;
                z_e += vz_e * timeStep_e;
                //                z_e += 0.5 * az_e * timeStep_e * timeStep_e + vz_e * timeStep_e;
                
                
                // get new velocities via acceleration
                //            const std::vector<Double_t> cP = crossProduct(vx_e, vy_e, vz_e, Bx, By, Bz);
                
                //            const Double_t ax_e = (- Ex * q - q * cP[0])/m_e_kg;
                //            const Double_t ay_e = (- Ey * q - q * cP[1])/m_e_kg;
                
                //            vx_e += ax_e * timeStep_e;
                //            vy_e += ay_e * timeStep_e;
                vz_e += az_e * timeStep_e;
                
                
                // collide with wrong detector -> not captured
                if (z_e > mcp_R) {
                    timeCount = -1;
                    
                    //                std::cout << "\t E =  " << lor_e.at(i).E() - lor_e.at(i).M() << std::endl;
                    //                std::cout << "\t vz = " << lor_e.at(i).Pz()/lor_e.at(i).M() * c << std::endl;
                    
                    break;
                }
                
                // collide with sides of chamber -> not captured
                if (2 * rad_e >= radius) {
                    timeCount = -1;
                    break;
                }
                
                if (z_e > mcp_e && timeCount <= maxTime) continue;
                
                impact = kTRUE;
                
                if (timeCount > maxTime) {
                    
                    timeCount = -1;
                    
                }
                break;
            }
            
            const Double_t TOFE = (Double_t)(timeCount) * timeStep_e;
            
            //            x_e_vec.push_back(x_e);
            //            y_e_vec.push_back(y_e);
            z_e_vec.push_back(z_e);
            vz_e_vec.push_back(vz_e);
            TOFE_vec.push_back(TOFE);
            rad_e_vec.push_back(rad_e);
        }
        
        // correct time of flight, as we likely overshot detector but can easily fix that
        // as fields are constant
        std::vector<Double_t> z_diff_e;
        std::vector<Double_t> kappa_e;
        std::vector<Double_t> delta_t_e;
        std::vector<Double_t> TOFE_c;
        
        const Double_t vz_R_p = vz_R - az_R * timeStep_R;
        const Double_t z_diff_R = z_R - mcp_R;
        const Double_t kappa_R = 0.5 * az_R * timeStep_R * timeStep_R + vz_R_p * timeStep_R - z_diff_R;
        const Double_t delta_t_R = abs((-vz_R_p + sqrt(vz_R_p * vz_R_p - 2. * az_R * kappa_R))/az_R);
        Double_t TOFR_c = TOFR - timeStep_R + delta_t_R;
        
        if (ERRORS) TOFR_c += TOF_gauss(rng);
        
        for (UInt_t i = 0; i < Nelectrons; i++) {
            
            const Double_t vz_e_p = vz_e_vec.at(i) - az_e * timeStep_R;
            z_diff_e.push_back(abs(z_e_vec.at(i) - mcp_e));
            kappa_e.push_back(0.5 * abs(az_e) * timeStep_e * timeStep_e + abs(vz_e_p) * timeStep_e - z_diff_e.at(i));
            delta_t_e.push_back(abs((-abs(vz_e_p) + sqrt(vz_e_p * vz_e_p - 2. * abs(az_e) * kappa_e.at(i)))/(abs(az_e))));
            Double_t TOFE_c_single = TOFE_vec.at(i) - timeStep_e + delta_t_e.at(i);
            if (ERRORS) TOFE_c_single += TOF_gauss(rng);
            TOFE_c.push_back(TOFE_c_single);
        }
        
        
        // get x and y using cyclotron frequency and time of flight
        
        const Double_t x_R_2 = rad_R * (-sin(phi_R_th) + sin(abs(w_R) * (TOFR_c) + phi_R_th));
        const Double_t y_R_2 = rad_R * ( cos(phi_R_th) - cos(abs(w_R) * (TOFR_c) + phi_R_th));
        
        // use analytical work for x and y (stepping subjects us to additional rounding errors -> avoid if possible)
        // also can include suggested errors in x and y for electron
        
        std::vector<Double_t> x_e_c;
        std::vector<Double_t> y_e_c;
        
        for (UInt_t i = 0; i < Nelectrons; i++) {
            
            Double_t x_e_c_tmp = rad_e_vec.at(i) * (-sin(phi_e_th.at(i)) + sin(abs(w_e) * TOFE_c.at(i) + phi_e_th.at(i)));
            Double_t y_e_c_tmp = rad_e_vec.at(i) * ( cos(phi_e_th.at(i)) - cos(abs(w_e) * TOFE_c.at(i) + phi_e_th.at(i)));
            
            if (ERRORS) {
                
                x_e_c_tmp += electron_gauss(rng);
                y_e_c_tmp += electron_gauss(rng);
            }
            
            x_e_c.push_back(x_e_c_tmp);
            y_e_c.push_back(y_e_c_tmp);
        }
        
        // currently, analytical work is not working for recoil ion - not sure why. Decreasing B field by factor of 100 fixes it,
        // so it is likely some kind of precision issue.
        //        const Double_t x_R_c = x_R_2;
        //        const Double_t y_R_c = y_R_2;
        
        Double_t x_R_c = x_R + rad_R * sin(abs(w_R) * TOFR_c + phi_R_th) - rad_R * sin(abs(w_R) * TOFR + phi_R_th);
        Double_t y_R_c = y_R - rad_R * cos(abs(w_R) * TOFR_c + phi_R_th) + rad_R * cos(abs(w_R) * TOFR + phi_R_th);
        
//                if (ERRORS) {
//        
//                    x_R_c += ion_gauss(rng);
//                    y_R_c += ion_gauss(rng);
//                }
        
        // if we don't capture the ion, skip this event (this shouldn't really happen)
        if (TOFR < 0.) continue;
        
        const Double_t Nturns_R_c = abs(floor(w_R * TOFR_c/(2. * M_PI)));
        
        std::vector<Double_t> Nturns_e_c;
        for (UInt_t i = 0; i < Nelectrons; i++) Nturns_e_c.push_back(abs(floor(w_e * TOFE_c.at(i)/(2. * M_PI))));
        
        const Double_t phi_R_c = -( - atan2(y_R_c, x_R_c) - (w_R * TOFR_c - 2. * M_PI * Nturns_R_c) / 2.);
        
        std::vector<Double_t> phi_e_c;
        for (UInt_t i = 0; i < Nelectrons; i++) phi_e_c.push_back(M_PI - (-atan2(y_e_c.at(i), x_e_c.at(i)) - (w_e * TOFE_c.at(i) - 2. * M_PI * Nturns_e_c.at(i)) / 2.));
        
        // for now, we'll fill tree info into here
        
        // gamma ray (presumably detected)
        event->type = 0;
        event->x = 0.;
        event->y = 0.;
        event->t = decay_times.at(counter);
        event->m = 0.;
        event->w = 0.;
        event->phi = 0.;
        event->gamma = lor_g;
        
        t1->Fill();
        
        //ion
        event->x = x_R_c;
        event->y = y_R_c;
        event->t = decay_times.at(counter) + TOFR_c;
        event->m = lor_R.M();
        event->w = w_R;
        event->phi = phi_R_c;
        event->type = 2;
        t1->Fill();
        
        //electron
        for (UInt_t i = 0; i < Nelectrons; i++) {
            
            if (TOFE_c.at(i) < 0.) continue; // don't fill the electrons we somehow miss (e.g. hit wrong detector, not captured by field sufficiently)
            event->x = x_e_c.at(i);
            event->y = y_e_c.at(i);
            event->t = decay_times.at(counter) + TOFE_c.at(i);
            event->m = lor_e.at(i).M();
            event->w = w_e;
            event->phi = phi_e_c.at(i);
            event->type = 1;
            t1->Fill();
        }
        
        // if we missed an electron, we're not interested in plotting the results for that run,
        // here we are making a histogram of the actual, complete events for comparison purposes
        
        Bool_t skip = kFALSE;
        for (UInt_t i = 0; i < Nelectrons; i++) {
            
            if (TOFE_c.at(i) < 0.) skip = kTRUE;
        }
        
        if (skip) continue;
        
        TLorentzVector lor_R_m_c;
        {
            const Double_t pmrx_c = q_ion * Bz * sqrt(x_R_c * x_R_c + y_R_c * y_R_c) * cos(phi_R_c) / (2. * abs(sin(w_R * TOFR_c / 2.))) * kgms_to_keV;
            const Double_t pmry_c = q_ion * Bz * sqrt(x_R_c * x_R_c + y_R_c * y_R_c) * sin(phi_R_c) / (2. * abs(sin(w_R * TOFR_c / 2.))) * kgms_to_keV;
            const Double_t pmrz_c = (m_R_kg * mcp_R/TOFR_c - q_ion * Ez * TOFR_c / 2.) * kgms_to_keV;
            const Double_t m_R = lor_R.M();
            lor_R_m_c.SetXYZM(pmrx_c,pmry_c,pmrz_c,m_R);
        }
        
        std::vector<TLorentzVector> lor_e_m_c;
        for (UInt_t i = 0; i < Nelectrons; i++) {
            
            TLorentzVector lor_e_m_c_single;
            const Double_t pmex_c = -(-q) * Bz * sqrt(x_e_c.at(i) * x_e_c.at(i) + y_e_c.at(i) * y_e_c.at(i)) * cos(phi_e_c.at(i)) / (2. * abs(sin(w_e * TOFE_c.at(i) / 2.))) * kgms_to_keV;
            const Double_t pmey_c = -(-q) * Bz * sqrt(x_e_c.at(i) * x_e_c.at(i) + y_e_c.at(i) * y_e_c.at(i)) * sin(phi_e_c.at(i)) / (2. * abs(sin(w_e * TOFE_c.at(i) / 2.))) * kgms_to_keV;
            const Double_t pmez_c = -(m_e_kg * (-mcp_e)/TOFE_c.at(i) + (-q) * Ez * TOFE_c.at(i) / 2.) * kgms_to_keV;
            
            lor_e_m_c_single.SetXYZM(pmex_c,pmey_c,pmez_c,m_e);
            lor_e_m_c.push_back(lor_e_m_c_single);
        }
        
        TLorentzVector lor_v_m_c = lor_Cs-lor_R_m_c-lor_g;
        for (UInt_t i = 0; i < Nelectrons; i++) lor_v_m_c = lor_v_m_c - lor_e_m_c.at(i);
        
        // for plotting, we can now fill RealEvent
        
        revent->lor_g = lor_g;
        revent->lor_R_m_c = lor_R_m_c;
        revent->lor_e_m_c = lor_e_m_c.at(0);
        revent->lor_v_m_c = lor_v_m_c;
        t2->Fill();
        
        const Double_t m_v_sq_c = lor_v_m_c.M2();
        
        // we're throwing out events that result in L-shell x-rays, so don't plot those
        if (lor_g.E() > 30.) histogram->Fill(m_v_sq_c);
        
        for (UInt_t blah = 0; blah < Nelectrons; blah++) {
            
            histogram3->Fill(TOFE_c.at(blah));
        }
        
        histogram2->Fill(TOFR_c);
        
        if (counter % 1000 != 0) continue;
        std::cout << "++ " << counter << " \\ " << m_v_sq_c << std::endl;
    }
    
    t1->Write();
    t2->Write();
    histogram->Write();
    histogram->SetDirectory(0);
    histogram->Draw();
    //    histogram2->Write();
    //    histogram2->SetDirectory(0);
    //    histogram2->Draw();
}

Double_t MyExp(Double_t *x, Double_t *par)
{
    Double_t arg = 0;
    if (par[2] != 0) arg = (x[0]-par[1])/par[2];
    
    const Double_t funcVal = par[0]*TMath::Exp(-0.5*TMath::Power(TMath::Abs(arg),par[3]));
    
    return funcVal;
}

Double_t MyExpSum(Double_t *x, Double_t *par)
{
    Double_t exp1 = MyExp(x,&par[0]);
    Double_t exp2 = MyExp(x,&par[4]);
    Double_t backg = par[8];
    
    return exp1+exp2+backg;
}

void GenHist(TString fname_root = "sim.root") {
    
    auto histogram = DefineHist1D("m_{#nu}^{2}","m_{#nu}^{2} [keV^{2}]",1000,-2E2,2E2);
    auto histogram_test = DefineHist1D("m_{#nu}^{2} test","m_{#nu}^{2} [keV^{2}]",500,0,2E2);
    
    TString tname = TString("tree1");
    
    auto f = TFile::Open(fname_root);
    if (!f || f->IsZombie())
    {
        Error("GenHist2","Cannot open file: %s",fname_root.Data());
        return;
    }
    
    TTreeReader myReader(tname, f);
    
    TTreeReaderValue<Event> myEvent(myReader, "event");
    UInt_t Nentries = myReader.GetEntries(kTRUE);
    
    Double_t gammaTime = 0.;
    
    Double_t x_R;
    Double_t y_R;
    Double_t TOFR;
    Double_t m_R;
    Double_t phi_R;
    Double_t w_R;
    
    Double_t q_ion;
    
    Double_t x_e;
    Double_t y_e;
    Double_t TOFE;
    Double_t phi_e;
    Double_t w_e;
    
    Double_t time_window_R_start = 1.5E-5;
    Double_t time_window_R_end = 7.5E-5;
    Double_t time_window_e_start = 1.5E-7;
    Double_t time_window_e_end = 4.E-7;
    
    Double_t twindow_R = 1.E-4;
    Double_t twindow_e = 8.E-7;
    
    TLorentzVector lor_Cs;
    lor_Cs.SetXYZM(0.,0.,0.,m_Cs);
    
    TLorentzVector lor_g;
    TLorentzVector lor_R;
    std::vector<TLorentzVector> lor_e;
    TLorentzVector lor_v;
    
    UInt_t Nelectrons;
    
    Double_t nx;
    Double_t ny;
    Double_t nz;
    
    Double_t nmag;
    Double_t m_v_sq;
    
    std::vector<Double_t> m_v_sq_vec;
    std::vector<Double_t> energies = {33.562,33.624,34.419};
    
    Bool_t fullEvent;
    
    for (UInt_t i = 0; i < Nentries; i++) {
        myReader.SetEntry(i);
        if (i % 5000 == 0) std::cout << "Entry [" << i << "] of " << Nentries << std::endl;
        // loop through the gammas
        if (myEvent->type == 0) {
            
            if (myEvent->t < gammaTime) continue;
            gammaTime = myEvent->t;
            
            // skip L-shell gammas
            if (myEvent->gamma.E() < 30.) continue;
            
            lor_e.clear();
            Nelectrons = 0;
            
            // get the direction of lor_g, don't precisely know energy
            
            nmag = sqrt(myEvent->gamma.Px() * myEvent->gamma.Px() + myEvent->gamma.Py() * myEvent->gamma.Py() + myEvent->gamma.Pz() * myEvent->gamma.Pz());
            
            nx = myEvent->gamma.Px() / nmag;
            ny = myEvent->gamma.Py() / nmag;
            nz = myEvent->gamma.Pz() / nmag;
            
            // get electron(s) that fall in time window
            
            fullEvent = kFALSE;
            
            for (UInt_t j = i; j < Nentries; j++) {
                
                myReader.SetEntry(j);
                
                if (myEvent->type != 1) continue;
                if (myEvent->t - gammaTime > time_window_e_end) break;
                if (myEvent->t - gammaTime < time_window_e_start) continue;
                
                x_e = myEvent->x;
                y_e = myEvent->y;
                TOFE = myEvent->t - gammaTime;
                w_e = myEvent->w;
                phi_e = myEvent->phi;
                
                // create a single electron TLorentzVector, to add to vector
                TLorentzVector lor_e_single;
                {
                    const Double_t m_e_kg = m_e * keV_to_kg; // electron mass [kg]
                    const Double_t pmex = -(-q) * Bz * sqrt(x_e * x_e + y_e * y_e) * cos(phi_e) / (2. * abs(sin(w_e * TOFE / 2.))) * kgms_to_keV;
                    const Double_t pmey = -(-q) * Bz * sqrt(x_e * x_e + y_e * y_e) * sin(phi_e) / (2. * abs(sin(w_e * TOFE / 2.))) * kgms_to_keV;
                    const Double_t pmez = -(m_e_kg * (-mcp_e)/TOFE + (-q) * Ez * TOFE / 2.) * kgms_to_keV;
                    lor_e_single.SetXYZM(pmex,pmey,pmez,m_e);
                    lor_e.push_back(lor_e_single);
                    Nelectrons++;
                }
                
                fullEvent = kTRUE;
            }
            
            // need to check fullEvent here...
            if (!fullEvent) continue;
            
            // get an ion
            fullEvent = kFALSE;
            
            for (UInt_t j = i; j < Nentries; j++) {
                
                myReader.SetEntry(j);
                
                if (myEvent->type != 2) continue;
                if (myEvent->t - gammaTime > time_window_R_end) break;
                if (myEvent->t - gammaTime < time_window_R_start) continue;
                
                x_R = myEvent->x;
                y_R = myEvent->y;
                TOFR = myEvent->t - gammaTime;
                m_R = myEvent->m;
                w_R = myEvent->w;
                phi_R = myEvent->phi;
                
                {
                    const Double_t m_R_kg = m_R * keV_to_kg; // ion mass[kg]
                    const Double_t q_ion = q * (Double_t)(Nelectrons);
                    
                    const Double_t pmrx = q_ion * Bz * sqrt(x_R * x_R + y_R * y_R) * cos(phi_R) / (2. * abs(sin(w_R * TOFR / 2.))) * kgms_to_keV;
                    const Double_t pmry = q_ion * Bz * sqrt(x_R * x_R + y_R * y_R) * sin(phi_R) / (2. * abs(sin(w_R * TOFR / 2.))) * kgms_to_keV;
                    const Double_t pmrz = (m_R_kg * mcp_R/TOFR - q_ion * Ez * TOFR / 2.) * kgms_to_keV;
                    lor_R.SetXYZM(pmrx,pmry,pmrz,m_R);
                }
                
                fullEvent = kTRUE;
                break;
            }
            
            // ... and here
            if (!fullEvent) continue;
            
            for (UInt_t j = 0; j < energies.size(); j++) {
                
                Double_t gamma_energy = energies.at(j);
                lor_g.SetPxPyPzE(nx * gamma_energy, ny * gamma_energy, nz * gamma_energy, gamma_energy);
                lor_v = lor_Cs-lor_R-lor_g;
                for (UInt_t k = 0; k < Nelectrons; k++) {
                    
                    lor_v = lor_v - lor_e.at(k);
                }
                
                m_v_sq = lor_v.M2();
                histogram->Fill(m_v_sq);
            }
        }
    }
    histogram->Draw();
    
    // fill the test histogram, an alternative to fitting (subtract left side from right, assuming rough symmetry in errors)
    for (Int_t i = 0; i < 500; i++) {
        
        histogram_test->SetBinContent(i, histogram->GetBinContent(500+i) - histogram->GetBinContent(500-i));
        if (histogram_test->GetBinContent(i) < 0) histogram_test->SetBinContent(i, 0);
    }
    
    // no histrogram title because the name is already in the stat box
    histogram->SetBit(TH1::kNoTitle);
    
    // Get guesses for the function parameters by fitting limited ranges
    
    Double_t m_v_sterile_2 = m_v_sterile * m_v_sterile;
    
    // First the m_v=0 peak. It is ok to fix its position at m_v^2=0
    
    TF1 *g1 = new TF1("g1",MyExp,-50,50,4);
    const Int_t binmax = histogram->GetMaximumBin();
    g1->SetParameters(histogram->GetBinContent(binmax),0.,5.,2.);
    g1->FixParameter(1,0.);
    TFitResultPtr g1_fitres = histogram->Fit(g1,"RS");
    
    // Use the exponent and width as was estimated for the m_v^2=0 peak,
    // because no reason to beleiev that detector resolution is m_v dependend
    
    TF1 *g2 = new TF1("g2",MyExp,m_v_sterile_2-10,m_v_sterile_2+10,4);
    g2->SetParameters(m_v_sterile,5,g1->GetParameter(2),g1->GetParameter(3));
    g2->FixParameter(2,g1->GetParameter(2));
    g2->FixParameter(3,g1->GetParameter(3));
    histogram->Fit(g2,"R");
    
    TF1 *l1 = new TF1("l1","pol0",-150,-100);
    histogram->Fit(l1,"QR");
    
    // Define a function that is the sum of Gauss + background and set
    // the initial parameter with the individual fit results
    
    Double_t par[9];
    TF1 *total = new TF1("total",MyExpSum,-100,m_v_sterile_2+100,9);
    total->SetParName(0,"h_0");
    total->SetParName(1,"m_0");
    total->SetParName(2,"s_0");
    total->SetParName(3,"e_0");
    total->SetParName(4,"h_1");
    total->SetParName(5,"m_1");
    total->SetParName(6,"s_1");
    total->SetParName(7,"e_1");
    total->SetParName(8,"backg const");
    total->SetLineColor(kRed);
    
    g1->GetParameters(&par[0]);
    g2->GetParameters(&par[4]);
    l1->GetParameters(&par[8]);
    total->SetParameters(par);
    
    // Also do not allow negative background
    //total->SetParLimits(8,0.,1000000.);
    total->FixParameter(8,0.);
    //Do not allow negative sigma
    total->SetParLimits(2,0.,100.);
    total->SetParLimits(6,0.,100.);
    
    // Show in the histogram the fit statistics
    // four digits: mode = pcev (default = 0111)
    //p = 1 print probability
    //c = 1 print Chi-square/number of degrees of freedom
    //e = 1 print errors (if e=1, v must be 1)
    //v = 1 print name/values of parameters
    
    TFitResultPtr total_fitres = histogram->Fit(total,"RS");
    
    // Estimate the number of events in m_v^2=100 peak by integrating
    // the m_v^2=0 function and the total function across m_1 +/- 2*s_1
    const Double_t xmin = total->GetParameter(5)-2*TMath::Abs(total->GetParameter(6));
    const Double_t xmax = total->GetParameter(5)+2*TMath::Abs(total->GetParameter(6));
    
    const Double_t cnts_backgr = g1->Integral(xmin,xmax);
    const Double_t err_cnts_backgr = g1->IntegralError(xmin,xmax,g1_fitres->GetParams(),
                                                       g1_fitres->GetCovarianceMatrix().GetMatrixArray());
    const Double_t cnts_total      = total->Integral(xmin,xmax);
    const Double_t err_cnts_total  = total->IntegralError(xmin,xmax,total_fitres->GetParams(),
                                                        total_fitres->GetCovarianceMatrix().GetMatrixArray());
    const Double_t cnts_v     = cnts_total-cnts_backgr;
    const Double_t err_cnts_v = TMath::Sqrt(err_cnts_backgr*err_cnts_backgr+err_cnts_total*err_cnts_total);
    
    std::cout << "counts_total : = " << cnts_total  << " +/- " << err_cnts_total  << std::endl;
    std::cout << "counts_backgr: = " << cnts_backgr << " +/- " << err_cnts_backgr << std::endl;
    std::cout << "counts_v     : = " << cnts_v      << " +/- " << err_cnts_v << std::endl;
    
    gStyle->SetOptStat(11);
    gStyle->SetOptFit();
    gStyle->SetStatW(0.4);
    gStyle->SetStatH(0.2);
    
    //    histogram_test->Draw();
}

void Read(TString fname = "sim.root") {
    
    // Set drawing style
    InitStyle();
    
    //Define electron canvas
    auto cv_e = new TCanvas("kinematics electron","kinematics electron",10,10,700,700);
    cv_e->Divide(2,2,0.01,0.03);
    
    cv_e->cd(1);
    auto hist_P_e     = DefineHist1D("p_{e}","p [keV]",50,0.0,50.0);
    cv_e->cd(2);
    auto hist_beta_e  = DefineHist1D("#beta_{e}","#beta",50,0.0,1.0);
    cv_e->cd(3);
    auto hist_theta_e = DefineHist1D("#theta_{e}","#theta [deg]",50,0.0,180.0);
    cv_e->cd(4);
    auto hist_phi_e   = DefineHist1D("#phi_{e}","#phi [deg]",50,-180.0,180.0);
    
    //Define Nucleus canvas
    auto cv_R = new TCanvas("kinematics nucleon","kinematics nucleon",60,60,700,700);
    cv_R->Divide(2,2,0.01,0.03);
    
    cv_R->cd(1);
    auto hist_P_R     = DefineHist1D("p_{R}","p [keV]",50,0.0,1000.0);
    cv_R->cd(2);
    auto hist_beta_R  = DefineHist1D("#beta_{R}","#beta",50,0.0,1.0);
    cv_R->cd(3);
    auto hist_theta_R = DefineHist1D("#theta_{R}","#theta [deg]",50,0.0,180.0);
    cv_R->cd(4);
    auto hist_phi_R   = DefineHist1D("#phi_{R}","#phi [deg]",50,-180.0,180.0);
    
    //Define neutrino canvas
    auto cv_v = new TCanvas("kinematics neutrino","kinematics neutrino",110,110,700,700);
    cv_v->Divide(2,2,0.01,0.03);
    
    cv_v->cd(1);
    auto hist_P_v     = DefineHist1D("p_{#nu}","p [keV]",50,0.0,1000.0);
    cv_v->cd(2);
    auto hist_M2_v    = DefineHist1D("M^{2}_{#nu}","M^{2} [keV^2]",50,-100.0,100.0);
    cv_v->cd(3);
    auto hist_theta_v = DefineHist1D("#theta_{#nu}","#theta [deg]",50,0.0,180.0);
    cv_v->cd(4);
    auto hist_phi_v   = DefineHist1D("#phi_{#nu}","#phi [deg]",50,-180.0,180.0);
    
    // Open the file containing the tree.
    TFile *myFile = TFile::Open(fname);
    if (!myFile || myFile->IsZombie()) {
        return;
    }
    // Create a TTreeReader for the tree, for instance by passing the
    // TTree's name and the TDirectory / TFile it is in.
    TTreeReader myReader("tree2",myFile);
    
    TTreeReaderValue<TLorentzVector> myLor_g(myReader,"lor_g");
    TTreeReaderValue<TLorentzVector> myLor_R_m_c(myReader,"lor_R_m_c");
    TTreeReaderValue<TLorentzVector> myLor_e_m_c(myReader,"lor_e_m_c");
    TTreeReaderValue<TLorentzVector> myLor_v_m_c(myReader,"lor_v_m_c");
    
    // Loop over all entries of the TTree or TChain.
    while (myReader.Next()) {
        
        hist_P_e->Fill(myLor_e_m_c->P());
        hist_beta_e->Fill(myLor_e_m_c->P()/myLor_e_m_c->M());
        hist_theta_e->Fill(myLor_e_m_c->Theta()*180/M_PI);
        hist_phi_e->Fill(myLor_e_m_c->Phi()*180/M_PI);
        
        hist_P_R->Fill(myLor_R_m_c->P());
        hist_beta_R->Fill(myLor_R_m_c->P()/myLor_R_m_c->M());
        hist_theta_R->Fill(myLor_R_m_c->Theta()*180/M_PI);
        hist_phi_R->Fill(myLor_R_m_c->Phi()*180/M_PI);
        
        hist_P_v->Fill(myLor_v_m_c->E());
        hist_M2_v->Fill(myLor_v_m_c->M2());
        hist_theta_v->Fill(myLor_v_m_c->Theta()*180/M_PI);
        hist_phi_v->Fill(myLor_v_m_c->Phi()*180/M_PI);
    }
    
    cv_e->Modified();
    cv_e->Update();
    cv_e->SaveAs("electron.png");
    
    cv_R->Modified();
    cv_R->Update();
    cv_R->SaveAs("nucleon.png");
    
    cv_v->Modified();
    cv_v->Update();
    cv_v->SaveAs("neutrino.png");
}
