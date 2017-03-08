/* Projet de 3A
 * Temigo
 */
#include "event.h"
#include "globals.h"

#include <stdio.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TEllipse.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TApplication.h>
#include <TFile.h>
#include <TVector.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TArrow.h>
#include <TSystem.h>
#include <TROOT.h>

#include <queue>
#include <iostream>
#include <fstream>

IntegralFunction::IntegralFunction(TF2 * f, const char * name):
    function(f),
    name(name) 
    {
        f12 = new TF12(name, function, 1.0, "y");
    }

double IntegralFunction::operator() (double * x, double * p) const
{
    if (*x < 0)
    {
        return 0.0;
    }
    else
    {
        function->SetParameter(0, p[0]);
        f12->SetXY(*x);
        //std::cout << function->GetNumberFreeParameters() << std::endl;
        return f12->Integral(phi(*x), TMath::Pi());
    }
}

Event::Event(Double_t rho, Double_t max_y, Double_t R, const char * lut_filename, TF1 * f, bool with_cutoff, bool raw_cutoff, bool minimal) :
    rho(rho),
    max_y(max_y),
    R(R),
    lut_filename(lut_filename),
    cutoff(f),
    WITH_CUTOFF(with_cutoff),
    RAW_CUTOFF(raw_cutoff),
    MINIMAL(minimal)
{
    if (WITH_CUTOFF && cutoff == NULL)
    {
        std::cout << "[event.cpp] Error initializing event : cutoff undefined." << std::endl;
    }
    if (WITH_CUTOFF && cutoff != NULL)
    {
        cutoff->SetParameter(1, R); // FIXME set or fix parameter ?

        integrand = new TF2("integrand", "cutoff / (x * (1. + x^2 - 2. * x * cos(y)))", 0.0, TMath::Infinity(), 0, TMath::Pi());   
        integralFunction = new IntegralFunction(integrand, "f_cutoff");
        f_cutoff = new TF1("f_cutoff", integralFunction, 0.0 , TMath::Infinity(), 1);
        f_cutoff->SetNpx(50); // FIXME Npx value
    }    
}

Event::~Event()
{
}

Double_t Event::f(Double_t * x, Double_t * parameters)
{
    Double_t r = *x;
    if (r <= 1. / 2. )
    {
        return TMath::Pi() / (r * (1. -r*r));
    }
    else
    {
        return 2. / (r * TMath::Abs(1. -r*r)) * TMath::ATan(TMath::Abs(1. -r) / (1. +r) * TMath::Sqrt((r+1. / 2. )/(r-1. / 2.)));
    }
}

Double_t Event::g(Double_t r)
{
    return 5. * TMath::Pi() / (3. * r * (1. + r * r));
}

// Boundary for theta approaching zero
Double_t phi(Double_t r)
{
    if (r > 1.0 / 2.0)
    {
        return TMath::ACos(1.0 / (2.0 * r));
    }
    return 0.0;
}

// Rejection method to generate radius r
Double_t Event::r_generate(Double_t x01, bool display_cutoff)
{
    if (WITH_CUTOFF && RAW_CUTOFF)
    {
        gRandom->SetSeed();
        f_cutoff->SetParameter(0, x01);
        //gROOT->GetListOfFunctions()->Print();

        // Cut-off on small sizes : we don't need dipoles of size < rho
        TF1 g2("g2", "5 * TMath::Pi() / 3 * 1/(x * (1+x^2)) * exp(- [0] / (2*[1]^2) * (2*x^2))", rho, 100);
        g2.SetParameter(0, 1.0);
        g2.SetParameter(1, R);

        TF1 g3("g3", "5 * TMath::Pi() / 3 * 1/(x * (1+x^2))", 0, 100);

        // Draw f(r)
        if (display_cutoff)
        {
            TCanvas C("C", "C", 0, 0, 2000, 1000);
            C.cd(1);
            gPad->SetLogy();
            f_cutoff->DrawF1(0.0, 10);

            TF1 f_old("f_old", f, 0.0, 100, 0);
            f_old.SetLineColor(42);
            f_old.DrawF1(0.0, 10, "same");

            TF1 g("g", "5 * TMath::Pi() / 3 * 1/(x * (1+x^2)) * exp(- [0] / (2*[1]^2) * (1 + 2*x^2-2*x))", 0, 100);
            g.SetParameter(0, 1.0);
            g.SetParameter(1, R);
            g.SetLineColor(46);
            g.DrawF1(0.0, 10, "same");

            // Closest to the gaussian cutoff
            g2.SetLineColor(38);
            g2.DrawF1(0.1, 10, "same");

            g3.SetLineColor(30);
            g3.DrawF1(0.0, 10, "same");

            gPad->Update();
            char c;
            std::cin>>c;
        }

        // Get random following f_cutoff distribution using rejection sampling
        Double_t rho2 = rho / x01;
        Double_t R = gRandom->Uniform(0., 1.); // FIXME conflict with R ?
        Double_t result = 1. / TMath::Sqrt(TMath::Exp((1. - R) * TMath::Log(1. + 1. /(rho2 * rho2)))- 1. );
        Double_t temp = gRandom->Uniform(0., 1.);
        while (temp > f_cutoff->Eval(result) / g(result))
        {
            R = gRandom->Uniform(0., 1.);
            result = 1. / TMath::Sqrt(TMath::Exp((1. - R) * TMath::Log(1. + 1. /(rho2 * rho2)))- 1.);
            temp = gRandom->Uniform(0., 1.);
        }
        return result;
    }
    else
    {
        Double_t rho2 = rho / x01;
        Double_t R = gRandom->Uniform(0., 1.);
        Double_t result = 1. / TMath::Sqrt(TMath::Exp((1. - R) * TMath::Log(1. + 1. /(rho2 * rho2)))- 1. );

        Double_t temp = gRandom->Uniform(0., 1.);
        while (temp > f(&result) / g(result) )//* 3. / (5. * TMath::Pi()))
        {
            R = gRandom->Uniform(0., 1.);
            result = 1. / TMath::Sqrt(TMath::Exp((1. - R) * TMath::Log(1. + 1. /(rho2 * rho2)))- 1.);
            temp = gRandom->Uniform(0., 1.);
        }
        return result;
    }
}

// Distribution for theta
Double_t Event::theta(Double_t r, Double_t x01)
{
    gRandom->SetSeed();
    if (WITH_CUTOFF && RAW_CUTOFF)
    {
        // FIXME r range for g_theta : begin at rho ?
        // FIXME is proportional factor important here ?
        //TF2 g_theta("g_theta", "1/[2] * 1/(x * (1. + x^2 - 2. * x * cos(y))) * cutoff", 0.0, TMath::Infinity(), 0.0, TMath::Pi());
        TF2 g_theta("g_theta", " 1/(1. + x^2 - 2. * x * cos(y)) * cutoff", 0.0, TMath::Infinity(), phi(r), TMath::Pi());
        g_theta.SetParameter(0, x01);
        TF12 g_theta12("g_theta12", &g_theta, r, "y");
        return g_theta12.GetRandom(phi(r), TMath::Pi());
    }
    else
    {
        Double_t a = TMath::Abs(1. - r)/(1. + r);
        Double_t R = gRandom->Uniform(0., 1.);
        if (r <= 1. / 2.) {
            return 2. * TMath::ATan((1. - r)/(1. + r) * 1. / TMath::Tan(R * TMath::Pi() / 2. ));
        }
        else {
            return 2. * TMath::ATan(a * 1. / TMath::Tan(R * TMath::ATan(a * TMath::Sqrt((r+1. / 2. )/(r-1. / 2. )))));
        }
    }
}

// For Rapidity distribution - lifetime
// Raw computation of lambda - for lookup tables see getLambda
Double_t Event::lambda(Double_t x01)
{
    if (WITH_CUTOFF)
    {
        f_cutoff->SetParameter(0, x01);
        //TF1 * integral_function = new TF1("integral_function", "f_cutoff * 2/TMath::Pi()", 0.0, TMath::Infinity());
        //f_cutoff->DrawF1(0.0, 2.);
        return 2.0 / TMath::Pi() * f_cutoff->Integral(rho / x01, TMath::Infinity());
    }
    else
    {
        Double_t resultat = 0. ;
        Double_t rho2 = rho / x01; // Scale first

        // Compute for rho <= 1/2
        if (rho2 < 1. / 2. ) {
            resultat = TMath::Log(1. / 3. * (1. / (rho2 * rho2) - 1. ));
        }
        // Compute for rho > 1/2
        Double_t borne = (rho2 < 1. / 2.) ? 1. / 2. : rho2;
        TF1 integral_function("integral_function", "1. /(x * TMath::Abs(1. -x*x)) * TMath::ATan(TMath::Abs(1. -x)/(1. +x) * TMath::Sqrt((x+1. /2. )/(x-1. /2. )))", borne, TMath::Infinity());
        //integral_function->DrawF1(1., 100.);
        resultat = resultat + 4. / TMath::Pi() * integral_function.Integral(borne, TMath::Infinity());

        return resultat;
    }
}

/* Creates and saves in __filename__ a lookup table
 * Logarithmic scale
 */
void Event::WriteLookupTable()
{
    int n = 320;

    Double_t x01 = 0.0;
    // We shouldn't need to go under the cutoff rho for small sizes
    Double_t step = 0.0000000000000000000000000000001;
    Double_t l;
    if (DEBUG) std::cout << "Writing lookup table... ";
    std::ofstream lut(lut_filename);
    if (lut.is_open())
    {
        lut << rho << " " << WITH_CUTOFF << "\n";
        for (int i = 1; i <= n; ++i)
        {
            if (i%10 == 0) step *= 10. ;
            x01 += step;
            l = lambda(x01);
            lut << x01 << " " << l << "\n";
        }
        lut.close();
    }
    if (DEBUG) std::cout << "Done." << std::endl;
}

/* Loads the lookup table from __filename__ and put it into Interpolator
 */
void Event::LoadLookupTable()
{
    if (DEBUG) std::cout << "Reading lookup table..." << std::endl;
    std::ifstream lut(lut_filename);
    Double_t rho_lut;
    bool cutoff;
    std::vector<Double_t> x01, l;
    if (lut.is_open())
    {
        lut >> rho_lut >> cutoff;
        if (rho_lut != rho || cutoff != WITH_CUTOFF)
        {
            if (rho_lut != rho) std::cout << RED << "Warning : rho is not the same in the Lookup Table. Rebuilding it..." << RESET << " \n Old rho : " << rho_lut << std::endl;
            if (cutoff != WITH_CUTOFF) std::cout << RED << "Warning : WITH_CUTOFF differs in the Lookup Table. Rebuilding it..." << RESET << " \n Old WITH_CUTOFF : " << cutoff << std::endl;
            lut.close();
            WriteLookupTable();
            lut.clear();
            lut.open(lut_filename);
            lut >> rho_lut >> cutoff;
        }
        Double_t a, b;
        while(lut >> a >> b)
        {
            lookup_table.insert(std::pair<Double_t, Double_t>(a, b));
            x01.push_back(a);
            l.push_back(b);
        }
        lut.close();
    }
    interpolator.SetData(x01, l);
    x01_min = x01.front();
    x01_max = x01.back();
    if (DEBUG) std::cout << "\033[1;32m Done. \033[0m" << std::endl;
}

void Event::PrintLookupTable()
{
    if (DEBUG) std::cout << "Printing lookup table..." << std::endl;
    for (auto const& x: lookup_table)
    {
        std::cout << x.first << " " << x.second << std::endl;
    }
    if (DEBUG) std::cout << "done." << std::endl;    
}

void Event::SetInterpolatorData()
{
    std::vector<Double_t> x01, l;
    for (auto const& x: lookup_table)
    {
        x01.push_back(x.first);
        l.push_back(x.second);
    }
    interpolator.SetData(x01, l);
}

Double_t Event::getLambda(Double_t x01)
{
    if (lookup_table.size() <= 0)
    {
        LoadLookupTable();
    }
    // Interpolate
    if (x01_min > x01)
    {
        if (DEBUG) std::cout << RED << "Warning : x01 out of lookup table range (lower values). This is not supposed to happen !" << RESET << std::endl;
    }
    if (x01 > x01_max)
    {
        if (DEBUG) std::cout << RED << "Warning : x01 out of lookup table range (upper values)." << RESET << " \n Adding entry for " << x01 << std::endl;
        // Add new points permanently in the LUT
        Double_t x01_from = x01_max;
        int p = TMath::Floor(TMath::Log10(x01_from));
        Double_t step = TMath::Power(10, p);
        Double_t l;
        std::ofstream lut(lut_filename, std::ofstream::out | std::ofstream::app);
        if (lut.is_open())
        {
            int i = (int) TMath::Floor(TMath::Exp(TMath::Log(x01_from)-p * TMath::Log(10))) + 1;
            while (x01_from < x01)
            {
                if (i%10 == 0) step *= 10. ;
                x01_from += step;
                l = lambda(x01_from);
                lookup_table.insert(std::pair<Double_t, Double_t>(x01_from, l));
                lut << x01_from << " " << l << "\n";
                ++i;
                if (DEBUG) std::cout << x01_from << " " << l << std::endl;
            }
            lut.close();
        }

        x01_max = x01_from;
        SetInterpolatorData();
    }
    return interpolator.Eval(x01);
}

// Generate rapidity
Double_t Event::y_generate(Double_t x01)
{
    gRandom->SetSeed();
    Double_t lbda = getLambda(x01);
    Double_t R = gRandom->Uniform(0., 1. / lbda);
    return -1. / lbda * TMath::Log(1. - lbda * R);
}

// To be used in tests - incorporates the cutoff directly
Double_t Event::r_generate2(Double_t x01)
{
    gRandom->SetSeed();
    Double_t r = r_generate(x01);
    Double_t t = theta(r, x01);
    Double_t temp = gRandom->Uniform(0.0, 1.0);
    cutoff->SetParameter(0, x01);
    while(temp > cutoff->Eval(r, t))
    {
        r = r_generate(x01);
        t = theta(r, x01);
        temp = gRandom->Uniform(0.0, 1.0);
    } 
    return r;   
}

// Draw a single gluon with an ellipse
void Event::draw(Double_t x, Double_t y, Double_t rapidity)
{
    TEllipse ellipse(x, y, 0.005 , 0.005);
    ellipse.SetFillColor(kCyan + (int) TMath::Ceil(10. * rapidity));
    ellipse.Draw();
}

bool Event::generate(Dipole * dipole, Dipole * dipole1, Dipole * dipole2, Double_t max_y)
{
    // First generate rho and theta in normalized referential
    gRandom->SetSeed();
    Double_t r = r_generate(dipole->radius);
    Double_t t = theta(r, dipole->radius);

    if(WITH_CUTOFF && !RAW_CUTOFF) // Reject with probability given by cutoff
    {
        Double_t temp = gRandom->Uniform(0.0, 1.0);
        cutoff->SetParameter(0, dipole->radius);
        while(temp > cutoff->Eval(r, t))
        {
            r = r_generate(dipole->radius);
            t = theta(r, dipole->radius);
            temp = gRandom->Uniform(0.0, 1.0);
        }
    }

    Double_t rapidity = y_generate(dipole->radius);   
    Double_t hb = gRandom->Uniform(0., 1.); // Up or down quadrant
    Double_t gd = gRandom->Uniform(0., 1.); // Left or right quadrant

    if (dipole->rapidity + rapidity > max_y) return false;

    dipole1->rapidity = dipole->rapidity + rapidity;
    dipole2->rapidity = dipole->rapidity + rapidity;

    // Set coordinates of centers in normalized referential
    // Use cartesian coordinates instead of SetMagPhi because of the four quadrants
    if (hb < 0.5)
    {
        dipole1->coord.Set(r * TMath::Cos(t), r * TMath::Sin(t));
        dipole1->phi = t;
        dipole1->radius = r / 2. ;
    }
    else
    {
        dipole1->coord.Set(r * TMath::Cos(t), - r * TMath::Sin(t));
        dipole1->phi = -t;
        dipole1->radius = r / 2. ;
    }

    //dipole1->coord.Print();
    // Get the right coordinates for the parameter impact
    dipole1->coord /= 2.;
    dipole2->coord = dipole1->coord + TVector2(0.5, 0.0);

    TVector2 temp = dipole2->coord - TVector2(1, 0);
    dipole2->phi = dipole2->coord.Phi_0_2pi(temp.Phi());
    dipole2->radius = temp.Mod();
    
    // Exchange
    if (gd < 0.5)
    {
        Double_t x = dipole2->coord.X();
        dipole2->coord.Set(- dipole1->coord.X(), dipole1->coord.Y());
        dipole2->coord += TVector2(1.0, 0.0);
        dipole1->coord.Set(1.0 - x, dipole2->coord.Y());
        Double_t phi = dipole1->phi;
        dipole1->phi = TMath::Pi() - dipole2->phi;
        dipole2->phi = - phi;
        Double_t radius = dipole1->radius;
        dipole1->radius = dipole2->radius;
        dipole2->radius = radius;    
    }

    // Now we have dipole1 and dipole2 in the normalized referential, get back
    // Scale
    Double_t scale_factor = dipole->radius;
    dipole1->coord *= scale_factor;
    dipole2->coord *= scale_factor;
    dipole1->radius *= scale_factor;
    dipole2->radius *= scale_factor;

    // Rotate and stay in 0 - 2pi
    dipole1->coord = TVector2(dipole1->coord.X()-dipole->radius/2., dipole1->coord.Y()).Rotate(dipole->phi) + TVector2(dipole->radius/2., 0);
    dipole2->coord = TVector2(dipole2->coord.X()-dipole->radius/2., dipole2->coord.Y()).Rotate(dipole->phi) + TVector2(dipole->radius/2., 0);

    //dipole1->phi = dipole1->coord.Phi_0_2pi(dipole1->phi + dipole->phi);
    //dipole2->phi = dipole2->coord.Phi_0_2pi(dipole2->phi + dipole->phi);
    dipole1->phi = dipole1->phi + dipole->phi;
    dipole2->phi = dipole2->phi + dipole->phi;

    // Translation
    TVector2 translation_vector = dipole->coord - TVector2(dipole->radius/2., 0.0);
    dipole1->coord += translation_vector;
    dipole2->coord += translation_vector;

    // Back from radius to size of dipoles
    dipole1->radius *= 2;
    dipole2->radius *= 2;

    if (DEBUG && (dipole1->radius < rho || dipole2->radius < rho)) std::cout << dipole1->radius << " " << dipole2->radius << std::endl;
    return true;
}

void Event::make_tree(TTree * tree, bool draw_dipole, bool draw_step_by_step)
{
    gRandom->SetSeed(); // /!\ IMPORTANT or we always get the same values
    Dipole dipole(0,0); // must be initialized before setting branch adresses

    std::queue<Dipole> dipoles;
    dipoles.push(dipole);

    //TTree * tree = new TTree(treename, "Dipole splitting");
    // Essential branches
    tree->Branch("rapidity", &dipole.rapidity, "rapidity/D");
    tree->Branch("radius", &dipole.radius, "radius/D");
    //tree->Branch("isLeaf", &dipole.isLeaf, "isLeaf/O");
    // For drawing purpose ?
    //tree->Branch("phi", &dipole.phi, "phi/D");
    //tree->Branch("coord", "TVector2", &dipole.coord);
    // For common ancestor
    //tree->Branch("depth", &dipole.depth, "depth/L");
    //tree->Branch("index_children", &dipole.index_children, "index_children/L");
    //tree->Branch("index_parent", &dipole.index_parent, "index_parent/L");

    int i = 0;
    Long64_t index = -1; // index of the dipole at current depth

    if (draw_dipole || draw_step_by_step)
    {
        gPad->DrawFrame(-1.0, -1, 2.0, 1.0, TString::Format("#splitline{Dipole splitting - rho = %.12g}{rapidity maximum = %.12g}", rho, max_y));        
    }

    Long64_t current_depth = 0;
    Long64_t current_depth_dipoles = 1; // Number of dipoles at current depth
    Long64_t current_depth_dipoles_split = 0;

    // Tree is filled in a Breadth first search (BFS) order
    while (!dipoles.empty())
    {
        dipole = dipoles.front();
        dipoles.pop();

        // Determine current depth
        if (dipole.depth > current_depth) // Beginning to fill a new level
        {
            current_depth = dipole.depth; 
            // Constant throughout the current depth
            current_depth_dipoles = 2 * current_depth_dipoles_split;
            // Incremented throughout the current depth
            current_depth_dipoles_split = 0;
            index = 0;
        }
        else // Increment index of the dipole in current depth
        {
            if (index < std::numeric_limits<Long64_t>::max())
            {
                index += 1;
            }
            else
            {
                std::cout << "Warning : overflow on index value" << std::endl;
                break;
            }
        }
        dipole.index = i;
        dipole.nb_left_brothers_split = current_depth_dipoles_split;
        dipole.nb_right_brothers = current_depth_dipoles - index - 1;

        Long64_t children_index = dipole.index + dipole.nb_right_brothers + 2 * dipole.nb_left_brothers_split + 1;
        //std::cout << dipole.index << " " << dipole.nb_left_brothers_split << " " << dipole.nb_right_brothers << " " << children_index << std::endl;
        
        Dipole dipole1(dipole.depth + 1, children_index), 
               dipole2(dipole.depth + 1, children_index + 1);
        bool success = generate(&dipole, &dipole1, &dipole2, max_y);
        if (success)
        {
            dipole1.index_parent = dipole.index;
            dipole2.index_parent = dipole.index;
            dipoles.push(dipole1);
            dipoles.push(dipole2);
            current_depth_dipoles_split += 1;    
            dipole.index_children = children_index;
        }
        else
        {
            dipole.isLeaf = true;
        }
        // remove the dipole.isLeaf to draw in blue the old dipoles
        if ((draw_dipole) || draw_step_by_step)
        {
            dipole.Draw();
        }
                  
        // FIXME for some reason the raw values differ slightly from the tree scan ! 
        //std::cout << dipole1.radius << " " << dipole1.coord.X() << " " << dipole1.coord.Y() <<  std::endl;
        ++i;
        if (dipole.isLeaf || !MINIMAL)
        {
            tree->Fill(); // Stores dipole (parent)
        }
        // if step by step
        if (draw_step_by_step)
        {
            gPad->Modified();
            gPad->Update();
            gPad->WaitPrimitive();
            std::cin.ignore(); // Wait for keypress before continuing
        }
        if (DEBUG) std::cout << "\r" << i << " dipôles générés";
    }
    if (draw_dipole || draw_step_by_step) 
    {
        gPad->Update();
        char c;
        std::cin>>c;
    }

    if (DEBUG) std::cout << std::endl;

    /*TVector * v = new TVector(2);
    //v->SetName("simulation_parameters");
    v[0] = rho;
    v[1] = max_y;
    tree->GetUserInfo()->Add(v);*/

    //tree->Print();
    // Print entries
    //tree->Scan("rapidity:radius:phi:coord.X():coord.Y():depth:index_children:index_parent");
    //tree->Draw("radius");
    //output->Write();
    tree->ResetBranchAddresses(); // all branches must set their adresses to zero
    //return tree;
}

std::vector<int> Event::make_histograms(std::vector<Double_t> r)
{
    gRandom->SetSeed(); // /!\ IMPORTANT or we always get the same values
    Dipole dipole(0,0); // must be initialized before setting branch adresses

    std::vector<int> n; // number of dipoles of size >= r
    for (int k = 0; k < r.size(); ++k)
    {
        n.push_back(0);
    }

    std::queue<Dipole> dipoles;
    dipoles.push(dipole);

    int i = 0;
    Long64_t index = -1; // index of the dipole at current depth

    Long64_t current_depth = 0;
    Long64_t current_depth_dipoles = 1; // Number of dipoles at current depth
    Long64_t current_depth_dipoles_split = 0;

    // Tree is filled in a Breadth first search (BFS) order
    while (!dipoles.empty())
    {
        dipole = dipoles.front();
        dipoles.pop();

        // Determine current depth
        if (dipole.depth > current_depth) // Beginning to fill a new level
        {
            current_depth = dipole.depth; 
            // Constant throughout the current depth
            current_depth_dipoles = 2 * current_depth_dipoles_split;
            // Incremented throughout the current depth
            current_depth_dipoles_split = 0;
            index = 0;
        }
        else // Increment index of the dipole in current depth
        {
            if (index < std::numeric_limits<Long64_t>::max())
            {
                index += 1;
            }
            else
            {
                std::cout << "Warning : overflow on index value" << std::endl;
                break;
            }
        }
        dipole.index = i;
        dipole.nb_left_brothers_split = current_depth_dipoles_split;
        dipole.nb_right_brothers = current_depth_dipoles - index - 1;

        Long64_t children_index = dipole.index + dipole.nb_right_brothers + 2 * dipole.nb_left_brothers_split + 1;
        //std::cout << dipole.index << " " << dipole.nb_left_brothers_split << " " << dipole.nb_right_brothers << " " << children_index << std::endl;
        
        Dipole dipole1(dipole.depth + 1, children_index), 
               dipole2(dipole.depth + 1, children_index + 1);
        bool success = generate(&dipole, &dipole1, &dipole2, max_y);
        if (success)
        {
            dipole1.index_parent = dipole.index;
            dipole2.index_parent = dipole.index;
            dipoles.push(dipole1);
            dipoles.push(dipole2);
            current_depth_dipoles_split += 1;    
            dipole.index_children = children_index;
        }
        else
        {
            dipole.isLeaf = true;
        }
                  
        // FIXME for some reason the raw values differ slightly from the tree scan ! 
        //std::cout << dipole1.radius << " " << dipole1.coord.X() << " " << dipole1.coord.Y() <<  std::endl;
        ++i;
        if (dipole.isLeaf)
        {
            for (int k = 0; k < r.size(); ++k)
            {
                if (dipole.radius >= r[k]) n[k] += 1;
            }
        }
        if (DEBUG) std::cout << "\r" << i << " dipôles générés";
    }

    if (DEBUG) std::cout << std::endl;
    return n;
}

/********************************************* NORMALIZED + FIT *************************/

// Split 1 dipole - to test the distribution
void Event::generate_normalized(Double_t * x, Double_t * y, Double_t * rapidity)
{
    // Generate r, theta and y
    Double_t r = r_generate();
    Double_t t = theta(r);

    if(WITH_CUTOFF && !RAW_CUTOFF) // Reject with probability given by cutoff
    {
        Double_t temp = gRandom->Uniform(0.0, 1.0);
        cutoff->SetParameter(0, 1.0);
        while(temp > cutoff->Eval(r, t))
        {
            r = r_generate();
            t = theta(r);
            temp = gRandom->Uniform(0.0, 1.0);
        }
    }

    *rapidity = y_generate(r);

    // Compute dP
    Double_t quadrant = gRandom->Uniform(0., 1.);
    if (quadrant < 0.25)
    {
        *x = r * TMath::Cos(t);
        *y = r * TMath::Sin(t);
    }
    else if (quadrant < 0.5)
    {
        *x = r * TMath::Cos(t);
        *y = - r * TMath::Sin(t);
    }
    else if (quadrant < 0.75)
    {
        *x = 1.0 - r * TMath::Cos(t);
        *y = r * TMath::Sin(t);            
    }
    else
    {
        *x = 1.0 - r * TMath::Cos(t);
        *y = - r * TMath::Sin(t);            
    }
    //std::cout << *x << " " << *y << std::endl;
}

// Test one generation : visualize probability distribution of the gluon at splitting
void Event::bare_distribution()
{
    bool DRAW_ELLIPSES = false; // If we want different colors, use ellipses
    bool DRAW_STEP_BY_STEP = false;
    int number_occurrences = 3000;

    Double_t x[number_occurrences], y[number_occurrences], rapidity[number_occurrences];

    TCanvas C("C", "C", 0, 0, 1024, 768);
    //if (!DRAW_STEP_BY_STEP) C->Divide(2, 1, 0.05, 0.05);

    C.cd(1);
    gPad->SetTitle("QCD");
    gPad->DrawFrame(-1.0, -0.8, 2.0, 0.8, "Gluons");
    
    TH1D hist("hist", "Plot X", 100, -1, 2);
    Double_t hist_y = 0.04;
    Double_t margin = 0.01;

    //TF1 * r_distribution = new TF1("r_distribution", r_distribution, rho, 1);
    //TF1 * theta_distribution = new TF1("theta_distribution", theta, 0, TMath::Pi());
    // Generate gluons according to dP
    for (int i = 0; i < number_occurrences; ++i)
    {
        generate_normalized(&x[i], &y[i], &rapidity[i]);
        if (DRAW_ELLIPSES || DRAW_STEP_BY_STEP) draw(x[i], y[i], rapidity[i]);
        if (DRAW_STEP_BY_STEP)
        {
            C.Modified();
            C.Update();
            C.WaitPrimitive();
            std::cin.ignore();
        }
        if (y[i] > (hist_y - margin) && y[i] < (hist_y + margin)) hist.Fill(x[i]);
        //std::cout << i << " " << x[i] << " " << y[i] << std::endl;
    }
    C.Update();

    if (!DRAW_ELLIPSES && !DRAW_STEP_BY_STEP)
    {
        TGraph gluons(number_occurrences, x, y);
        gluons.SetTitle("P3A");
        gluons.Draw("AP");
        gluons.SetMarkerColor(2);
        gluons.GetXaxis()->SetRangeUser(-1., 2.);
        gluons.GetYaxis()->SetRangeUser(-0.8, 0.8);
        gluons.Draw("AP");
        C.Update();
        char c;
        std::cin>>c;        
    }

    //C.Update();

    /*C.cd(2);
    gPad->SetLogy();
    hist.Draw("E1");

    TF1 p("p", "[0] / ((x*x + [1]*[1])*((x-1)*(x-1)+[1]*[1]))", -0.5, 1.5);
    p.FixParameter(1, hist_y);
    hist.Fit("p", "R");

    Double_t chi2 = p.GetChisquare();
    std::cout << "Chi2 value : " << chi2 << std::endl;  */ 
}

void Event::fit_r()
{
    int number_occurrences = 10000;

    TCanvas C("C", "Fit r", 0, 0, 1024, 768);
    TH1D hist("hist", "hist", 400, 0, 1);
    gPad->SetLogy();

    for (int i = 0 ; i < number_occurrences ; ++i)
    {
        hist.Fill(r_generate());
    }
    hist.Draw("E1");

    TF1 fr1("fr1", "[0] / (x*(1. - x*x))", rho, 0.5);
    if (!WITH_CUTOFF) hist.Fit("fr1", "R");

    TF1 fr2("fr2", "2. * [0] / (x*TMath::Abs(1. - x*x)) * TMath::ATan(TMath::Abs(1. - x)/(1. + x) * TMath::Sqrt((x+1./2.)/(x-1./2.)))", 0.5, 2);
    if (!WITH_CUTOFF) hist.Fit("fr2", "R+");

    if (!WITH_CUTOFF)
    {
        Double_t chi2_1 = fr1.GetChisquare();
        Double_t chi2_2 = fr2.GetChisquare();
        hist.SetTitle(TString::Format("#splitline{Distribution de la taille r du dipole (#rho = %.12g)}{%d bins - chi2 = %.12g et %.12g}", rho, hist.GetSize()-2, chi2_1, chi2_2));
    }
    C.Update();
    char c;
    std::cin>>c;
}

void Event::fit_y()
{
    int number_occurrences = 1000000;
    Double_t x01 = 1.0;

    TCanvas C("C", "Fit y", 0, 0, 1024, 768);
    TH1D hist("hist", "hist", 100, 0, 1);
    gPad->SetLogy();

    for (int i = 0 ; i < number_occurrences ; ++i)
    {
        hist.Fill(y_generate(x01));
    }
    hist.Draw("E1");

    TF1 fy("fy", "[0] * exp(- [1] * x)", 0, 1);
    fy.FixParameter(1, getLambda(x01));
    if (!WITH_CUTOFF)
    {
        hist.Fit("fy", "R");
        Double_t chi2 = fy.GetChisquare();
        hist.SetTitle(TString::Format("#splitline{Distribution de la rapidite y du dipole (#rho = %.12g)}{%d bins - chi2 = %.12g}", rho, hist.GetSize()-2, chi2));  
    }
    C.Update();
    char c;
    std::cin>>c;
}

void Event::fit_x()
{
    int number_occurrences = 1000000;
    Double_t x, y, rapidity;

    TCanvas C("C", "C", 0, 0, 1024, 768);
    
    Double_t hist_y = 0.2;
    Double_t margin = 0.005;
    TH1D hist("hist", "hist", 100, -0.5, 1.5);

    //TF1 * r_distribution = new TF1("r_distribution", r_distribution, rho, 1);
    //TF1 * theta_distribution = new TF1("theta_distribution", theta, 0, TMath::Pi());
    // Generate gluons according to dP
    for (int i = 0; i < number_occurrences; ++i)
    {
        generate_normalized(&x, &y, &rapidity);
        if (y > (hist_y - margin) && y < (hist_y + margin)) hist.Fill(x);
    }

    gPad->SetLogy();
    hist.Draw("E1");

    TF1 p("p", "[0] / ((x*x + [1]*[1])*((x-1)*(x-1)+[1]*[1]))", -0.3, 1.3);
    p.FixParameter(1, hist_y);
    hist.Fit("p", "R");

    Double_t chi2 = p.GetChisquare();
    hist.SetTitle(TString::Format("#splitline{Distribution of projection X with Y = %.12g +/- %.12g - #rho = %.12g}{%d bins - chi2 = %.12g}", hist_y, margin, rho, hist.GetSize()-2, chi2)); 
}