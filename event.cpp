/* Projet de 3A
 * Temigo
 */
#include "event.h"

#include <stdio.h>
#include <TMath.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TEllipse.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TApplication.h>
#include <TFile.h>
#include <queue>
#include <TVector.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TArrow.h>
#include <iostream>
#include <fstream>
#include <TSystem.h>

Double_t RHO = 0.;
Double_t LAMBDA = 0.;
Double_t MAX_Y = 2. ;

Event::Event(Double_t rho, Double_t max_y) :
    rho(rho),
    max_y(max_y)
{}

/*Event::Event(TTree * tree)
{

}*/

Double_t Event::f(Double_t r)
{
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

// Rejection method to generate radius r
Double_t Event::r_generate(Double_t rho)
{
    //TRandom r(0);
    //TF1 * g = new TF1("g", g, rho, TMath::Infinity());
    //Double_t result = g->GetRandom();
    Double_t R = gRandom->Uniform(0., 1.);
    Double_t result = 1. / TMath::Sqrt(TMath::Exp((1. - R) * TMath::Log(1. + 1. /(rho * rho)))- 1. );

    Double_t temp = gRandom->Uniform(0., 1.);
    while (temp > f(result) / g(result) )//* 3. / (5. * TMath::Pi()))
    {
        //result = g->GetRandom();
        R = gRandom->Uniform(0., 1.);
        result = 1. / TMath::Sqrt(TMath::Exp((1. - R) * TMath::Log(1. + 1. /(rho * rho)))- 1.);
        temp = gRandom->Uniform(0., 1.);
    }
    return result;
}

// Boundary for theta approaching zero
Double_t Event::phi(Double_t r)
{
    if (r > 1/2)
    {
        return TMath::ACos(1/(2*r));
    }
    return 0;
}

// Distribution for theta
Double_t Event::theta(Double_t r)
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

// For Rapidity distribution - lifetime
Double_t Event::lambda(Double_t x01)
{
    Double_t resultat = 0. ;
    Double_t rho2 = rho / x01; // Scale first

    // Compute for rho <= 1/2
    if (rho2 < 1. / 2. ) {
        resultat = TMath::Log(1. / 3. * (1. / (rho2 * rho2) - 1. ));
    }
    // Compute for rho > 1/2
    Double_t borne = (rho2 < 1. / 2.) ? 1. / 2. : rho2;
    TF1 * integral_function = new TF1("integral_function", "1. /(x * TMath::Abs(1. -x*x)) * TMath::ATan(TMath::Abs(1. -x)/(1. +x) * TMath::Sqrt((x+1. /2. )/(x-1. /2. )))", borne, TMath::Infinity());
    //integral_function->DrawF1(1., 100.);
    resultat = resultat + 4. / TMath::Pi() * integral_function->Integral(borne, TMath::Infinity());

    // Update cache
    RHO = rho;
    LAMBDA = resultat;

    return resultat;
}

/* Creates and saves in __filename__ a lookup table
 * Logarithmic scale
 */
void Event::WriteLookupTable(const char * filename)
{
    int n = 320;

    Double_t x01 = 0.0;
    Double_t step = 0.0000000000000000000000000000001;
    Double_t l;
    std::ofstream lut(filename);
    if (lut.is_open())
    {
        lut << rho << " " << n << "\n";
        for (int i = 1; i <= n; ++i)
        {
            if (i%10 == 0) step *= 10. ;
            x01 += step;
            l = lambda(x01);
            lut << x01 << " " << l << "\n";
        }
        lut.close();
    }
}

/* Loads the lookup table from __filename__ and put it into Interpolator
 */
void Event::LoadLookupTable(const char * filename)
{
    std::cerr << "Reading lookup table..." << std::endl;
    std::ifstream lut(filename);
    Double_t rho;
    int n;
    std::vector<Double_t> x01, l;
    if (lut.is_open())
    {
        lut >> rho >> n;
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
    std::cerr << "\033[1;32m Done. \033[0m" << std::endl;
}

void Event::PrintLookupTable()
{
    std::cerr << "Printing lookup table..." << std::endl;
    for (auto const& x: lookup_table)
    {
        std::cerr << x.first << " " << x.second << std::endl;
    }
    std::cerr << "done." << std::endl;    
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
        LoadLookupTable("lookup_table");
    }
    // Interpolate
    if (x01_min > x01 || x01 > x01_max)
    {
        std::cerr << "\033[1;31mWarning : x01 out of lookup table range.\033[0m \n Adding entry for " << x01 << std::endl;
        // Add new points
        lookup_table.insert(std::pair<Double_t, Double_t>(x01, lambda(x01)));
        x01_min = TMath::Min(x01_min, x01);
        x01_max = TMath::Max(x01_max, x01);
        SetInterpolatorData();
    }
    //std::cerr << rho << " " << x01 << " " << interpolator.Eval(x01) << std::endl;
    return interpolator.Eval(x01);
}

// Generate rapidity
Double_t Event::y_generate(Double_t x01)
{
    gRandom->SetSeed();
    Double_t lbda = getLambda(x01);
    Double_t R = gRandom->Uniform(0., 1. / lbda);
    // std::cerr << R << " " << lbda << std::endl;
    return -1. / lbda * TMath::Log(1. - lbda * R);
    // FIXME borne sup de la fonction ? TMath::Infinity() ?
    /*TF1 * y_distribution = new TF1("y_distribution", "exp(- [0] * x)", 0. , 1000.);
    y_distribution->SetParameter(0, lbda);
    return y_distribution->GetRandom();*/
}

// Draw a single gluon
void Event::draw(Double_t x, Double_t y, Double_t rapidity)
{
    TEllipse * ellipse = new TEllipse(x, y, 0.005 , 0.005);
    //Float_t transparency = 1. - rapidity / 10.;

    //std::cerr << rapidity << " " << transparency << std::endl;
    ellipse->SetFillColor(kCyan + (int) TMath::Ceil(10. * rapidity));
    //ellipse->SetFillColorAlpha(4, transparency);
    ellipse->Draw();
    //C->Update();
}

void Event::draw_tree(TTree * tree)
{
    //TFile * f = new TFile("tree.root");
    //tree->MakeClass("");
    //tree->Draw("rapidity", "", "");
}

bool Event::generate(Double_t rho, Dipole * dipole, Dipole * dipole1, Dipole * dipole2, Double_t max_y)
{
    // FIXME use SetSeed() here or not ?
    //gRandom->SetSeed();
    // First generate rho and theta in normalized referential
    Double_t r = r_generate(rho);
    Double_t t = theta(r);
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
    //std::cerr << dipole1->phi << " " << dipole2->phi << std::endl;
    
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
    Double_t scale_factor = dipole->radius * 2.;
    //std::cerr << dipole->radius << " " << scale_factor << std::endl;
    dipole1->coord *= scale_factor;
    dipole2->coord *= scale_factor;
    dipole1->radius *= scale_factor;
    dipole2->radius *= scale_factor;

    // Rotate and stay in 0 - 2pi
    dipole1->coord = TVector2(dipole1->coord.X()-dipole->radius, dipole1->coord.Y()).Rotate(dipole->phi) + TVector2(dipole->radius, 0);
    dipole2->coord = TVector2(dipole2->coord.X()-dipole->radius, dipole2->coord.Y()).Rotate(dipole->phi) + TVector2(dipole->radius, 0);

    //dipole1->phi = dipole1->coord.Phi_0_2pi(dipole1->phi + dipole->phi);
    //dipole2->phi = dipole2->coord.Phi_0_2pi(dipole2->phi + dipole->phi);
    dipole1->phi = dipole1->phi + dipole->phi;
    dipole2->phi = dipole2->phi + dipole->phi;

    // Translation
    TVector2 translation_vector = dipole->coord - TVector2(dipole->radius, 0.0);
    //translation_vector.Print();
    dipole1->coord += translation_vector;
    dipole2->coord += translation_vector;

    return true;
}

TTree * Event::make_tree(const char * filename, const char * treename, bool draw, bool draw_step_by_step)
{
    gRandom->SetSeed(); // /!\ IMPORTANT or we always get the same values
    Dipole dipole(0,0);

    std::queue<Dipole> dipoles;
    dipoles.push(dipole);

    TTree * tree = new TTree(treename, "Dipole splitting");
    tree->Branch("rapidity", &dipole.rapidity, "rapidity/D");
    tree->Branch("radius", &dipole.radius, "radius/D");
    tree->Branch("phi", &dipole.phi, "phi/D");
    tree->Branch("coord", "TVector2", &dipole.coord);
    tree->Branch("depth", &dipole.depth, "depth/L");
    tree->Branch("index_children", &dipole.index_children, "index_children/L");
    tree->Branch("index_parent", &dipole.index_parent, "index_parent/L");
    tree->Branch("isLeaf", &dipole.isLeaf, "isLeaf/O");

    int i = 0;
    Long64_t index = -1; // index of the dipole at current depth

    TCanvas * C = new TCanvas("C", "C", 0, 0, 3000, 2000);
    C->cd(1);
    if (draw || draw_step_by_step)
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
                std::cerr << "Warning : overflow on index value" << std::endl;
                break;
            }
        }
        dipole.index = i;
        dipole.nb_left_brothers_split = current_depth_dipoles_split;
        dipole.nb_right_brothers = current_depth_dipoles - index - 1;

        Long64_t children_index = dipole.index + dipole.nb_right_brothers + 2 * dipole.nb_left_brothers_split + 1;
        //std::cerr << dipole.index << " " << dipole.nb_left_brothers_split << " " << dipole.nb_right_brothers << " " << children_index << std::endl;
        
        Dipole dipole1(dipole.depth + 1, children_index), 
               dipole2(dipole.depth + 1, children_index + 1);
        bool success = generate(rho, &dipole, &dipole1, &dipole2, max_y);
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

        if ((draw && dipole.isLeaf) || draw_step_by_step)
        {
            dipole.Draw();
        }
                  
        // FIXME for some reason the raw values differ slightly from the tree scan ! 
        //std::cerr << dipole1.radius << " " << dipole1.coord.X() << " " << dipole1.coord.Y() <<  std::endl;
        ++i;
        //max_depth = depth;
        tree->Fill(); // Stores dipole (parent)
        // if step by step
        if (draw_step_by_step)
        {
            gPad->Modified();
            gPad->Update();
            C->WaitPrimitive();
            //gSystem->Sleep(100);
            //gSystem->ProcessEvents();
            std::cin.ignore(); // Wait for keypress before continuing
        }
    }
    if (draw || draw_step_by_step) C->Update();

    std::cerr << i << " dipôles générés" << std::endl;

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
    return tree;
}

/********************************************* NORMALIZED + FIT *************************/

// Split 1 dipole - to test the distribution
void Event::generate_normalized(Double_t rho, Double_t * x, Double_t * y, Double_t * rapidity)
{
    // Generate r, theta and y
    Double_t r = r_generate(rho);
    Double_t t = theta(r);
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
    //std::cerr << *x << " " << *y << std::endl;
}

// Test one generation : visualize probability distribution of the gluon at splitting
void Event::bare_distribution()
{
    bool DRAW_ELLIPSES = false; // If we want different colors, use ellipses
    bool DRAW_STEP_BY_STEP = true;
    int number_occurrences = 300000;
    Double_t rho = 0.2;

    Double_t x[number_occurrences], y[number_occurrences], rapidity[number_occurrences];

    TCanvas * C = new TCanvas("C", "C", 0, 0, 1024, 768);
    if (!DRAW_STEP_BY_STEP) C->Divide(2, 1, 0.05, 0.05);

    C->cd(1);
    gPad->SetTitle("QCD");
    gPad->DrawFrame(-1.0, -0.8, 2.0, 0.8, "Gluons");
    
    TH1D * hist = new TH1D("hist", "Plot X", 100, -1, 2);
    Double_t hist_y = 0.04;
    Double_t margin = 0.01;

    //TF1 * r_distribution = new TF1("r_distribution", r_distribution, rho, 1);
    //TF1 * theta_distribution = new TF1("theta_distribution", theta, 0, TMath::Pi());
    // Generate gluons according to dP
    for (int i = 0; i < number_occurrences; ++i)
    {
        generate_normalized(rho, &x[i], &y[i], &rapidity[i]);
        if (DRAW_ELLIPSES || DRAW_STEP_BY_STEP) draw(x[i], y[i], rapidity[i]);
        if (DRAW_STEP_BY_STEP)
        {
            C->Modified();
            C->Update();
            C->WaitPrimitive();
            std::cin.ignore();
        }
        if (y[i] > (hist_y - margin) && y[i] < (hist_y + margin)) hist->Fill(x[i]);
    }
    C->Update();

    if (!DRAW_ELLIPSES && !DRAW_STEP_BY_STEP)
    {
        TGraph * gluons = new TGraph(number_occurrences, x, y);
        gluons->SetTitle("P3A");
        gluons->Draw("AP");
        gluons->SetMarkerColor(2);
        gluons->GetXaxis()->SetRangeUser(-1., 2.);
        gluons->GetYaxis()->SetRangeUser(-0.8, 0.8);
        gluons->Draw("AP");
        C->Update();
    }

    C->cd(2);
    gPad->SetLogy();
    hist->Draw("E1");

    TF1 * p = new TF1("p", "[0] / ((x*x + [1]*[1])*((x-1)*(x-1)+[1]*[1]))", -0.5, 1.5);
    p->FixParameter(1, hist_y);
    hist->Fit("p", "R");

    Double_t chi2 = p->GetChisquare();
    std::cerr << "Chi2 value : " << chi2 << std::endl;   
}

void Event::fit_r()
{
    int number_occurrences = 1000000;
    Double_t rho = 0.01;

    TCanvas * C = new TCanvas("C", "Fit r", 0, 0, 1024, 768);
    TH1D * hist = new TH1D("hist", "hist", 400, 0, 2);
    gPad->SetLogy();

    for (int i = 0 ; i < number_occurrences ; ++i)
    {
        hist->Fill(r_generate(rho));
    }
    hist->Draw("E1");

    TF1 * fr1 = new TF1("fr1", "[0] / (x*(1. - x*x))", 0, 0.5);
    hist->Fit("fr1", "R");

    TF1 * fr2 = new TF1("fr2", "2. * [0] / (x*TMath::Abs(1. - x*x)) * TMath::ATan(TMath::Abs(1. - x)/(1. + x) * TMath::Sqrt((x+1./2.)/(x-1./2.)))", 0.5, 2);
    hist->Fit("fr2", "R+");

    Double_t chi2_1 = fr1->GetChisquare();
    Double_t chi2_2 = fr2->GetChisquare();
    hist->SetTitle(TString::Format("#splitline{Distribution de la taille r du dipole (#rho = %.12g)}{%d bins - chi2 = %.12g et %.12g}", rho, hist->GetSize()-2, chi2_1, chi2_2));
}

void Event::fit_y()
{
    int number_occurrences = 1000000;
    Double_t rho = 0.01;

    TCanvas * C = new TCanvas("C", "Fit y", 0, 0, 1024, 768);
    TH1D * hist = new TH1D("hist", "hist", 100, 0, 1);
    gPad->SetLogy();

    for (int i = 0 ; i < number_occurrences ; ++i)
    {
        hist->Fill(y_generate(1.0));
    }
    hist->Draw("E1");

    TF1 * fy = new TF1("fy", "[0] * exp(- [1] * x)", 0, 1);
    fy->FixParameter(1, LAMBDA);
    hist->Fit("fy", "R");

    Double_t chi2 = fy->GetChisquare();
    hist->SetTitle(TString::Format("#splitline{Distribution de la rapidite y du dipole (#rho = %.12g)}{%d bins - chi2 = %.12g}", rho, hist->GetSize()-2, chi2));  
}

void Event::fit_x()
{
    int number_occurrences = 1000000;
    Double_t rho = 0.01;

    Double_t x, y, rapidity;

    TCanvas * C = new TCanvas("C", "C", 0, 0, 1024, 768);
    
    Double_t hist_y = 0.2;
    Double_t margin = 0.005;
    TH1D * hist = new TH1D("hist", "hist", 100, -0.5, 1.5);

    //TF1 * r_distribution = new TF1("r_distribution", r_distribution, rho, 1);
    //TF1 * theta_distribution = new TF1("theta_distribution", theta, 0, TMath::Pi());
    // Generate gluons according to dP
    for (int i = 0; i < number_occurrences; ++i)
    {
        generate_normalized(rho, &x, &y, &rapidity);
        if (y > (hist_y - margin) && y < (hist_y + margin)) hist->Fill(x);
    }

    gPad->SetLogy();
    hist->Draw("E1");

    TF1 * p = new TF1("p", "[0] / ((x*x + [1]*[1])*((x-1)*(x-1)+[1]*[1]))", -0.3, 1.3);
    p->FixParameter(1, hist_y);
    hist->Fit("p", "R");

    Double_t chi2 = p->GetChisquare();
    hist->SetTitle(TString::Format("#splitline{Distribution of projection X with Y = %.12g +/- %.12g - #rho = %.12g}{%d bins - chi2 = %.12g}", hist_y, margin, rho, hist->GetSize()-2, chi2)); 
}