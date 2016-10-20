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

Double_t RHO = 0.;
Double_t LAMBDA = 0.;
Double_t MAX_Y = 2. ;

Event::Event() :
    rho(RHO),
    max_y(MAX_Y)
{}

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
    //gRandom->SetSeed();
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
Double_t Event::lambda(Double_t rho)
{
    if (rho == RHO)
    {
        return LAMBDA;
    }
    Double_t resultat = 0. ;
    // Compute for rho <= 1/2
    
    if (rho < 1. / 2. ) {
        resultat = TMath::Log(1. / 3. * (1. / (rho * rho) - 1. ));
    }
    // Compute for rho > 1/2
    Double_t borne = (rho < 1. / 2.) ? 1. / 2. : rho;
    TF1 * integral_function = new TF1("integral_function", "1. /(x * TMath::Abs(1. -x*x)) * TMath::ATan(TMath::Abs(1. -x)/(1. +x) * TMath::Sqrt((x+1. /2. )/(x-1. /2. )))", borne, TMath::Infinity());
    //integral_function->DrawF1(1., 100.);
    resultat = resultat + 4. / TMath::Pi() * integral_function->Integral(borne, TMath::Infinity());

    // Update cache
    RHO = rho;
    LAMBDA = resultat;

    return resultat;
}

// Generate rapidity
Double_t Event::y_generate(Double_t rho)
{
    Double_t lbda = lambda(rho);
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
    Double_t x, y;
    Double_t r = r_generate(rho);
    Double_t t = theta(r);
    Double_t rapidity = y_generate(rho);   
    Double_t quadrant = gRandom->Uniform(0., 1.);

    if (rapidity > max_y) return false;

    dipole1->rapidity = dipole->rapidity + rapidity;
    dipole2->rapidity = dipole->rapidity + rapidity;
    dipole1->radius = r / 2. ;

    // Set coordinates of centers in normalized referential
    // Use cartesian coordinates instead of SetMagPhi because of the four quadrants
    if (quadrant < 0.25)
    {
        dipole1->coord.Set(r * TMath::Cos(t), r * TMath::Sin(t));
    }
    else if (quadrant < 0.5)
    {
        dipole1->coord.Set(r * TMath::Cos(t), - r * TMath::Sin(t));
    }
    else if (quadrant < 0.75)
    {
        dipole1->coord.Set(1.0 - r * TMath::Cos(t), r * TMath::Sin(t));           
    }
    else
    {
        dipole1->coord.Set(1.0 - r * TMath::Cos(t), - r * TMath::Sin(t));           
    }
    TVector2 temp(-1.0, 0.0);
    dipole2->coord = temp + dipole1->coord;
    dipole2->radius = dipole2->coord.Mod() / 2.;
    dipole1->coord = dipole1->coord / 2. + temp / 2.;
    dipole2->coord = dipole2->coord / 2. - temp / 2.;
    dipole1->phi = t;
    dipole2->phi = dipole2->coord.Phi();

    // Now we have dipole1 and dipole2 in the normalized referential, get back
    // Translation
    dipole1->coord += dipole->coord;
    dipole2->coord += dipole->coord;

    // Scale
    dipole1->coord *= dipole->radius * 2.;
    dipole2->coord *= dipole->radius * 2.;

    // Rotate and stay in 0 - 2pi
    dipole1->coord = dipole1->coord.Rotate(dipole->phi);
    dipole2->coord = dipole2->coord.Rotate(dipole->phi);
    dipole1->phi = dipole1->coord.Phi_0_2pi(dipole1->phi + dipole->phi);
    dipole2->phi = dipole2->coord.Phi_0_2pi(dipole2->phi + dipole->phi);

    // Set references parent/child
    /*dipole->child1 = dipole1;
    dipole->child2 = dipole2;
    dipole1->parent = dipole;
    dipole2->parent = dipole;*/
    return true;
}

// Split 1 dipole - to test the distribution
void Event::generate_normalized(Double_t rho, Double_t * x, Double_t * y, Double_t * rapidity)
{
    // Generate r, theta and y
    Double_t r = r_generate(rho);
    Double_t t = theta(r);
    *rapidity = y_generate(rho);

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
    int number_occurrences = 300000;
    Double_t rho = 0.01;

    Double_t x[number_occurrences], y[number_occurrences], rapidity[number_occurrences];

    TCanvas * C = new TCanvas("C", "C", 0, 0, 1024, 768);
    C->Divide(2, 1, 0.05, 0.05);

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
        if (DRAW_ELLIPSES) draw(x[i], y[i], rapidity[i]);
        if (y[i] > (hist_y - margin) && y[i] < (hist_y + margin)) hist->Fill(x[i]);
    }
    C->Update();

    if (!DRAW_ELLIPSES)
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
        hist->Fill(y_generate(rho));
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

void Event::make_tree(const char * filename)
{
    Dipole dipole(0,0);

    std::queue<Dipole> dipoles;
    dipoles.push(dipole);

    TFile * output = new TFile(filename, "recreate");
    TTree * tree = new TTree("T", "Gluons generated");
    tree->Branch("rapidity", &dipole.rapidity, "rapidity/D");
    tree->Branch("radius", &dipole.radius, "radius/D");
    tree->Branch("phi", &dipole.phi, "phi/D");
    tree->Branch("coord", "TVector2", &dipole.coord);
    tree->Branch("depth", &dipole.depth, "depth/L");
    tree->Branch("index", &dipole.index, "index/L");
    tree->Branch("isLeaf", &dipole.isLeaf, "isLeaf/O");
    //tree->Branch("parent", &dipole.parent, "parent/D");
    //tree->Branch("child1", "Dipole", &dipole.child1);
    //tree->Branch("child2", "Dipole", &dipole.child2);

    Double_t rho = 0.01; // Cut-off
    Double_t max_y = 1.5; // Maximal rapidity
    int i = 0;

    Long64_t depth;
    Long64_t index;
    Long64_t max_depth;

    // FIXME remove number_occurrences and keep only max_y
    while (!dipoles.empty())
    {
        //std::cerr << dipoles.size() << std::endl;
        dipole = dipoles.front();
        dipoles.pop();

        if (dipole.depth < std::numeric_limits<Long64_t>::max())
        {
            depth = dipole.depth + 1;
        }
        else
        {
            std::cerr << "Warning : overflow on depth value" << std::endl;
            break;
        }

        if (dipole.index < std::numeric_limits<Long64_t>::max())
        {
            index = dipole.index * 2;
        }
        else
        {
            std::cerr << "Warning : overflow on index value" << std::endl;
            break;
        }
        
        Dipole dipole1(depth, index), dipole2(depth, index+1);
        bool success = generate(rho, &dipole, &dipole1, &dipole2, max_y);
        if (success)
        {
            dipoles.push(dipole1);
            dipoles.push(dipole2);
        }
        else
        {
            dipole.isLeaf = true;
        }
        // FIXME for some reason the raw values differ slightly from the tree scan ! 
        //std::cerr << dipole1.radius << " " << dipole1.coord.X() << " " << dipole1.coord.Y() <<  std::endl;
        ++i;
        max_depth = depth;
        tree->Fill(); // Stores dipole
    }
    --max_depth;// because the last depth is for leaves that we didn't stored

    TVector * v = new TVector(1);
    //v->SetName("max_depth");
    v[0] = max_depth;
    tree->GetUserInfo()->Add(v);
    //tree->Print();
    // Print entries
    //tree->Scan("rapidity:radius:phi:coord.X():coord.Y():depth:index", "depth == 34");
    //tree->Draw("radius");
    output->Write();
}