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
        return 2. / (r * std::abs(1. -r*r)) * TMath::ATan(TMath::Abs(1. -r) / (1. +r) * TMath::Sqrt((r+1. / 2. )/(r-1. / 2.)));
    }
}

Double_t Event::g(Double_t r)
{
    return 5. * TMath::Pi() / (3. * r * (1. + r * r));
}

// Rejection method to generate radius r
Double_t Event::r_generate(Double_t rho)
{
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

// Draw gluon
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

void Event::generate(Double_t rho, Dipole * dipole, Dipole * dipole1, Dipole * dipole2)
{
    Double_t x, y;
    Double_t r = r_generate(rho);
    Double_t t = theta(r);
    Double_t rapidity = y_generate(rho);    
    Double_t quadrant = gRandom->Uniform(0., 1.);

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
    dipole->child1 = dipole1;
    dipole->child2 = dipole2;
    dipole1->parent = dipole;
    dipole2->parent = dipole;
}

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

void Event::bare_distribution()
{
    bool DRAW_ELLIPSES = false; // If we want different colors, use ellipses
    int number_occurrences = 20000;
    Double_t rho = 0.05;

    Double_t x[number_occurrences], y[number_occurrences], rapidity[number_occurrences];

    TCanvas * C = new TCanvas("C", "C", 0, 0, 1024, 768);
    gPad->SetTitle("QCD");
    gPad->DrawFrame(-1.0, -0.8, 2.0, 0.8, "Gluons");
    
    //TF1 * r_distribution = new TF1("r_distribution", r_distribution, rho, 1);
    //TF1 * theta_distribution = new TF1("theta_distribution", theta, 0, TMath::Pi());
    // Generate gluons according to dP
    for (int i = 0; i < number_occurrences; ++i)
    {
        generate_normalized(rho, &x[i], &y[i], &rapidity[i]);
        if (DRAW_ELLIPSES) draw(x[i], y[i], rapidity[i]);
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
}

void Event::make_tree(const char * filename)
{
    Dipole dipole;

    std::queue<Dipole> dipoles;
    dipoles.push(dipole);

    TFile * output = new TFile(filename, "recreate");
    TTree * tree = new TTree("T", "Gluons generated");
    tree->Branch("rapidity", &dipole.rapidity, "rapidity/D");
    tree->Branch("radius", &dipole.radius, "radius/D");
    tree->Branch("phi", &dipole.phi, "phi/D");
    tree->Branch("coord", "TVector2", &dipole.coord);
    //tree->Branch("parent", &dipole.parent, "parent/D");
    //tree->Branch("child1", "Dipole", &dipole.child1);
    //tree->Branch("child2", "Dipole", &dipole.child2);

    int number_occurrences = 20000;
    Double_t rho = 0.05; // Cut-off
    Double_t max_y = 10.0; // Maximal rapidity
    int i = 0;

    // FIXME remove number_occurrences and keep only max_y
    while (i < number_occurrences && !dipoles.empty())
    {
        //std::cerr << dipoles.size() << std::endl;
        dipole = dipoles.front();
        dipoles.pop();
        tree->Fill();

        Dipole dipole1, dipole2;
        generate(rho, &dipole, &dipole1, &dipole2);
        if (dipole1.rapidity <= max_y) dipoles.push(dipole1);
        if (dipole2.rapidity <= max_y) dipoles.push(dipole2);
        // FIXME for some reason the raw values differ slightly from the tree scan ! 
        //std::cerr << dipole1.radius << " " << dipole1.coord.X() << " " << dipole1.coord.Y() <<  std::endl;
        ++i;
    }

    //tree->Print();
    // Print entries
    //tree->Scan("rapidity:radius:phi:coord.X():coord.Y():child1:child2:parent");
    //tree->Draw("radius");
    output->Write();
}


