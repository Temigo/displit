#include <stdio.h>
#include <TMath.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TEllipse.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TApplication.h>

struct Dipole {
    Double_t x; // x coordinate of the middle
    Double_t y; // y coordinate of the middle
    Double_t phi; // angle of the dipole in its own referential
    Double_t rapidity; 
};

Double_t f(Double_t r)
{
    if (r <= 1. / 2. ) {
        return TMath::Pi() / (r * (1. -r*r));
    }
    else {
        return 2. / (r * std::abs(1. -r*r)) * TMath::ATan(TMath::Abs(1. -r) / (1. +r) * TMath::Sqrt((r+1. / 2. )/(r-1. / 2.)));
    }
}

Double_t g(Double_t r)
{
    return 5. * TMath::Pi() / (3. * r * (1. + r * r));
}

// Rejection method to generate radius r
Double_t r_generate(Double_t rho)
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
Double_t phi(Double_t r)
{
    if (r > 1/2)
    {
        return TMath::ACos(1/(2*r));
    }
    return 0;
}

// Distribution for theta
Double_t theta(Double_t r)
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
Double_t lambda(Double_t rho)
{
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
    return resultat;
}

// Generate rapidity
Double_t y_generate(Double_t rho)
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

void p3a()
{
    int number_occurrences = 20000;
    Double_t rho = 0.05;

    Double_t r, t, rapidity;
    Double_t x[number_occurrences], y[number_occurrences];

    // Draw
    TCanvas * C = new TCanvas("C", "C", 0, 0, 1024, 768);
    gPad->SetTitle("QCD");
    gPad->DrawFrame(-1.0, -0.8, 2.0, 0.8, "Gluons");
    

    //TF1 * r_distribution = new TF1("r_distribution", r_distribution, rho, 1);
    //TF1 * theta_distribution = new TF1("theta_distribution", theta, 0, TMath::Pi());
    // Generate gluons according to dP
    for (int i = 0; i < number_occurrences; ++i)
    {
        // Generate r, theta and y
        r = r_generate(rho);
        t = theta(r);
        rapidity = y_generate(rho);
        //r = r_distribution->GetRandom();
        //t = theta_distribution->GetRandom();
        //Double_t p = dP(r, t, y, rho);
        // Compute dP
        Double_t quadrant = gRandom->Uniform(0., 1.);
        if (quadrant < 0.25)
        {
            x[i] = r * TMath::Cos(t);
            y[i] = r * TMath::Sin(t);
        }
        else if (quadrant < 0.5)
        {
            x[i] = r * TMath::Cos(t);
            y[i] = - r * TMath::Sin(t);
        }
        else if (quadrant < 0.75)
        {
            x[i] = 1.0 - r * TMath::Cos(t);
            y[i] = r * TMath::Sin(t);            
        }
        else
        {
            x[i] = 1.0 - r * TMath::Cos(t);
            y[i] = - r * TMath::Sin(t);            
        }
        //std::cerr << r << " " << t << " " << R << std::endl;
        //std::cerr << x << " " << y << std::endl;
        
        TEllipse * ellipse = new TEllipse(x[i] , y[i] , 0.005 , 0.005);
        //Float_t transparency = 1. - rapidity / 10.;

        //std::cerr << rapidity << " " << transparency << std::endl;
        ellipse->SetFillColor(kCyan + (int) TMath::Ceil(10. * rapidity));
        //ellipse->SetFillColorAlpha(4, transparency);
        ellipse->Draw();
        //C->Update();
    }
    C->Update();

    /*TGraph * gluons = new TGraph(number_occurrences, x, y);
    gluons->SetTitle("P3A");
    gluons->Draw("AP");
    gluons->SetMarkerColor(2);
    gluons->GetXaxis()->SetRangeUser(-1., 2.);
    gluons->GetYaxis()->SetRangeUser(-0.8, 0.8);
    gluons->Draw("AP");
    C->Update();*/
}
int main( int argc, const char* argv[] )
{
    TApplication * myapp = new TApplication("myapp", 0, 0);
    //connect(myapp,SIGNAL(lastWindowClosed()),TQtRootSlot::CintSlot(),SLOT(TerminateAndQuit());
    p3a();
    myapp->Run();
    return 0;
}
