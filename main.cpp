#include "main.h"
#include "utils.h"
#include "event.h"

#include <TApplication.h>
#include<TMath.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLine.h>
#include <sstream>
#include <iostream>

int main( int argc, char* argv[] )
{
    TApplication * myapp = new TApplication("myapp", &argc, argv);
    
    //connect(myapp,SIGNAL(lastWindowClosed()),TQtRootSlot::CintSlot(),SLOT(TerminateAndQuit());

    gRandom = new TRandom3(0);

    std::istringstream iss( argv[1] );
    int val = 0; // By default only read file
    int nb_events = 10; // Number of events to read/generate
    Double_t rho = 0.01;
    Double_t max_y = 3.0;

    if (!(iss >> val))
    {
        std::cerr << "Not valid argument." << std::endl;
        return EXIT_FAILURE;
    }

    if (val)
    {
        TF2 * cutoff = new TF2("cutoff", "exp(-[0] / (2 * [1]^2) * (1 + 2*x^2 -2*x*cos(y)))", 0, TMath::Infinity(), 0, TMath::Pi());
        try
        {
            generate_events(nb_events, rho, max_y, true, cutoff, false);
        }
        catch (...)
        {
            return EXIT_FAILURE;
        }
        delete cutoff;
    }

    // General fit
    //general_plot(myapp);
    //stat_events(myapp, max_y, 1.0);
    //CommonAncestorPlot(myapp, nb_events, max_y, 1.0, rho);
    //fluctuations(myapp, nb_events, max_y, 1.0, rho);

    // Compute biggest children
    /*TFile f("tree.root");
    TTree * tree;
    f.GetObject("tree0", tree); 
    //tree->MakeClass("EventTree");
    compute_biggest_child(tree, myapp);*/

    // Compare 3 methods
    Event e1(rho, max_y, "lookup_table");
    TF2 * cutoff = new TF2("cutoff", "exp(-[0] / (2 * [1]^2) * (1 + 2*x^2 -2*x*cos(y)))", 0, TMath::Infinity(), 0, TMath::Pi());
    Event e2(rho, max_y, "lookup_table", cutoff, true, true);
    Event e3(rho, max_y, "lookup_table", cutoff, true, false);

    TCanvas canvas("canvas", "Sizes", 1080, 780);
    gPad->SetLogy();
    gStyle->SetPadTickY(1); // y ticks on the right too
    gStyle->SetOptStat(0); // get rid of the statistics box
    
    Double_t max_r = 20;
    int N = 1000000;
    TH1D hist1("hist1", "hist1", 100, 0, max_r);
    TH1D hist2("hist2", "hist2", 100, 0, max_r);
    TH1D hist3("hist3", "hist3", 100, 0, max_r);
    for (int i = 0; i <N; ++i)
    {
        hist1.Fill(e1.r_generate());
        hist2.Fill(e2.r_generate());
        hist3.Fill(e3.r_generate());
    }
    hist1.SetLineWidth(2);
    hist1.SetMarkerStyle(8);
    hist1.SetMarkerSize(0.7);
    hist1.Draw("E1");
    TF1 fr1("fr1", "[0] / (x*(1. - x*x))", 0.02, 0.5);
    TF1 fr2("fr2", "2. * [0] / (x*TMath::Abs(1. - x*x)) * TMath::ATan(TMath::Abs(1. - x)/(1. + x) * TMath::Sqrt((x+1./2.)/(x-1./2.)))", 0.5, 2);
    fr1.SetLineColor(kGreen-6);
    fr2.SetLineColor(kGreen-6);
    hist1.Fit("fr1", "R");
    hist1.Fit("fr2", "R+");

    hist2.SetLineColor(kRed);
    hist2.Draw("E1 same");
    hist3.SetLineColor(kViolet);
    hist3.Draw("E1 same");

    TLegend legend(0.4, 0.6, 0.8, 0.8);
    legend.SetFillColor(0); // white bg
    legend.SetBorderSize(0); // get rid of the box
    legend.SetTextSize(0.045);
    legend.AddEntry(&hist1,"Without cutoff", "L");
    legend.AddEntry(&hist2,"With cutoff - method 1", "L");
    legend.AddEntry(&hist3,"With cutoff - method 2", "L");
    legend.AddEntry(&fr1, "Fit (formula without cutoff)", "L");
    legend.Draw();

    TLine line(rho, 0, rho, hist1.GetMaximum());
    line.SetLineColor(kBlack);
    line.SetLineStyle(7);
    line.Draw();

    //TString::Format("#splitline{Distribution de la taille r du dipole (#rho = %.12g)}{%d bins - chi2 = %.12g et %.12g}", rho, hist.GetSize()-2, chi2_1, chi2_2)
    hist1.SetTitle(TString::Format("%d generations of r with x01 = %.12g and rho = %.12g;Size r of dipoles;Number of dipoles", N, 1.0, rho));
    canvas.Update();

    //Event e(rho, max_y, "lookup_table");
    //e->WriteLookupTable("lookup_table");
    //e->LoadLookupTable();
    //e->PrintLookupTable();
    //std::cerr << e.getLambda(35) << std::endl;
    //e.make_tree("tree", false);
    //e.bare_distribution();
    /*Double_t r = e.r_generate(true);
    std::cerr << r << std::endl;
    std::cerr << e.theta(r) << std::endl;
    std::cerr << e.lambda(1.0) << std::endl;*/
    //e.fit_r();
    //e.theta(1.0, 0.02);

    myapp->Run();
    delete cutoff;
    delete myapp;
    return EXIT_SUCCESS;
}