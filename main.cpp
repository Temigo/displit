#include "main.h"
#include "utils.h"
#include "event.h"
#include "graphics.h"

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
    int nb_events = 1000; // Number of events to read/generate
    Double_t rho = 0.001;
    Double_t max_y = 3.0;

    if (!(iss >> val))
    {
        std::cerr << "Not valid argument." << std::endl;
        return EXIT_FAILURE;
    }

    if (val)
    {
        // Gaussian cutoff
        //TF2 * cutoff = new TF2("cutoff", "exp(-[0] / (2 * [1]^2) * (1 + 2*x^2 -2*x*cos(y)))", 0, TMath::Infinity(), 0, TMath::Pi());
        // Lorentzian cutoff
        TF2 * cutoff = new TF2("cutoff", "1 / (1 + ([0]^2 / (2 * [1]^2) * (1 + 2*x^2 -2*x*cos(y))))", 0, TMath::Infinity(), 0, TMath::Pi());
        // Maxwellian cutoff
        //TF2 * cutoff = new TF2("cutoff", "1 / (1 + exp([0]^2 / (2 * [1]^2) * (1 + 2*x^2 -2*x*cos(y))))", 0, TMath::Infinity(), 0, TMath::Pi());
        // Heaviside
        //TF2 * cutoff = new TF2("cutoff", "( (-1+x^2 -x*cos(y))> 0 ? 0.0 : 1.0) * y/y", 0, TMath::Infinity(), 0, TMath::Pi());

        /*cutoff->SetParameter(0, 1.0);
        cutoff->SetParameter(1, 2.0);
        cutoff->Draw("surf1");*/
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
    draw_cutoffs(myapp);

    // Compute biggest children
    /*TFile f("tree.root");
    TTree * tree;
    f.GetObject("tree0", tree); 
    //tree->MakeClass("EventTree");
    compute_biggest_child(tree, myapp);*/



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

    //myapp->Run();
    //delete cutoff;

    //fit_bare_r(rho, max_y, myapp);
    
    delete myapp;
    return EXIT_SUCCESS;
}