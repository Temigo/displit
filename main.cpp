#include "main.h"
#include "utils.h"
#include "event.h"

#include <TApplication.h>
#include<TMath.h>
#include <TROOT.h>
#include <sstream>
#include <iostream>

int main( int argc, const char* argv[] )
{
    TApplication * myapp = new TApplication("myapp", 0, 0);
    
    //connect(myapp,SIGNAL(lastWindowClosed()),TQtRootSlot::CintSlot(),SLOT(TerminateAndQuit());
    //p3a();

    std::istringstream iss( argv[1] );
    int val = 0; // By default only read file
    int nb_events = 1000; // Number of events to read/generate
    Double_t rho = 0.000001;
    Double_t max_y = 1.0;

    /*int nthreads = 4;
    ROOT::EnableImplicitMT(nthreads);
    ROOT::EnableThreadSafety();*/

    if (!(iss >> val))
    {
        std::cerr << "Not valid argument." << std::endl;
        return EXIT_FAILURE;
    }

    if (val)
    {
        generate_events(nb_events, rho, max_y);
    }

    /*TMinuit minuit(1); // Instantiate Minuit for 3 parameters
    // Or use TFitter ?
    minuit.SetFCN(minuitFunction);

    // Define parameters
    minuit->SetParameter(0, "X", 3, 1, 0, 0);

    // Run minimizer
    // SIMPLEX : Robustly finds the nearest local minimum, but precision isn't very high
    // MIGRAD : Precisely find the minimum, but you need to be careful since it assumes the inputs are close to the minimum
    minuit->ExecuteCommand("SIMPLEX", 0, 0);
    Double_t bestX = minuit->GetParameter(0);*/

    // General fit
    //general_plot(myapp);
    //stat_events(myapp, nb_events, max_y, 1.0);
    //CommonAncestorPlot(myapp, nb_events, max_y, 1.0, rho);
    //fluctuations(myapp, nb_events, max_y, 1.0, rho);

    // Compute biggest children
    /*TFile f("tree.root");
    TTree * tree;
    f.GetObject("tree0", tree); 
    //tree->MakeClass("EventTree");
    compute_biggest_child(tree, myapp);*/

    //TF2 * cutoff = new TF2("cutoff", "exp(-[0] / (2 * [1]^2) * (1 + 2*x^2 -2*x*cos(y)))", 0, TMath::Infinity(), 0, TMath::Pi());
    //Event * e = new Event(rho, max_y, "lookup_table", cutoff, true);
    Event * e = new Event(rho, max_y, "lookup_table");
    //e->WriteLookupTable("lookup_table");
    //e->LoadLookupTable();
    //e->PrintLookupTable();
    std::cerr << e->getLambda(25) << std::endl;
    //e->make_tree("tree.root", "tree", false, true);
    //e->bare_distribution();
    //std::cerr << e->getLambda(1.0) << std::endl;
    /*Double_t r = e->r_generate(true);
    std::cerr << r << std::endl;
    std::cerr << e->theta(r) << std::endl;
    std::cerr << e->lambda(1.0) << std::endl;*/

    myapp->Run();

    return EXIT_SUCCESS;
}