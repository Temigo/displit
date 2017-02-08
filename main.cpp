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
#include <RooStats/Heaviside.h>
#include <mpi.h>

Double_t heaviside(Double_t * x, Double_t * p)
{
    std::cerr << x[0] << " " << x[1] << std::endl;
    return (x[0] > 2 ? 0.0 : 1.0);
}

int main( int argc, char* argv[] )
{
    TApplication * myapp = new TApplication("myapp", &argc, argv);
    
    //connect(myapp,SIGNAL(lastWindowClosed()),TQtRootSlot::CintSlot(),SLOT(TerminateAndQuit());

    gRandom = new TRandom3(0);

    std::istringstream iss( argv[1] );
    //std::istringstream iss2 ( argv[2] );
    int val = 0; // By default only read file
    //char * tree_file = argv[2];
    int nb_events = 10; // Number of events to read/generate
    Double_t rho = 0.01;
    Double_t max_y = 3.0;

    //std::cerr << argv[2];
    if (!(iss >> val))
    {
        std::cerr << "Not valid argument." << std::endl;
        return EXIT_FAILURE;
    }

    if (val)
    {
        int rank, size; // rank of the process and number of processes
        MPI_Init(NULL, NULL);

        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        std::cerr << "I am process " << rank << std::endl;
        if (rank == 0) std::cerr << size << " processes" << std::endl;
        //std::cerr << tree_file << std::endl;
        // Gaussian cutoff
        TF2 * cutoff = new TF2("cutoff", "exp(-[0] / (2 * [1]^2) * (1 + 2*x^2 -2*x*cos(y)))", 0, TMath::Infinity(), 0, TMath::Pi());
        // Lorentzian cutoff
        //TF2 * cutoff = new TF2("cutoff", "1 / (1 + ([0]^2 / (2 * [1]^2) * (1 + 2*x^2 -2*x*cos(y))))", 0, TMath::Infinity(), 0, TMath::Pi());
        // Maxwellian cutoff
        //TF2 * cutoff = new TF2("cutoff", "1 / (1 + exp([0]^2 / (2 * [1]^2) * (1 + 2*x^2 -2*x*cos(y))))", 0, TMath::Infinity(), 0, TMath::Pi());
        // Heaviside
        //TF2 * cutoff = new TF2("cutoff", heaviside, 0, TMath::Infinity(), 0, TMath::Pi(), 2);
        //TF1 * cutoff = new TF1("cutoff", "(x > 2) ? 0.0 : 1.0", 0, TMath::Infinity());
        //std::cerr << cutoff->Eval(3, 2.5);
        //std::cerr << cutoff->Integral(0, 77, 0, TMath::Pi());
        //cutoff->SetParameter(0, 1.0);
        //cutoff->SetParameter(1, 2.0);
        //cutoff->Draw("surf2");

        if (rank >= 2 && rank <= 6)
        {
            rho = TMath::Power(10, - rank);
            std::cerr << "Rank : " << rank << " with rho " << rho << std::endl;
        
            std::string s = "mpi_tree_" + std::to_string(nb_events) + "events_cutoff10-" + std::to_string(rank) + "_ymax" + std::to_string(max_y) + "_gaussian.root";
            const char * tree_file = s.c_str();
            std::string s2 = "lookup_table_gaussian_cutoff10-" + std::to_string(rank);
            const char * lut_file = s2.c_str();
            try
            {
                generate_events(nb_events, rho, max_y, true, cutoff, false, tree_file, lut_file);
            }
            catch (...)
            {
                return EXIT_FAILURE;
            }
        }
        delete cutoff;

        MPI_Finalize();
    }

    // General fit
    //general_plot(myapp);
    //stat_events(myapp, max_y, 1.0);
    //CommonAncestorPlot(myapp, nb_events, max_y, 1.0, rho);
    //fluctuations(myapp, nb_events, max_y, 1.0, rho);
    //draw_cutoffs(myapp);

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