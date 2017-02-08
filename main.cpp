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
#include <fstream>

Double_t heaviside(Double_t * x, Double_t * p)
{
    std::cerr << x[0] << " " << x[1] << std::endl;
    return (x[0] > 2 ? 0.0 : 1.0);
}

int main( int argc, char* argv[] )
{
    TApplication * myapp = new TApplication("myapp", &argc, argv);
    gRandom = new TRandom3(0);

    std::istringstream iss( argv[1] );
    //std::istringstream iss2 ( argv[2] );
    int val = 0; // By default only read file
    //char * tree_file = argv[2];

    /*std::vector<std::tuple<int, Double_t, Double_t, std::string > > parameters = {
        std::make_tuple(0, 0.01, 2.0, "gaussian"),
    };*/

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

        int nb_events; // Number of events to read/generate
        Double_t rho;
        Double_t max_y;
        std::string cutoff_type;
        unsigned int len; // string length

        if (rank == 0)
        {
            std::cerr << size << " processes" << std::endl;
            std::ifstream params("parameters");
            std::cout << "Hi";
            if (params.is_open())
            {
                int i = 0;
                while (params >> nb_events >> rho >> max_y >> cutoff_type)
                {
                    ++i;
                    len = cutoff_type.length();
                    std::cout << "Sending to process " << i << std::endl;
                    MPI_Send(&len, 1, MPI_UNSIGNED, i, 0, MPI_COMM_WORLD);
                    MPI_Send(&nb_events, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                    MPI_Send(&rho, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
                    MPI_Send(&max_y, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
                    MPI_Send(cutoff_type.c_str(), cutoff_type.length(), MPI_CHAR, i, 0, MPI_COMM_WORLD);
                }
                params.close();
            }            
        }
        else
        {
            unsigned int len;
            MPI_Recv(&len, 1, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            std::vector<char> cutoff_type_temp(len);
            MPI_Recv(&nb_events, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&rho, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&max_y, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(cutoff_type_temp.data(), len, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            cutoff_type.assign(cutoff_type_temp.begin(), cutoff_type_temp.end());
            std::cout << "Process " << rank << " with " << nb_events << " events, rho = " << rho << ", ymax = " << max_y << ", cutoff " << cutoff_type << std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);

        // Gaussian cutoff
        TF2 * cutoff_gaussian = new TF2("cutoff_gaussian", "exp(-[0] / (2 * [1]^2) * (1 + 2*x^2 -2*x*cos(y)))", 0, TMath::Infinity(), 0, TMath::Pi());
        // Lorentzian cutoff
        TF2 * cutoff_lorentzian = new TF2("cutoff_lorentzian", "1 / (1 + ([0]^2 / (2 * [1]^2) * (1 + 2*x^2 -2*x*cos(y))))", 0, TMath::Infinity(), 0, TMath::Pi());
        // Maxwellian cutoff
        //TF2 * cutoff = new TF2("cutoff", "1 / (1 + exp([0]^2 / (2 * [1]^2) * (1 + 2*x^2 -2*x*cos(y))))", 0, TMath::Infinity(), 0, TMath::Pi());
        // Heaviside
        TF2 * cutoff_heaviside = new TF2("cutoff_heaviside", heaviside, 0, TMath::Infinity(), 0, TMath::Pi(), 2);
        //TF1 * cutoff = new TF1("cutoff", "(x > 2) ? 0.0 : 1.0", 0, TMath::Infinity());
        //std::cerr << cutoff->Eval(3, 2.5);
        //std::cerr << cutoff->Integral(0, 77, 0, TMath::Pi());
        //cutoff->SetParameter(0, 1.0);
        //cutoff->SetParameter(1, 2.0);
        //cutoff->Draw("surf2");

        std::map<std::string, TF2 *> cutoffs = {
            {"gaussian", cutoff_gaussian},
            {"lorentzian", cutoff_lorentzian},
            {"rigid", cutoff_heaviside}
        };
        cutoffs[cutoff_type]->SetName("cutoff");

        if (rank > 0)
        {
            // Set up Filenames
            std::string s = "mpi_tree_" + std::to_string(nb_events) + "events_cutoff" + std::to_string(rho) + "_ymax" + std::to_string(max_y) + "_" + cutoff_type + ".root";
            const char * tree_file = s.c_str();
            std::string s2 = "lookup_table_gaussian_cutoff10-" + std::to_string(rank);
            const char * lut_file = s2.c_str();

            try
            {
                generate_events(nb_events, rho, max_y, true, cutoffs[cutoff_type], false, tree_file, lut_file);
            }
            catch (...)
            {
                return EXIT_FAILURE;
            }
        }
        delete cutoff_gaussian;
        delete cutoff_lorentzian;
        delete cutoff_heaviside;

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