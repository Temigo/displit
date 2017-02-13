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
#include <TFile.h>

#include <sstream>
#include <iostream>
#include <RooStats/Heaviside.h>
#include <mpi.h>
#include <fstream>
#include <regex>

Double_t heaviside(Double_t * x, Double_t * p)
{
    return (x[0] > 2 ? 0.0 : 1.0);
}

std::string encode_parameters(int nb_events, Double_t rho, Double_t max_y, std::string cutoff_type)
{
    return "mpi_tree_" + std::to_string(nb_events) + "events_cutoff" + std::to_string(rho) + "_ymax" + std::to_string(max_y) + "_" + cutoff_type + ".root";
}

void decode_parameters(std::string filename, Double_t * rho, Double_t * max_y, std::string * cutoff_type, int * nb_events)
{
    std::string double_regex = "[-+]?[0-9]*\.?[0-9]+";
    std::regex r("mpi_tree_([[:digit:]]+)events_cutoff("+double_regex+")_ymax("+double_regex+")_(\\w+).root");
    std::smatch m;
    if (std::regex_match(filename, m, r))
    {
        *nb_events = std::stoi(m[1]);
        *rho = std::stod(m[2]);
        *max_y = std::stod(m[3]);
        *cutoff_type = m[4];
    }
    else
    {
        std::cout << "no match" << std::endl;
    }
}

int main( int argc, char* argv[] )
{
    //TApplication * myapp = new TApplication("myapp", &argc, argv);
    TApplication * myapp = new TApplication("myapp", NULL, NULL); // do not pass command line arguments to ROOT
    gRandom = new TRandom3(0);
    std::string val;

    if (argc <= 1)
    {
        std::cerr << "Not valid argument." << std::endl;
        return EXIT_FAILURE;
    }
    else
    {
        val = argv[1];
    }

    if (val == "-h" || val == "--help")
    {
        std::cout << "Usage: ./main generate [parameters file] | fluctuations [file] | check [file] | fit-bare-r [rho] [max_y] | ancestors [nb_events] [rho] [max_y] | draw-cutoffs" << std::endl;
        std::cout << "Options: \n \t --help : Print usage" << std::endl;
    }
    else if (val == "generate")
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
            std::ifstream params(argv[2]);
            if (params.is_open())
            {
                int i = 0;
                while (params >> nb_events >> rho >> max_y >> cutoff_type)
                {
                    ++i;
                    if (i >= size)
                    {
                        std::cout << "Warning : not enough processes. I stop reading parameters here." << std::endl;
                        break;
                    }
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
        else // TODO fix it if process rank too big
        {
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
        TF2 * cutoff_lorentzian1 = new TF2("cutoff_lorentzian1", "1 / (1 + ([0]^2 / (2 * [1]^2) * (1 + 2*x^2 -2*x*cos(y))))", 0, TMath::Infinity(), 0, TMath::Pi());
        TF2 * cutoff_lorentzian2 = new TF2("cutoff_lorentzian2", "1 / (1 + ([0]^2 / (2 * [1]^2) * (1 + 2*x^2 -2*x*cos(y)))^2 )", 0, TMath::Infinity(), 0, TMath::Pi());
        TF2 * cutoff_lorentzian3 = new TF2("cutoff_lorentzian3", "1 / (1 + ([0]^2 / (2 * [1]^2) * (1 + 2*x^2 -2*x*cos(y)))^3 )", 0, TMath::Infinity(), 0, TMath::Pi());
        // Maxwellian cutoff
        //TF2 * cutoff = new TF2("cutoff", "1 / (1 + exp([0]^2 / (2 * [1]^2) * (1 + 2*x^2 -2*x*cos(y))))", 0, TMath::Infinity(), 0, TMath::Pi());
        // Heaviside
        //TF2 * cutoff_heaviside = new TF2("cutoff_heaviside", heaviside, 0, TMath::Infinity(), 0, TMath::Pi(), 2);
        TF1 * cutoff_heaviside = new TF1("cutoff_heaviside", "(x > 2) ? 0.0 : 1.0", 0, TMath::Infinity());
        //std::cerr << cutoff->Eval(3, 2.5);
        //std::cerr << cutoff->Integral(0, 77, 0, TMath::Pi());
        //cutoff->SetParameter(0, 1.0);
        //cutoff->SetParameter(1, 2.0);
        //cutoff->Draw("surf2");

        std::map<std::string, TF1 *> cutoffs = {
            {"gaussian", cutoff_gaussian},
            {"lorentzian1", cutoff_lorentzian1},
            {"lorentzian2", cutoff_lorentzian2},
            {"lorentzian3", cutoff_lorentzian3},
            {"rigid", cutoff_heaviside}
        };
        cutoffs[cutoff_type]->SetName("cutoff");

        if (rank > 0)
        {
            // Set up Filenames
            std::string s = encode_parameters(nb_events, rho, max_y, cutoff_type);
            const char * tree_file = s.c_str();
            std::string s2 = "lookup_table_" + cutoff_type + "_cutoff" + std::to_string(rho) + "_rank" + std::to_string(rank);
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
        delete cutoff_lorentzian1;
        delete cutoff_lorentzian2;
        delete cutoff_lorentzian3;
        delete cutoff_heaviside;

        MPI_Finalize();
    }
    else if (val == "fluctuations")
    {
        int rank, size; // rank of the process and number of processes
        MPI_Init(NULL, NULL);

        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        std::string filename;
        unsigned int len;

        if (rank == 0)
        {
            std::cerr << size << " processes" << std::endl;
            std::ifstream files(argv[2]);
            if (files.is_open())
            {
                int i = 0;
                while (files >> filename)
                {
                    ++i;
                    if (i >= size)
                    {
                        std::cout << "Warning : not enough processes. I stop reading filenames here." << std::endl;
                        break;
                    }
                    len = filename.length();
                    std::cout << "Sending to process " << i << std::endl;
                    MPI_Send(&len, 1, MPI_UNSIGNED, i, 0, MPI_COMM_WORLD);
                    MPI_Send(filename.c_str(), filename.length(), MPI_CHAR, i, 0, MPI_COMM_WORLD);
                }
                files.close();
            }            
        }
        else // TODO fix it if process rank too big
        {
            MPI_Recv(&len, 1, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            std::vector<char> filename_temp(len);
            MPI_Recv(filename_temp.data(), len, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            filename.assign(filename_temp.begin(), filename_temp.end());
            std::cout << "Process " << rank << " with " << filename << std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);

        if (rank > 0)
        {
            std::string filename_hist = filename + "hist";
            Double_t max_y, rho;
            std::string cutoff_type;
            int nb_events;
            decode_parameters(filename, &rho, &max_y, &cutoff_type, &nb_events);
            Double_t r = 0.05;
            fluctuations(max_y, 1.0, rho, r, filename.c_str(), filename_hist.c_str(), true);
        }

        MPI_Finalize();
    }
    else if (val == "check")
    {
        Double_t max_y, x01, rho;
        std::string cutoff_type;
        int nb_events;
        std::string filename = argv[2];
        TFile f(filename.c_str(), "READ");
        int nb_events_real = f.GetListOfKeys()->GetSize();

        decode_parameters(filename, &rho, &max_y, &cutoff_type, &nb_events);
        std::cout << "File " << filename << std::endl;
        std::cout << "Parameters:\n" << "\tRho = " << rho << "\n\tY_max = " << max_y << "\n\tCutoff type = " << cutoff_type;
        std::cout << "\nNumber of events: " << nb_events_real;
        if (nb_events != nb_events_real) std::cout << " -- WARNING should be " << nb_events;
        std::cout << std::endl;

        //f.Close(); // segmentation violation ??
    }
    else if (val == "fit-bare-r")
    {
        Double_t rho = std::stod(argv[2]);
        Double_t max_y = std::stod(argv[3]);
        fit_bare_r(rho, max_y, myapp);
    }
    else if (val == "ancestors")
    {
        int nb_events = std::stoi(argv[2]);
        Double_t rho = std::stoi(argv[3]);
        Double_t max_y = std::stoi(argv[4]);
        CommonAncestorPlot(myapp, nb_events, max_y, 1.0, rho);
    }
    else if (val == "draw-cutoffs")
    {
        draw_cutoffs(myapp);
    }
    // General fit
    //general_plot(myapp);
    //stat_events(myapp, max_y, 1.0);

    // Compute biggest children
    /*TFile f("tree.root");
    TTree * tree;
    f.GetObject("tree0", tree); 
    //tree->MakeClass("EventTree");
    compute_biggest_child(tree, myapp);*/

    //myapp->Run();
    //delete myapp;
    return EXIT_SUCCESS;
}