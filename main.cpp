#include "main.h"
#include "event.h"

#include <TApplication.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TMinuit.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <sstream>
#include <TVector.h>
#include <TString.h>
#include <TRandom.h>
#include <TList.h>
#include <TGraph.h>

/*
Double_t myFunction(Double_t par)
{
    return par /
}

void minuitFunction(int& nDim, double* gout, double& result, double par[], int flg)
{
    result = myFunction(par[0]);
}
*/
/* Nombre moyen de dipôles de taille r partant d'un dipôle 
 * de taille x01 en allant jusqu'à la rapidité y
 */
Double_t n(Double_t r, Double_t x01, Double_t y)
{
    // FIXME manque alpha barre : paramètre à ajuster ??
    return TMath::BesselI0(2 * TMath::Sqrt(y * TMath::Log(x01*x01 / (r*r))));
}

void general_plot(TApplication * myapp)
{
    TCanvas c;
    c.Divide(2,2,0.05,0.05);

    TFile f("tree.root");
    TTree * tree;
    f.GetObject("T", tree);
    tree->Print();
    //TF1 distrib("distrib", "sqrt([0])", 0, 10);

    TVector * v = (TVector *) tree->GetUserInfo()->At(0);
    Long64_t max_depth = v[0][0];
    std::cerr << "Max depth : " << max_depth << std::endl;
    tree->GetUserInfo()->Print();

    c.cd(1);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetTickx();
    //tree->Fit("distrib", "radius");
    //tree->Draw("radius", TString::Format("depth == %lld", max_depth));
    tree->Draw("radius", "isLeaf");
    TH1F *htemp = (TH1F*)gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetTitle("Radius r");
    htemp->GetYaxis()->SetTitle("Number of dipoles of radius r");
    
    c.cd(2);
    gPad->SetLogy();
    tree->Draw("rapidity", "isLeaf && radius < 1.5 && radius > 0.5");
    TH1F *htemp2 = (TH1F*)gPad->GetPrimitive("htemp");
    TF1 * rapidity = new TF1("rapidity", "[1] * exp([0]*x)", 1, 1.4);
    htemp2->Fit("rapidity", "R");

    c.cd(3);
    gPad->SetLogy();
    //tree->Draw("coord.X():coord.Y()", "isLeaf && coord.X() > -1. && coord.X() < 1. && coord.Y() < 1. && coord.Y() > -1.");
    tree->Draw("radius", "isLeaf && radius < 0.2");
    TH1F *htemp3 = (TH1F*)gPad->GetPrimitive("htemp");
    TF1 * n = new TF1("n", "[0] * exp([1] * sqrt(- log(x)))", 0, 0.2);
    htemp3->Fit("n", "R");
    //tree->Draw("radius:rapidity", "isLeaf && radius < 0.5", "surf2");

    c.cd(4);
    /*tree->Draw("phi", "isLeaf");
    TH1F *htemp4 = (TH1F*)gPad->GetPrimitive("htemp");
    TF1 * pi = new TF1("pi", "[0]", 0, 6);
    htemp4->Fit("pi", "R");    */
    gPad->SetLogy();
    gPad->SetGridy();    
    tree->Draw("coord.X()", "isLeaf && coord.X() < 2. && coord.X() > -2.");
    /*tree->Draw("radius", "isLeaf");
    TH1F *htemp4 = (TH1F*)gPad->GetPrimitive("htemp");
    TF1 * n2 = new TF1("n2", "TMath::BesselI0([0] * sqrt(log(x)))", 0, 1);
    htemp4->Fit("n2", "R"); */   

    //tree->MakeSelector("DrawEvent");
    myapp->Run();
}

void stat_events(TApplication * myapp, int nb_events, Double_t max_y)
{
    TCanvas c;
    c.cd(1);
    gPad->SetLogy();
    TH1F * hf = new TH1F("hf", "Total", 100, 0, 1);
    TFile f("tree.root");
    int N = 500;
    Double_t r[N], nb[N];
    Double_t radius = 0.02;
    Double_t pas = 0.001;
    for (int i = 0; i < N; ++i)
    {
        radius += pas;
        int number = 0;
        for (int j = 0; j < nb_events; ++j)
        {
            //TFile f("tree_0.root");
            TTree * tree;
            f.GetObject(TString::Format("tree%d", j), tree);
            //tree->Print();
            tree->Draw(TString::Format("radius>>+h%d", i), TString::Format("isLeaf && radius > %.12g", radius));
            TH1F *htemp = (TH1F*)gPad->GetPrimitive(TString::Format("h%d", i));
            //std::cerr << htemp->GetEntries() << std::endl;
            //std::cerr << htemp->GetMean() << std::endl;
            number = htemp->GetEntries();
        }
        nb[i] = number;
        r[i] = radius;
        //std::cerr << radius << " " << number << std::endl;
    }

    //TH1F * htemp = (TH1F*)gPad->GetPrimitive("h");
    //std::cerr << htemp->GetEntries() << std::endl;
    //h->SetDirectory(0);
    //h->DrawCopy();
    //TH1F *htemp3 = (TH1F*)gPad->GetPrimitive("h");
    TF1 * nf = new TF1("nf", "[2] * TMath::BesselI0(2 * sqrt([0] * log([1]*[1] / (x*x))))", 0, 1);
    nf->FixParameter(0, max_y);
    nf->FixParameter(1, 1.0); // x01

    c.cd(2);
    TGraph * n = new TGraph(N, r, nb);
    n->SetTitle(TString::Format("#splitline{Average number of dipoles starting with a size x01=%.12g}{until rapidity %.12g with a size over r}", 1.0, max_y));
    n->Draw("A*");
        //n->GetXaxis()->SetRangeUser(0.,0.12);
        //n->GetYaxis()->SetRangeUser(-0.8, 0.8);    
    n->Draw("A*");
    c.Update();
    n->Fit("nf");

    myapp->Run();
}

void draw_tree(TApplication * myapp, TTree * tree)
{
    //gPad->GetCanvas()->Divide(1,2);
    //gPad->GetCanvas()->cd(2);
    tree->Draw("rapidity:coord.X()", "coord.X() > -0.1 && coord.X() < 0.1");
    myapp->Run();
}

void generate_events(int nb_events, Double_t rho, Double_t max_y)
{
    TFile * output = new TFile("tree.root", "recreate");
    std::cerr << "Rho = " << rho << " ; Maximum rapidity = " << max_y << std::endl;
    for (int j = 0; j < nb_events; ++j)
    {
        Event * e = new Event(rho, max_y);
        TTree * tree = e->make_tree("tree.root", TString::Format("tree%d", j), false);
    }        
    output->Write();
    std::cerr << "Generated and saved " << nb_events << " events in tree.root." << std::endl;
}

int main( int argc, const char* argv[] )
{
    TApplication * myapp = new TApplication("myapp", 0, 0);
    
    //connect(myapp,SIGNAL(lastWindowClosed()),TQtRootSlot::CintSlot(),SLOT(TerminateAndQuit());
    //p3a();

    std::istringstream iss( argv[1] );
    int val = 0; // By default only read file
    int nb_events = 1; // Number of events to read/generate
    Double_t rho = 0.01;
    Double_t max_y = 1.6;

    if (!(iss >> val))
    {
        std::cerr << "Not valid argument." << std::endl;
    }

    if (val)
    {
        generate_events(nb_events);
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

    //general_plot(myapp);
    //stat_events(myapp, nb_events, max_y);
    TFile f("tree.root");
    TTree * tree;
    f.GetObject("tree0", tree);    
    tree->Draw("coord.Y():coord.X()", "isLeaf");
    //draw_tree(myapp, tree);
    myapp->Run();

    return 0;
}