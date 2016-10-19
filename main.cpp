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
    tree->Draw("rapidity", "isLeaf && radius < 1.5 && radius > 0.5");
    TH1F *htemp2 = (TH1F*)gPad->GetPrimitive("htemp");
    TF1 * rapidity = new TF1("rapidity", "[1] * exp([0]*x)", 1, 1.4);
    htemp2->Fit("rapidity", "R");

    c.cd(3);
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

    //tree->MakeClass("MyClass");
    myapp->Run();
}

int main( int argc, const char* argv[] )
{
    TApplication * myapp = new TApplication("myapp", 0, 0);
    //connect(myapp,SIGNAL(lastWindowClosed()),TQtRootSlot::CintSlot(),SLOT(TerminateAndQuit());
    //p3a();

    std::istringstream iss( argv[1] );
    int val = 0; // By default only read file
    if (!(iss >> val))
    {
        std::cerr << "Not valid argument." << std::endl;
    }

    if (val)
    {
        Event * e = new Event();
        Long64_t max_depth = e->make_tree("tree.root");
    }

    //draw_tree(tree);
    //tree->Draw("radius");

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
    Event e;
    e.bare_distribution();
    myapp->Run();

    return 0;
}