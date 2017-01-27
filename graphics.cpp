#include "graphics.h"
#include "event.h"

#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLine.h>
#include <TH1D.h>

void draw_cutoffs(TApplication * myapp)
{
    Double_t r_max = 4;
    // Gaussian cutoff
    TF2 * cutoff = new TF2("cutoff", "exp(-[0] / (2 * [1]^2) * (1 + 2*x^2 -2*x*cos(y)))", 0, r_max, 0, TMath::Pi());
    cutoff->SetParameter(0, 1.0);
    cutoff->SetParameter(1, 2.0);
    // Lorentzian cutoff
    TF2 * cutoff_lorentzian = new TF2("cutoff_lorentzian", "1 / (1 + ([0]^2 / (2 * [1]^2) * (1 + 2*x^2 -2*x*cos(y))))", 0, r_max, 0, TMath::Pi());
    cutoff_lorentzian->SetParameter(0, 1.0);
    cutoff_lorentzian->SetParameter(1, 2.0);
    cutoff_lorentzian->SetLineColor(kBlue);
    // Maxwellian cutoff
    TF2 * cutoff_maxwellian = new TF2("cutoff_maxwellian", "1 / (1 + exp([0]^2 / (2 * [1]^2) * (1 + 2*x^2 -2*x*cos(y))))", 0, r_max, 0, TMath::Pi());
    cutoff_maxwellian->SetParameter(0, 1.0);
    cutoff_maxwellian->SetParameter(1, 2.0);
    cutoff_maxwellian->SetLineColor(kBlack);
    // Heaviside
    TF2 * cutoff_rigid = new TF2("cutoff_rigid", "( (-1+x^2 -x*cos(y))> 0 ? 0.0 : 1.0) * y/y", 0, r_max, 0, TMath::Pi());
    //cutoff->SetParameter(0, 1.0);
    //cutoff->SetParameter(1, 2.0);

    cutoff->Draw(); 
    cutoff_lorentzian->Draw("same"); 

    myapp->Run();  
}

void fit_bare_r(Double_t rho, Double_t max_y, TApplication * myapp)
{
    // Event 1 : without cutoff
    // Event 2 : with cutoff, with raw computation
    // Event 3 : with cutoff, no raw computation
    Event e1(rho, max_y, "lookup_table");
    // Gaussian cutoff
    //TF2 * cutoff = new TF2("cutoff", "exp(-[0] / (2 * [1]^2) * (1 + 2*x^2 -2*x*cos(y)))", 0, TMath::Infinity(), 0, TMath::Pi());
    // Lorentzian cutoff
    //TF2 * cutoff = new TF2("cutoff", "1 / (1 + ([0]^2 / (2 * [1]^2) * (1 + 2*x^2 -2*x*cos(y))))", 0, TMath::Infinity(), 0, TMath::Pi());
    // Maxwellian cutoff
    //TF2 * cutoff = new TF2("cutoff", "1 / (1 + exp([0]^2 / (2 * [1]^2) * (1 + 2*x^2 -2*x*cos(y))))", 0, TMath::Infinity(), 0, TMath::Pi());
    // Heaviside
    TF2 * cutoff = new TF2("cutoff", "(x > (2*[0]) ? 0.0 : 1.0) * y/y", 0, TMath::Infinity(), 0, TMath::Pi());
    Event e2(rho, max_y, "lookup_table", cutoff, true, true);
    Event e3(rho, max_y, "lookup_table", cutoff, true, false);

    TCanvas canvas("canvas", "Sizes", 1080, 780);
    gPad->SetLogy();
    gStyle->SetPadTickY(1); // y ticks on the right too
    gStyle->SetOptStat(0); // get rid of the statistics box
    
    Double_t max_r = 4;
    int N = 100000;
    TH1D hist1("hist1", "hist1", 100, 0, max_r);
    TH1D hist2("hist2", "hist2", 100, 0, max_r);
    TH1D hist3("hist3", "hist3", 100, 0, max_r);
    for (int i = 0; i <N; ++i)
    {
        hist1.Fill(e1.r_generate());
        hist2.Fill(e2.r_generate());
        hist3.Fill(e3.r_generate2());
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
    myapp->Run();
}

void fit_fluctuations(Double_t rho, Double_t max_y, TApplication * myapp)
{
    // Compare 3 methods
    Event e1(rho, max_y, "lookup_table");
    TF2 * cutoff = new TF2("cutoff", "exp(-[0] / (2 * [1]^2) * (1 + 2*x^2 -2*x*cos(y)))", 0, TMath::Infinity(), 0, TMath::Pi());
    Event e2(rho, max_y, "lookup_table", cutoff, true, true);
    Event e3(rho, max_y, "lookup_table", cutoff, true, false);

    TCanvas canvas("canvas", "Sizes", 1080, 780);
    gPad->SetLogy();
    gStyle->SetPadTickY(1); // y ticks on the right too
    gStyle->SetOptStat(0); // get rid of the statistics box
    
    Double_t max_r = 1000;
    Double_t r = 0.1;
    int N = 10;
    TH1D hist1("hist1", "hist1", 100, 0, max_r);
    TH1D hist2("hist2", "hist2", 100, 0, max_r);
    TH1D hist3("hist3", "hist3", 100, 0, max_r);
    for (int i = 0; i <N; ++i)
    {
        TTree * tree1 = new TTree(TString::Format("tree1%d", i), "Dipole splitting");
        e1.make_tree(tree1, false);
        int n = tree1->Draw("radius", TString::Format("isLeaf && radius >= %.12g", r), "goff");
        hist1.Fill(n);
        std::cerr << n << std::endl;
    }
    for (int i = 0; i <N; ++i)
    {
        TTree * tree2 = new TTree(TString::Format("tree2%d", i), "Dipole splitting");
        e2.make_tree(tree2, false);
        int n = tree2->Draw("radius", TString::Format("isLeaf && radius >= %.12g", r), "goff");
        hist2.Fill(n);
        std::cerr << n << std::endl;
    }    
    for (int i = 0; i <N; ++i)
    {
        TTree * tree3 = new TTree(TString::Format("tree3%d", i), "Dipole splitting");
        e3.make_tree(tree3, false);
        int n = tree3->Draw("radius", TString::Format("isLeaf && radius >= %.12g", r), "goff");
        hist3.Fill(n);
        std::cerr << n << std::endl;
    }

    hist1.SetLineWidth(2);
    hist1.SetMarkerStyle(8);
    hist1.SetMarkerSize(0.7);
    hist1.Draw("E1");
    TF1 fr1("fr1", "[0] / (x*(1. - x*x))", 0.02, 0.5);
    TF1 fr2("fr2", "2. * [0] / (x*TMath::Abs(1. - x*x)) * TMath::ATan(TMath::Abs(1. - x)/(1. + x) * TMath::Sqrt((x+1./2.)/(x-1./2.)))", 0.5, 2);
    fr1.SetLineColor(kGreen-6);
    fr2.SetLineColor(kGreen-6);
    //hist1.Fit("fr1", "R");
    //hist1.Fit("fr2", "R+");

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
    //legend.AddEntry(&fr1, "Fit (formula without cutoff)", "L");
    legend.Draw();

    TLine line(rho, 0, rho, hist1.GetMaximum());
    line.SetLineColor(kBlack);
    line.SetLineStyle(7);
    line.Draw();

    //TString::Format("#splitline{Distribution de la taille r du dipole (#rho = %.12g)}{%d bins - chi2 = %.12g et %.12g}", rho, hist.GetSize()-2, chi2_1, chi2_2)
    hist1.SetTitle(TString::Format("%d events with x01 = %.12g, rho = %.12g, y_max = %.12g;Number n of dipoles (leaves) of size over %.12g;Fluctuations p_n", N, 1.0, rho, max_y, r));
    canvas.Update();
    myapp->Run();
}