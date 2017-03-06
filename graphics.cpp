#include "graphics.h"
#include "event.h"
#include "utils.h"

#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLine.h>
#include <TH1D.h>
#include <TFile.h>
#include <TImage.h>
#include <THStack.h>
#include <TGraphErrors.h>

#include <fstream>
#include <regex>

void draw_cutoffs(TApplication * myapp, std::map<std::string, TF1 *> cutoffs, 
                    std::map<std::string, int> colors)
{
    Double_t r_max = 4;
    int i = 0;
    Double_t xmin = 0.0;
    Double_t xmax = 5.0;

    TLegend legend(0.6, 0.45, 0.8, 0.8);
    //legend.SetEntrySeparation(3.0);
    legend.SetFillColor(0); // white bg
    legend.SetBorderSize(0); // get rid of the box
    legend.SetTextSize(0.045);
    legend.SetHeader("Infrared Cutoffs");

    for (auto const &cutoff : cutoffs)
    {
        std::cout << cutoff.first << std::endl;
        cutoff.second->SetParameter(0, 1.0);
        cutoff.second->SetParameter(1, 2.0);
        cutoff.second->SetRange(xmin, xmax);
        std::string s = cutoff.first + "12";

        if (cutoff.second->GetNdim() > 1)
        {
            TF12 * cutoff12 = new TF12(s.c_str(), (TF2*) cutoff.second, 0.0, "x");
            cutoff12->SetLineColor(colors[cutoff.first]);
            cutoff12->SetLineWidth(1);
            //cutoff12->GetHistogram()->GetXaxis()->SetTitle("Size of dipole r");
            cutoff12->SetLineStyle(1);
            cutoff12->SetTitle("");
            cutoff12->DrawCopy((i == 0) ? "" : "same");

            legend.AddEntry(cutoff12, cutoff.first.c_str(), "L"); 

            /*cutoff12->SetLineWidth(1);
            cutoff12->SetXY(TMath::Pi()/3.0);
            cutoff12->SetLineStyle(2);
            cutoff12->DrawCopy("same");
            cutoff12->SetXY(TMath::Pi()*2.0/3.0);
            cutoff12->SetLineStyle(8);
            cutoff12->DrawCopy("same");*/

            cutoff12->SetXY(TMath::Pi());
            cutoff12->SetLineStyle(9);
            cutoff12->SetLineWidth(1);
            cutoff12->DrawCopy("same"); 

            //cutoff12->GetXaxis()->SetTitle("r");
            //gPad->Modified();
            //gPad->Update();
        }
        else
        {
            cutoff.second->SetTitle("");
            cutoff.second->SetLineColor(colors[cutoff.first]);
            cutoff.second->SetLineWidth(1);
            cutoff.second->SetLineStyle(1);
            cutoff.second->Draw((i == 0) ? "" : "same");
            legend.AddEntry(cutoff.second, cutoff.first.c_str(), "L"); 
            //cutoff.second->GetXaxis()->SetTitle("r");
        }

        ++i;         
    }

    legend.Draw();

    gPad->Update();

    gPad->Print("cutoffs.eps");
    /*TImage * img = TImage::Create();
    img->FromPad(gPad->GetCanvas());
    img->WriteImage("cutoffs.png");
    delete img;*/

    myapp->Run();  
}

void fit_bare_r(Double_t rho, Double_t max_y, TApplication * myapp)
{
    // Event 1 : without cutoff
    // Event 2 : with cutoff, with raw computation
    // Event 3 : with cutoff, no raw computation
    Event e1(rho, max_y, 2.0, "lookup_table");
    // Gaussian cutoff
    //TF2 * cutoff = new TF2("cutoff", "exp(-[0] / (2 * [1]^2) * (1 + 2*x^2 -2*x*cos(y)))", 0, TMath::Infinity(), 0, TMath::Pi());
    // Lorentzian cutoff
    //TF2 * cutoff = new TF2("cutoff", "1 / (1 + ([0]^2 / (2 * [1]^2) * (1 + 2*x^2 -2*x*cos(y))))", 0, TMath::Infinity(), 0, TMath::Pi());
    // Maxwellian cutoff
    //TF2 * cutoff = new TF2("cutoff", "1 / (1 + exp([0]^2 / (2 * [1]^2) * (1 + 2*x^2 -2*x*cos(y))))", 0, TMath::Infinity(), 0, TMath::Pi());
    // Heaviside
    //TF2 * cutoff = new TF2("cutoff", "(x > (2*[0]) ? 0.0 : 1.0) * y/y", 0, TMath::Infinity(), 0, TMath::Pi());
    TF1 * cutoff = new TF1("cutoff", "(x > 2) ? 0.0 : 1.0", 0, TMath::Infinity());
    Event e2(rho, max_y, 2.0, "lookup_table", cutoff, true, true);
    Event e3(rho, max_y, 2.0, "lookup_table", cutoff, true, false);

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
    Event e1(rho, max_y, 2.0, "lookup_table");
    TF2 * cutoff = new TF2("cutoff", "exp(-[0] / (2 * [1]^2) * (1 + 2*x^2 -2*x*cos(y)))", 0, TMath::Infinity(), 0, TMath::Pi());
    Event e2(rho, max_y, 2.0, "lookup_table", cutoff, true, true);
    Event e3(rho, max_y, 2.0, "lookup_table", cutoff, true, false);

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

void compare_histo(TApplication * myapp, std::string histofiles)
{
    TCanvas canvas("canvas", "Sizes", 1080, 780);
    gPad->SetLogy();
    gStyle->SetOptStat(0); // get rid of the statistics box

    std::ifstream files(histofiles);
    int lineColor, lineWidth;
    std::string title, filename;

    TLegend legend(0.2, 0.2, 0.4, 0.4);
    legend.SetFillColor(0); // white bg
    legend.SetBorderSize(0); // get rid of the box
    legend.SetTextSize(0.045);

    THStack * stack = new THStack("hs", "");

    if (files.is_open())
    {
        int i = 0;
        //files >> title;
        //legend.SetHeader(title.c_str());
        while (files >> lineColor >> lineWidth >> title >> filename)
        {
            TFile f(filename.c_str(), "READ");
            TH1F * h = n_to_nbar((TH1F*) f.Get("hfluct"));
            h->SetLineColor(lineColor);
            h->SetLineWidth(lineWidth);
            h->SetDirectory(0);
            //h->Draw((i == 0) ? "E1" : "E1 same");
            legend.AddEntry(h, title.c_str(), "L");
            std::cout << filename << std::endl;

            std::string cutoff_name = "pn_cutoff" + std::to_string(i);
            TF1 pn_cutoff(cutoff_name.c_str(), "[0] * [1]^2 / [2]^2 * exp(-[1]^2/(2. * [2]^2)) * exp(-x/[3])", 1, 4);
            pn_cutoff.FixParameter(1, 1.0);
            pn_cutoff.FixParameter(2, 2.0); // R
            pn_cutoff.SetParLimits(3, 0.01, 20);
            h->Fit(cutoff_name.c_str(), "*IR0+", "goff");

            stack->Add(h);
            ++i;
        }
        files.close();
    } 
    stack->Draw("nostack");
    stack->GetXaxis()->SetTitle("Multiplicity n/#bar{n}");
    stack->GetYaxis()->SetTitle("Number of events with multiplicity n");
    //legend.AddEntry(&fr1, "Fit (formula without cutoff)", "L");
    legend.Draw();

    gPad->Print("compare_histo.eps");
    // Parameter 0 : proportionality | 3 : c 
    //canvas.SetTitle(TString::Format("Chi2 : %.12g", pn_cutoff.GetChisquare()));

    //canvas.Update();
    myapp->Run();
}

void compare_c(TApplication * myapp, std::string histofiles)
{
    TCanvas canvas("canvas", "Sizes", 1080, 780);
    gPad->SetLogy();
    gStyle->SetOptStat(0); // get rid of the statistics box

    std::ifstream files(histofiles);   
    std::string filename;

    std::string double_regex = "[-+]?[0-9]*\.?[0-9]+";
    std::regex r("_r(" + double_regex + ")_FINAL");
    std::smatch m;

    std::vector<Double_t> rvalues, c, ex, ey;

    if (files.is_open())
    {
        int i = 0;
        //files >> title;
        //legend.SetHeader(title.c_str());
        while (files >> filename)
        {
            if (std::regex_search(filename, m, r))
            {           
                Double_t r = std::stod(m[1]);
                TFile f(filename.c_str(), "READ");
                TH1F * h = n_to_nbar((TH1F*) f.Get("hfluct"));
                h->SetDirectory(0);
                //h->Draw((i == 0) ? "E1" : "E1 same");
                //legend.AddEntry(h, title.c_str(), "L");
                std::cout << filename << std::endl;

                std::string cutoff_name = "pn_cutoff" + std::to_string(i);
                TF1 pn_cutoff(cutoff_name.c_str(), "[0] * [1]^2 / [2]^2 * exp(-[1]^2/(2. * [2]^2)) * exp(-x/[3])", 1, 4);
                pn_cutoff.FixParameter(1, 1.0);
                pn_cutoff.FixParameter(2, 2.0); // R
                pn_cutoff.SetParLimits(3, 0.01, 20);
                h->Fit(cutoff_name.c_str(), "*IR0Q+", "goff");
                std::cerr << r << " " << pn_cutoff.GetParameter(3) << " " << pn_cutoff.GetParError(3) << std::endl;
                rvalues.push_back(r);
                c.push_back(pn_cutoff.GetParameter(3));
                ex.push_back(0.0);
                ey.push_back(pn_cutoff.GetParError(3));
                ++i;
            }
            else
            {
                std::cerr << "Failed to decode r value." << std::endl;
            }
        }
        files.close();
    }    
    TGraphErrors * g = new TGraphErrors(c.size(), &(rvalues[0]), &(c[0]), &(ex[0]), &(ey[0]));
    g->Draw("AP");  
    myapp->Run();
}