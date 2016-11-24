#include "utils.h"
#include "event.h"
#include "EventTree.h"

#include <TFile.h>
#include <TMath.h>
#include <TMinuit.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TH3D.h>
#include <TVector.h>
#include <TString.h>
#include <TRandom.h>
#include <TList.h>
#include <TGraph.h>
#include <TThread.h>
#include <TROOT.h>
#include <TEntryList.h>

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

void BinLogX(TH1 *h)
{
    TAxis *axis = h->GetXaxis();
    int bins = axis->GetNbins();

    Axis_t from = axis->GetXmin();
    Axis_t to = axis->GetXmax();
    Axis_t width = (to - from) / bins;
    Axis_t *new_bins = new Axis_t[bins + 1];

    for (int i = 0; i <= bins; i++) {
        new_bins[i] = TMath::Power(10, from + i * width);
        std::cerr << new_bins[i] << std::endl;
    }
    axis->Set(bins, new_bins);
    delete new_bins;
}

/* - To set ROOT environment variables add to your ~/.rootrc file for instance :
 * Hist.Binning.2D.x: 500
 * - The binning information for saving the result of a Draw to an histogram is
 * taken from the environment variables Hist.Binning.?D.? (see TTree.Draw doc)
 * - Note for the TH1.Rebin (quote from the doc) :
 * NOTE: The bin edges specified in xbins should correspond to bin edges in the
 * original histogram. If a bin edge in the new histogram is in the middle of a 
 * bin in the original histogram, all entries in the split bin in the original 
 * histogram will be transfered to the lower of the two possible bins in the new 
 * histogram. This is probably not what you want.
 */
Axis_t * BinLogX2(TH1 *h, int bins)
{
    TAxis *axis = h->GetXaxis();
    int nbins = bins;
    if (bins == -1) nbins = axis->GetNbins();

    std::cerr << axis->GetXmin() << " " << TMath::Log(axis->GetXmin()) << std::endl;
    Axis_t from = TMath::Max(TMath::Log(axis->GetXmin()), -3.0);
    Axis_t to = TMath::Log(axis->GetXmax());
    Axis_t width = (to - from) / nbins;
    Axis_t *new_bins = new Axis_t[nbins + 1];

    for (int i = 0; i <= nbins; i++) {
        new_bins[i] = TMath::Power(10, from + i * width);
        std::cerr << new_bins[i] << std::endl;
    }
    axis->Set(bins, new_bins);
    return new_bins;
    //return h->Rebin(bins, "hnew", new_bins);
}

Axis_t * BinLogX3(int nbins)
{
    Axis_t from = -5.0;
    Axis_t to = 0.0;
    Axis_t width = (to - from) / nbins;
    Axis_t *new_bins = new Axis_t[nbins + 1];

    for (int i = 0; i <= nbins; i++) {
        new_bins[i] = TMath::Power(10, from + i * width);
        std::cerr << new_bins[i] << std::endl;
    }
    //axis->Set(bins, new_bins);
    return new_bins;
    //return h->Rebin(bins, "hnew", new_bins);
}

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

/* Compute fluctuations p_n
* = probability to have n dipoles of size >= r from a dipole of size x01 until
* rapidity y_max
*/
void fluctuations(TApplication * myapp, int nb_events, Double_t max_y, Double_t x01, Double_t rho)
{
    Double_t r = 0.5;
    bool logX = true;

    TCanvas c;

    c.cd(1);
    if (logX) gPad->SetLogx();
    gPad->SetLogy();
    TH1F * hfluct = new TH1F("hfluct", TString::Format("#splitline{Fluctuations over %d events}{with y_max = %.12g, r = %.12g and rho = %.12g}", nb_events, max_y, r, rho), 100, 0, 100);
    hfluct->GetXaxis()->SetTitle("Number of events n(r, x01, y_max)");
    hfluct->GetYaxis()->SetTitle("p_n(r, x01, y)");
    if (logX) BinLogX2(hfluct, 100);

    TFile f("tree.root");

    for (int j = 0; j < nb_events; ++j)
    {
        TTree * tree;
        f.GetObject(TString::Format("tree%d", j), tree);
        //tree->Print();
        std::cerr << j << std::endl;
        int n = tree->Draw("radius", TString::Format("isLeaf && radius >= %.12g", r), "goff");
        hfluct->Fill(n);
    }
    if (logX) hfluct->Sumw2();
    if (logX) hfluct->Scale(1, "width");
    hfluct->Draw();

    TF1 * pn = new TF1("pn", "[0] / x * [1] * [1] / ([2] * [2]) * exp(- log(x)*log(x)/(4*[3]))", 3, 50);
    pn->FixParameter(1, x01);
    pn->FixParameter(2, r);
    pn->FixParameter(3, max_y);
    hfluct->Fit("pn", "IR");

    c.SetTitle(TString::Format("Chi2 : %.12g", pn->GetChisquare()));
    c.Update();

    myapp->Run();
}


/* Compute and fit average number of dipoles with size >= r, x01 and max_y
 * (fit with Bessel function)
 */
void stat_events(TApplication * myapp, int nb_events, Double_t max_y, Double_t x01)
{
    TCanvas c;
    //c.Divide(2,1,0.05,0.05);
    c.cd(1);
    gPad->SetLogy();

    //TH1F * hf0 = new TH1F("hf0", "Total", 200, 0, 1);
    //Axis_t * new_bins = BinLogX3(100);
    //TH1F * hf = new TH1F("hf", "Total 2", 100, new_bins);
    //hf->Rebin(200, "hf", new_bins);

    TFile f("tree.root");

    for (int j = 0; j < nb_events; ++j)
    {
        TTree * tree;
        f.GetObject(TString::Format("tree%d", j), tree);
        tree->Print();
        tree->Draw("radius>>+hf", "isLeaf", "goff");
        //TH1F *htemp = (TH1F*)gPad->GetPrimitive("hf");
        tree->GetHistogram()->Print();
        //std::cerr << tree->GetHistogram()->GetXaxis()->GetNbins() << std::endl;
        //tree->GetHistogram()->SetDirectory(0);
        //std::cerr << hf->GetEntries() << std::endl;
        //std::cerr << htemp->GetMean() << std::endl;
        //number = htemp->GetEntries();
    }

    TH1F * htemp = (TH1F*)gDirectory->Get("hf");
    htemp->Print();
    //htemp->SetBit(TH1::klsAverage);
    std::cerr << htemp->GetEntries() << std::endl;
    // TH1::AddDirectory(kFALSE);   sets a global switch disabling the reference
    htemp->SetDirectory(0);
    TH1F * htemp_cum = (TH1F*) htemp->GetCumulative(false);
    htemp_cum->Draw(); // E1 for error bars
    htemp_cum->SetTitle(TString::Format("#splitline{Average number of dipoles of radius >= r over %d events}{with x01=%.12g and y_max=%.12g}", nb_events, 1.0, max_y));

    // First limit : rho >> y_max
    Double_t max_r = x01 / 10. * TMath::Exp(max_y / 2.);
    std::cerr << max_r << std::endl;
    TF1 * nf = new TF1("nf", "[2] * TMath::BesselI0(2. * sqrt([0] * log([1]*[1] / (x*x))))", 0.01, max_r);
    nf->FixParameter(0, max_y);
    nf->FixParameter(1, x01);

    htemp_cum->Fit("nf", "R");

    Double_t min_r = x01 / 2. * TMath::Exp(max_y / 2.);
    std::cerr << min_r << std::endl;
    TF1 * nf2 = new TF1("nf2", "[0] * [1] * [1] / (x*x) * exp(2. * sqrt([2] * log([1] * [1] / (x * x))))", 1.0, 1.5);
    nf2->FixParameter(1, x01);
    nf2->FixParameter(2, max_y);
    //htemp_cum->Fit("nf2", "R+");

    // Second limit : rho 
    // Log binning on x axis : not accurate because new_bins doesn't overlap old bins
    // See comment over BinLogX2
    /*c.cd(2);
    gPad->SetLogy();
    gPad->SetLogx();
    Axis_t * new_bins = BinLogX2(htemp_cum, 100);
    TH1F * htemp_cum_rebin = (TH1F*) htemp_cum->Rebin(100, "hnew", new_bins);
    htemp_cum_rebin->Draw();*/

    c.SetTitle(TString::Format("Chi2 : %.12g", nf->GetChisquare()));
    c.Update();

    myapp->Run();
}

// Interesting view !
void draw_tree(TApplication * myapp, TTree * tree)
{
    tree->Draw("rapidity:coord.X()", "coord.X() > -0.1 && coord.X() < 0.1");
    myapp->Run();
}

// Generate *nb_events* events with same parameters rho and max_y
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

// Find common ancestor of two dipoles given their index in TTree
Long64_t GetCommonAncestors(TTree * tree, Long64_t i1, Long64_t i2)
{
    Long64_t nentries = tree->GetEntries();
    EventTree * t = new EventTree(tree);

    // Get depths
    t->GetEntry(i1);
    Long64_t depth1 = t->depth;
    t->GetEntry(i2);
    Long64_t depth2 = t->depth;

    // Assumption : depth(i3) >= depth(i4)
    Long64_t i3 = i2; // if depth1 <= depth2
    Long64_t i4 = i1;
    if (depth1 > depth2)
    {
        i3 = i1;
        i4 = i2;
    }

    // Put i3 and i4 at same depth
    if (depth1 != depth2)
    {
        for (int k = 0; k < depth2 - depth1; ++k)
        {
            t->GetEntry(i3);
            i3 = t->index_parent;
        }
    }

    Long64_t p1 = i4;
    Long64_t p2 = i3;
    while(p1 != p2)
    {
        //std::cerr << p1 << " " << p2 << std::endl;
        t->GetEntry(p1);
        p1 = t->index_parent;
        t->GetEntry(p2);
        p2 = t->index_parent;
    }
    return p1;    
}

/* Select randomly k leaves
* Return true if it was possible
*/
bool RandomSelectkLeaves(TTree * tree, Long64_t indexes[], int k)
{
    //EventTree * t = new EventTree(tree);
    Long64_t nleaves = tree->GetEntries("isLeaf");
    //Long64_t nentries = tree->GetEntries();
    //std::cerr << nleaves << " leaves" << std::endl;
    if (nleaves < k) return false;

    //tree->Scan("depth:index_children:index_parent");
    tree->Draw(">>leaves", "isLeaf", "entrylist");
    TEntryList * elist = (TEntryList *) gDirectory->Get("leaves");
    for (Long64_t j = 0; j < k; ++j)
    {
        bool new_leaf = false;
        Long64_t i;
        while(!new_leaf)
        {
            new_leaf = true;
            i = gRandom->Integer(nleaves);
            // Check whether we already picked this index
            for (int l = 0; l < j; ++l)
            {
                if (indexes[l] == i) new_leaf = false;
            }
        }
        indexes[j] = elist->GetEntry(i);
    }
    return true;
}

void CommonAncestorPlot(TApplication * myapp, int nb_events, Double_t max_y, Double_t x01, Double_t rho)
{
    TCanvas c;
    c.cd(1);
    gPad->SetLogy();

    TH1F * hist = new TH1F("hancestor", TString::Format("#splitline{Common ancestor over %d events}{with y_max = %.12g and rho = %.12g}", nb_events, max_y, rho), 100, 0, max_y);
    hist->GetXaxis()->SetTitle("Rapidity");
    hist->GetYaxis()->SetTitle("Number of common ancestors");

    TFile f("tree.root");
    Long64_t indexes[2];

    for (int j = 0; j < nb_events; ++j)
    {
        TTree * tree;
        f.GetObject(TString::Format("tree%d", j), tree);
        //tree->Print();
        std::cerr << j << std::endl;
        if (RandomSelectkLeaves(tree, indexes, 2))
        {
            //std::cerr << indexes[0] << " " << indexes[1] << std::endl;
            Long64_t a = GetCommonAncestors(tree, indexes[0], indexes[1]);
            if (a < 0) a = 0;
            //std::cerr << a << std::endl;

            // Get rapidity of common ancestor
            EventTree * t = new EventTree(tree);
            t->GetEntry(a);
            hist->Fill(t->rapidity);
        }

        //tree->GetHistogram()->Print();
        //std::cerr << tree->GetHistogram()->GetXaxis()->GetNbins() << std::endl;
        //tree->GetHistogram()->SetDirectory(0);
        //std::cerr << hf->GetEntries() << std::endl;
        //std::cerr << htemp->GetMean() << std::endl;
        //number = htemp->GetEntries();
    }

    hist->Draw();

    myapp->Run();
}

void compute_biggest_child(TTree * tree, TApplication * myapp)
{
    TFile * f = new TFile("treeFriend.root", "recreate");
    TTree * treeFriend = tree->CopyTree("");
    treeFriend->SetName("treeFriend");
    //treeFriend->BuildIndex("radius", "isLeaf");

    Double_t max_children_radius, max_children_radius2;
    EventTree * t = new EventTree(treeFriend);

    TBranch * biggestChildSize = treeFriend->Branch("max_children_radius_inv", &max_children_radius, "max_children_radius_inv/D");
    TBranch * biggestChildSize2 = treeFriend->Branch("max_children_radius", &max_children_radius2, "max_children_radius/D");
    
    Long64_t nentries = tree->GetEntries(); // read the number of entries in tree

    std::cerr << nentries << std::endl;
    for (Long64_t i = nentries-1; i >= 0; --i)
    {
        t->GetEntry(i);
        if (t->isLeaf)
        {
            max_children_radius = t->radius;
        }
        else
        {
            max_children_radius = t->radius;
            Long64_t index_children = t->index_children;
            // FIXME binary tree is not complete => formula not exact
            t->GetEntry(index_children);
            max_children_radius = TMath::Max(max_children_radius, t->radius);
            t->GetEntry(index_children+1);
            max_children_radius = TMath::Max(max_children_radius, t->radius);
            if (max_children_radius < 0.)
            {
                std::cerr << i << " Negative" << std::endl;
            }
        }
        //std::cerr << i << " " << max_children_radius << std::endl;
        biggestChildSize->Fill();
    }
    for (Long64_t i = 0; i < nentries; ++i)
    {
        biggestChildSize->GetEntry(nentries-1-i);
        max_children_radius2 = max_children_radius;
        biggestChildSize2->Fill();
    }
    //treeFriend->GetListOfBranches()->Remove(biggestChildSize);

    treeFriend->Write();
    treeFriend->Print();
    //treeFriend->Scan("max_children_radius_inv:max_children_radius");

    tree->AddFriend(treeFriend);
    //tree->StartViewer();
    TCanvas C;
    gPad->SetLogz();
    //tree->Draw("radius:treeFriend.max_children_radius", "radius < 0.2 && treeFriend.max_children_radius < 0.2");
    tree->Draw("radius:treeFriend.max_children_radius>>hmax", "radius<0.1 && treeFriend.max_children_radius < 0.1", "colz");
    TH3D * htemp = (TH3D*) gDirectory->Get("hmax");
    std::cerr << htemp->GetEntries() << std::endl;
    htemp->SetDirectory(0); // Don't forget this one !!
    //htemp->Draw("box");
    //htemp->DrawPanel();
    //htemp->FitPanel();
    htemp->Print();

    //C.Draw();
    //C.Update();
    // treeFriend->Write("", TObject::kOverwrite());
    delete f;
    myapp->Run();
}