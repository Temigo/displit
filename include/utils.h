#include "color_codes.h"

#include <TApplication.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF2.h>

void BinLogX(TH1 *h);
Axis_t * BinLogX2(TH1 *h, int bins);
Axis_t * BinLogX3(int nbins);
Double_t n(Double_t r, Double_t x01, Double_t y);

void general_plot(TApplication * myapp);
void fluctuations(TApplication * myapp, int nb_events, Double_t max_y, Double_t x01, Double_t rho);
void stat_events(TApplication * myapp, int nb_events, Double_t max_y, Double_t x01);
void draw_tree(TApplication * myapp, TTree * tree);
void generate_events(int nb_events, Double_t rho, Double_t max_y, bool with_cutoff = false, TF2 * cutoff = NULL);
Long64_t GetCommonAncestors(TTree * tree, Long64_t i1, Long64_t i2);
bool RandomSelectkLeaves(TTree * tree, Long64_t indexes[], int k);
void CommonAncestorPlot(TApplication * myapp, int nb_events, Double_t max_y, Double_t x01, Double_t rho);
void compute_biggest_child(TTree * tree, TApplication * myapp);
