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
void fluctuations(Double_t max_y, Double_t x01, Double_t rho, Double_t r, 
                  const char * filename, 
                  const char * filename_hist, 
                  bool with_cutoff);
void draw_fluctuations(TApplication * myapp, const char * filename, 
                        bool logX, bool with_cutoff, 
                        Double_t x01, Double_t r, Double_t max_y);
void stat_events(TApplication * myapp, Double_t max_y, Double_t x01, const char * filename);
void draw_tree(TApplication * myapp, TTree * tree);
void generate_events(int nb_events, Double_t rho, Double_t max_y, Double_t R, bool with_cutoff = false, TF1 * cutoff = NULL, bool raw_cutoff = false, const char * tree_file = "tree.root", const char * lut_file = "lookup_table");
Long64_t GetCommonAncestors(TTree * tree, Long64_t i1, Long64_t i2);
bool RandomSelectkLeaves(TTree * tree, Long64_t indexes[], int k);
void CommonAncestorPlot(TApplication * myapp, int nb_events, Double_t max_y, Double_t x01, Double_t rho, const char * filename);
void compute_biggest_child(TTree * tree, TApplication * myapp);
