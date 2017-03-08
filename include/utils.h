#include "color_codes.h"

#include <TApplication.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF2.h>

void decode_parameters(std::string filename, Double_t * rho, Double_t * max_y, Double_t * R, std::string * cutoff_type, int * nb_events);
TH1F * n_to_nbar(TH1F * old_hist);
void BinLogX(TH1 *h);
Axis_t * BinLogX2(TH1 *h, int bins);
Axis_t * BinLogX3(int nbins);
Double_t n(Double_t r, Double_t x01, Double_t y);

void general_plot(TApplication * myapp);
void fluctuations(Double_t max_y, Double_t x01, Double_t rho, Double_t r, 
                  const char * filename, 
                  const char * filename_hist, 
                  bool with_cutoff, bool minimal);
void draw_fluctuations(TApplication * myapp, const char * filename, 
                        bool logX, bool with_cutoff, 
                        Double_t x01, Double_t r, Double_t max_y, Double_t R);
void stat_events(TApplication * myapp, Double_t max_y, Double_t x01, const char * filename, bool minimal);
void draw_tree(TApplication * myapp, TTree * tree);
void generate_events(int nb_events, Double_t rho, Double_t max_y, Double_t R, bool with_cutoff = false, TF1 * cutoff = NULL, bool raw_cutoff = false, const char * tree_file = "tree.root", const char * lut_file = "lookup_table");
void generate_histograms(int nb_events, Double_t rho, Double_t max_y, Double_t R, bool with_cutoff, TF1 * cutoff, bool raw_cutoff, const char * histogram_file, const char * lut_file, std::vector<Double_t> r);

Long64_t GetCommonAncestors(TTree * tree, Long64_t i1, Long64_t i2);
bool RandomSelectkLeaves(TTree * tree, Long64_t indexes[], int k, bool minimal);
void CommonAncestorPlot(TApplication * myapp, int nb_events, Double_t max_y, Double_t x01, Double_t rho, const char * filename, bool minimal);
void compute_biggest_child(TTree * tree, TApplication * myapp);
