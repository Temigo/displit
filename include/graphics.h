#include <TApplication.h>
#include <TF2.h>

void draw_cutoffs(TApplication * myapp, std::map<std::string, TF1 *> cutoffs, 
                    std::map<std::string, int> colors);
void fit_bare_r(Double_t rho, Double_t max_y, TApplication * myapp);
void fit_fluctuations(Double_t rho, Double_t max_y, TApplication * myapp);
void compare_histo(TApplication * myapp, std::string histofiles,
                    std::map<std::string, TF1 *> cutoffs);
void compare_c(TApplication * myapp, std::string histofiles);