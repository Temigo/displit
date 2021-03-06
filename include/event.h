#include "dipole.h"
#include "color_codes.h"

#include <TTree.h>
#include <Math/Interpolator.h>
#include <TF2.h>
#include <TF12.h>
#include <TMath.h>
#include <TH1F.h>

#include <map>
#include <vector>

Double_t phi(Double_t r);

struct IntegralFunction
{
    TF2 * function;
    TF12 * f12;
    const char * name;
    IntegralFunction(TF2 * f, const char * name);
    double operator() (double * x, double * p) const;
};

class Event
{
    Double_t rho; // Cut-off ultraviolet
    Double_t max_y; // Maximal rapidity
    TF1 * cutoff;
    TF2 * integrand;
    IntegralFunction * integralFunction;
    TF1 * f_cutoff; // f(r) distribution for size of dipoles
    TTree * tree;
    static constexpr char * filename = "tree.root";
    const char * lut_filename; // Lookup Table filename
    bool WITH_CUTOFF;
    bool RAW_CUTOFF;
    bool MINIMAL;
    Double_t R; // cutoff for large sizes

    // Lookup table
    ROOT::Math::Interpolator interpolator;
    std::map<Double_t, Double_t> lookup_table;
    Double_t x01_min, x01_max;

    public:
        Event(Double_t rho, Double_t max_y, Double_t R, const char * lut_filename, TF1 * f = NULL, bool with_cutoff = false, bool raw_cutoff=false, bool minimal=false);
        //Event(const char * filename, TTree * tree);
        ~Event();

        static Double_t f(Double_t  * r, Double_t * parameters = NULL);
        Double_t g(Double_t r);
        Double_t r_generate(Double_t x01 = 1.0, bool display_cutoff = false);
        Double_t theta(Double_t r, Double_t x01 = 1.0);
        Double_t lambda(Double_t x01);
        Double_t getLambda(Double_t x01);
        void WriteLookupTable();
        void LoadLookupTable();
        void PrintLookupTable();
        void SetInterpolatorData();
        Double_t y_generate(Double_t x01);
        Double_t r_generate2(Double_t x01 = 1.0);

        void draw(Double_t x, Double_t y, Double_t rapidity);
        bool generate(Dipole * dipole, Dipole * dipole1, Dipole * dipole2, Double_t max_y);

        void generate_normalized(Double_t * x, Double_t * y, Double_t * rapidity);
        void bare_distribution();
        void fit_r();
        void fit_y();
        void fit_x();
        
        void make_tree(TTree * tree, bool draw_dipole = false, bool draw_step_by_step = false);
        std::vector<int> make_histograms(std::vector<Double_t> r);
};