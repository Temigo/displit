#include "dipole.h"
#include "color_codes.h"

#include <TTree.h>
#include <map>
#include <Math/Interpolator.h>
#include <TF2.h>

Double_t phi(Double_t r);

class Event
{
    Double_t rho; // Cut-off ultraviolet
    Double_t max_y; // Maximal rapidity
    TF2 * cutoff;
    TF1 * f_cutoff; // f(r) distribution for size of dipoles
    TTree * tree;
    static constexpr char * filename = "tree.root";
    const char * lut_filename; // Lookup Table filename
    bool WITH_CUTOFF;
    Double_t R = 2.0; // cutoff for large sizes

    // Lookup table
    ROOT::Math::Interpolator interpolator;
    std::map<Double_t, Double_t> lookup_table;
    Double_t x01_min, x01_max;

    public:
        Event(Double_t rho, Double_t max_y, const char * lut_filename, TF2 * f = NULL, bool with_cutoff = false);
        //Event(const char * filename, TTree * tree);

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

        void draw(Double_t x, Double_t y, Double_t rapidity);
        void draw_tree(TTree * tree);
        bool generate(Dipole * dipole, Dipole * dipole1, Dipole * dipole2, Double_t max_y);

        void generate_normalized(Double_t * x, Double_t * y, Double_t * rapidity);
        void bare_distribution();
        void fit_r();
        void fit_y();
        void fit_x();
        
        TTree * make_tree(const char * treename = "T", bool draw_dipole = false, bool draw_step_by_step = false);
};