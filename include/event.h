#include "dipole.h"

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
    TTree * tree;
    static constexpr char * filename = "tree.root";

    // Lookup table
    ROOT::Math::Interpolator interpolator;
    std::map<Double_t, Double_t> lookup_table;
    Double_t x01_min, x01_max;

    public:
        Event(Double_t rho, Double_t max_y, TF2 * f);
        //Event(const char * filename, TTree * tree);

        static Double_t f(Double_t  * r, Double_t * parameters = NULL);
        Double_t g(Double_t r);
        Double_t r_generate();
        Double_t r_generate_cutoff();
        Double_t theta(Double_t r);
        Double_t lambda(Double_t x01);
        Double_t getLambda(Double_t x01);
        void WriteLookupTable(const char * filename);
        void LoadLookupTable(const char * filename);
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
        
        TTree * make_tree(const char * filename = filename, const char * treename = "T", bool draw = false, bool draw_step_by_step = false);
};