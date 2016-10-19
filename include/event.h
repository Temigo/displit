#include "dipole.h"

#include <TTree.h>

class Event
{
    Double_t rho;
    Double_t max_y;
    static constexpr char * filename = "tree.root";

    public:
        Event();

        Double_t f(Double_t r);
        Double_t g(Double_t r);
        Double_t r_generate(Double_t rho);
        Double_t phi(Double_t r);
        Double_t theta(Double_t r);
        Double_t lambda(Double_t rho);
        Double_t y_generate(Double_t rho);

        void draw(Double_t x, Double_t y, Double_t rapidity);
        void draw_tree(TTree * tree);
        void generate(Double_t rho, Dipole * dipole, Dipole * dipole1, Dipole * dipole2);
        void generate_normalized(Double_t rho, Double_t * x, Double_t * y, Double_t * rapidity);
        void bare_distribution();
        Long64_t make_tree(const char * filename = filename);
};