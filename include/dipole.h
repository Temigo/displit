#include <TVector2.h>

class Dipole {
    public:
        TVector2 coord; // vector for the dipole's center
        Double_t phi; // angle of the dipole in its own referential
        Double_t rapidity; 
        Double_t radius; // half length of dipole

        Dipole * parent;
        Dipole * child1;
        Dipole * child2;

        Dipole();
};