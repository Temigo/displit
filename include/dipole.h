#include <TVector2.h>

class Dipole {
    public:
        TVector2 coord; // vector for the dipole's center
        Double_t phi; // angle of the dipole in its own referential
        Double_t rapidity; 
        Double_t radius; // length of dipole

        /*Dipole * parent;
        Dipole * child1;
        Dipole * child2;*/

        // Because we cannot store custom data types in TTree
        // Binary tree => use indexes
        // Long64_t : from 0 to 2^64-1
        Long64_t depth;
        Long64_t index;
        Long64_t index_parent;
        Long64_t index_children;
        
        Long64_t nb_left_brothers_split;
        Long64_t nb_right_brothers;

        Bool_t isLeaf;

        Dipole(Long64_t depth, Long64_t index);
        Long64_t GetParentDepth();
        Long64_t GetParentIndex();
        void Draw();
};