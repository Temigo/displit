/* Implementation of dipole
*/
#include "dipole.h"

#include <iostream>

#include <TMath.h>
#include <TArrow.h>

Dipole::Dipole(Long64_t depth, Long64_t index) : 
    coord(0.5, 0.0), 
    phi(0.0), 
    rapidity(0.0), 
    radius(1.0), 
    /*parent(nullptr), 
    child1(nullptr), 
    child2(nullptr) */
    depth(depth),
    index(index),
    index_parent(0),
    index_children(-1),
    nb_left_brothers_split(0),
    nb_right_brothers(0),
    isLeaf(false)
{}

// FIXME : DEPRECATED
Long64_t Dipole::GetParentDepth()
{
    return depth-1;
}

// FIXME : DEPRECATED
Long64_t Dipole::GetParentIndex()
{
    return index / 2;
}

void Dipole::Draw()
{
    TArrow * ar = new TArrow(coord.X()-radius/2. *TMath::Cos(phi), coord.Y()-radius/2. *TMath::Sin(phi), coord.X() + radius/2. * TMath::Cos(phi), coord.Y() + radius/2. * TMath::Sin(phi), 0.05, "|");            
    ar->SetLineColor(2);
    if (!isLeaf) ar->SetLineColor(4);
    ar->Draw();
}