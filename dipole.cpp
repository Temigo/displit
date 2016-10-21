#include "dipole.h"

#include <iostream>
#include <TMath.h>
#include <TArrow.h>

Dipole::Dipole(Long64_t depth, Long64_t index) : 
    coord(0.5, 0.0), 
    phi(0.0), 
    rapidity(0.0), 
    radius(0.5), 
    /*parent(nullptr), 
    child1(nullptr), 
    child2(nullptr) */
    depth(depth),
    index(index),
    isLeaf(false)
{}

Long64_t Dipole::GetParentDepth()
{
    return depth-1;
}

Long64_t Dipole::GetParentIndex()
{
    return index / 2;
}

void Dipole::Draw()
{
    TArrow * ar = new TArrow(coord.X()-radius*TMath::Cos(phi), coord.Y()-radius*TMath::Sin(phi), coord.X() + radius * TMath::Cos(phi), coord.Y() + radius * TMath::Sin(phi), 0.05, "|");            
    ar->SetLineColor(2);
    ar->Draw();
}