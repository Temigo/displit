#include "dipole.h"

Dipole::Dipole(Long64_t depth, Long64_t index) : 
    coord(0.5, 0.5), 
    phi(0.0), 
    rapidity(0.0), 
    radius(1.0), 
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