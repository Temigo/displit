//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Nov  8 23:36:38 2016 by ROOT version 6.09/01
// from TTree tree0/Dipole splitting
// found on file: tree.root
//////////////////////////////////////////////////////////

#ifndef EventTree_h
#define EventTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TVector2.h"

class EventTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        rapidity;
   Double_t        radius;
   Double_t        phi;
   TVector2        *coord;
   Long64_t        depth;
   Long64_t        index_children;
   Long64_t        index_parent;
   Bool_t          isLeaf;

   // List of branches
   TBranch        *b_rapidity;   //!
   TBranch        *b_radius;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_coord;   //!
   TBranch        *b_depth;   //!
   TBranch        *b_index_children;   //!
   TBranch        *b_index_parent;   //!
   TBranch        *b_isLeaf;   //!

   EventTree(TTree *tree=0);
   virtual ~EventTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef EventTree_cxx
EventTree::EventTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("tree.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("tree.root");
      }
      f->GetObject("tree0",tree);

   }
   Init(tree);
}

EventTree::~EventTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t EventTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t EventTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void EventTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   coord = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("rapidity", &rapidity, &b_rapidity);
   fChain->SetBranchAddress("radius", &radius, &b_radius);
   fChain->SetBranchAddress("phi", &phi, &b_phi);
   fChain->SetBranchAddress("coord", &coord, &b_coord);
   fChain->SetBranchAddress("depth", &depth, &b_depth);
   fChain->SetBranchAddress("index_children", &index_children, &b_index_children);
   fChain->SetBranchAddress("index_parent", &index_parent, &b_index_parent);
   fChain->SetBranchAddress("isLeaf", &isLeaf, &b_isLeaf);
   Notify();
}

Bool_t EventTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void EventTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t EventTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef EventTree_cxx
