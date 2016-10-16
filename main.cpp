#include "main.h"
#include "event.h"

#include <TApplication.h>
#include <TFile.h>
#include <TTree.h>

int main( int argc, const char* argv[] )
{
    TApplication * myapp = new TApplication("myapp", 0, 0);
    //connect(myapp,SIGNAL(lastWindowClosed()),TQtRootSlot::CintSlot(),SLOT(TerminateAndQuit());
    //p3a();

    Event * e = new Event();
    e->make_tree("tree.root");
    //draw_tree(tree);
    //tree->Draw("radius");

    TFile f("tree.root");
    TTree * tree;
    f.GetObject("T", tree);
    tree->Print();
    //tree->Draw("coord.X():coord.Y()", "coord.X()>-2. && coord.X()<2. && coord.Y() < 2. && coord.Y() > -2.");
    //tree->Draw("log(rapidity):coord.X()", "coord.X() > -2. && coord.X() < 2.");
    tree->Draw("radius");

    myapp->Run();
    return 0;
}