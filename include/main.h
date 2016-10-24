#include <TApplication.h>

void general_plot(TApplication * myapp);
void stat_events(TApplication * myapp, int nb_events, Double_t max_y);
void generate_events(int nb_events, Double_t rho, Double_t max_y);
int main( int argc, const char* argv[] );