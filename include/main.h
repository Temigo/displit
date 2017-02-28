#include <TROOT.h>

std::string encode_parameters(int nb_events, Double_t rho, Double_t max_y, Double_t R, std::string cutoff_type, std::string optional);
void decode_parameters(std::string filename, Double_t * rho, Double_t * max_y, Double_t * R, std::string * cutoff_type, int * nb_events);
int main( int argc, char* argv[] );