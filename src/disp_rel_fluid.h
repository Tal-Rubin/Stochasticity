#ifndef DISP_REL_FLUID_
#define DISP_REL_FLUID_

#include <vector>

struct Ion
{
    int charge;     // in elementary charge units
    double mass;    // in atomic mass units
    double density; // in units of 1e20 [1/m^3]
};

class StixDispRel {
  
  public:   
      StixDispRel(double omega, double B, std::vector<Ion> ions);
      ~StixDispRel();
      void GetParams(double &S, double &D, double P);

  protected:
    std::vector<Ion> ions_;
    double S_;
    double D_;
    double P_;
};
#endif