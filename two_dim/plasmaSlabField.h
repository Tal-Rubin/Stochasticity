#ifndef PLASMASLABFIELD_
#define PLASMASLABFIELD_
#include <vector>
#include "cartField.h"
#include "disp_rel_fluid.h"


//! Rotating magnetic field in cylindrical geometry.
/*!
   User must subclass, and implement / initialize all protected variables.
*/
class PlasmaSlabField : public CartField {

   public:   
      PlasmaSlabField (double alpha, double dimless_wave_freq, double B0, int species_index, std::vector<Ion>& ions);
      ~PlasmaSlabField();
      void getVectorPotential(const double t, const double *xyz, double *A_xyz);
    /*  double f(double z);
      double fprime(double z);
      double getbeta();
      bool isevan();*/

   protected:
      void getFields(const double t, const double *xyz, double *E_xyz, double *B_xyz);
      //double getTimeScaleRTZ_ ( const double t, const double *rtz);

      double alpha_, dimless_wave_freq_,S_, D_,P_ ;

      double Ny_squared_, s_, d_;
      bool prop_;
};

#endif  // PLASMASLABFIELD_