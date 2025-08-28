#include "plasmaSlabField.h"

#include "geomUtils.h"

#include <math.h>
#include <cstdio>
#include <vector>

PlasmaSlabField::PlasmaSlabField (double alpha, double dimless_wave_freq, double B0, int species_index, std::vector<Ion>& ions)
    :alpha_(alpha),
    dimless_wave_freq_(dimless_wave_freq) {
   
  const double proton_charge_over_mass = 95788331.55943637;
  const double speed_of_light = 299792458.0;

  double reference_ion_cyclotron_frequency = B0*ions[species_index].charge*proton_charge_over_mass / ions[species_index].mass;
  
  double omega = dimless_wave_freq_ *  reference_ion_cyclotron_frequency;  // dimensional wave frequency is -k v gamma  = -v_ omeag_{ca} gamma

  StixDispRel dispersion(omega,  B0, ions);
  dispersion.GetParams(S_,D_,P_);
  s_ = S_/sqrt(S_*S_ + D_*D_);
  d_ = D_/sqrt(S_*S_ + D_*D_);
  
  fprintf(stderr,"S = %f, D = %f\n", S_,D_);

  if (S_ != 0){
    Ny_squared_ = S_ - pow(D_,2)/S_;
      fprintf(stderr,"N_y^2 = %f\n", Ny_squared_);
      if (Ny_squared_ > 0) {
        fprintf(stderr,"k_y = %f\n", omega*sqrt(Ny_squared_)/speed_of_light);
        prop_ = true;
      } else {
        fprintf(stderr,"Evanescent wave\n");
      }
  } else {
    fprintf(stderr,"Infinite N_y\n");
  }
}

PlasmaSlabField::~PlasmaSlabField() {}
/*
double PlasmaSlabField::f(double z) {
  if (z>= Lz_/2)    return 1.;
  if (z >= -Lz_/2)  return 0.5 + z/Lz_;
  return 0.;
}

double PlasmaSlabField::fprime(double z) {
   if (z>= Lz_/2)    return 0.;
   if (z >= -Lz_/2)  return 1/Lz_;
   return 0.;
}*/

void PlasmaSlabField::getFields(const double t, const double *xyz, double *E_xyz, double *B_xyz) {

   double x = xyz[0];
   double y = xyz[1];
   double z = xyz[2];
   
   if(prop_){
      B_xyz[0] = 0;
      B_xyz[1] = 0;
      B_xyz[2] = 1 + alpha_ * s_ * sin (y - dimless_wave_freq_ * t)/dimless_wave_freq_;


      E_xyz[0] = - alpha_ * s_ * sin(y - dimless_wave_freq_ * t);
      E_xyz[1] = alpha_ * d_ * cos(y - dimless_wave_freq_ * t);
      E_xyz[2] = 0;
   } 
}


void PlasmaSlabField::getVectorPotential( const double t, const double *xyz, double *A_xyz) {
   double x = xyz[0];
   double y = xyz[1];
   double z = xyz[2];

   A_xyz[0] = alpha_ * s_ * cos(y - dimless_wave_freq_ * t)/dimless_wave_freq_;
   A_xyz[1] = x+alpha_ * d_ * sin(y - dimless_wave_freq_ * t)/dimless_wave_freq_;
   A_xyz[2] = 0;

}

