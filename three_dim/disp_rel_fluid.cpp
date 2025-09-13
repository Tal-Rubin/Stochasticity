#include "disp_rel_fluid.h"
#include <cmath>


StixDispRel::StixDispRel(double omega, double B, std::vector<Ion> ions): ions_(ions) {
  S_ = 1;
  D_ = 0;
  P_ = 1;

  const double proton_charge_over_mass = 95788331.55943637;
  const double electron_charge_over_mass = -175882001077.2163;

  const double proton_plasma_frequency_squared = 1.733302137690269e+20;
  const double electron_plasma_frequency_squared = 3.182607353999257e+23;

  double electron_gyrofrequency = electron_charge_over_mass * B;

   if (omega == 0){
      double gamma_A = 0;
      
      for (Ion& ion :ions_) {
        double electron_density = ion.density*ion.charge;
        double ion_plasma_frequency_squared = pow(ion.charge, 2)*ion.density * proton_plasma_frequency_squared/ion.mass;
        double ion_gyrofrequency = ion.charge * proton_charge_over_mass * B / ion.mass;

        gamma_A += ion_plasma_frequency_squared * pow(ion_gyrofrequency,-2);
        gamma_A += electron_density * electron_plasma_frequency_squared * pow(electron_gyrofrequency,-2);
      }

      S_ += gamma_A;
   } else {
    for (Ion& ion :ions_) {
      double electron_density = ion.density*ion.charge;
      double ion_plasma_frequency_squared = pow(ion.charge, 2)*ion.density * proton_plasma_frequency_squared/ion.mass;
      double ion_gyrofrequency = ion.charge * proton_charge_over_mass * B / ion.mass;

      S_ += -ion_plasma_frequency_squared /(pow(omega,2)- pow(ion_gyrofrequency,2));
      S_ += -electron_density * electron_plasma_frequency_squared /(pow(omega,2)- pow(electron_gyrofrequency,2));

      double Dtemp = (ion_gyrofrequency * ion_plasma_frequency_squared /(pow(omega,2)- pow(ion_gyrofrequency,2))) 
                            + (electron_gyrofrequency * electron_density * electron_plasma_frequency_squared /(pow(omega,2)- pow(electron_gyrofrequency,2)));
      D_ += Dtemp / omega;
      
      //D_ += ion_gyrofrequency * ion_plasma_frequency_squared /(pow(omega,2)- pow(ion_gyrofrequency,2))/omega;
      //D_ += electron_gyrofrequency * electron_density * electron_plasma_frequency_squared /(pow(omega,2)- pow(electron_gyrofrequency,2))/omega;

      P_ += -ion_plasma_frequency_squared /pow(omega,2);
      P_ += -electron_density * electron_plasma_frequency_squared /pow(omega,2);
    }
  }
}


StixDispRel::~StixDispRel() {}

void StixDispRel::GetParams(double &S, double &D, double P) {
  S = S_;
  D = D_;
  P = P_;
}