#include <fstream>
#include <string>
#include <cmath>
#include <vector>

// LOOPP files
#include "BZUPusher.h"
#include "particle.h"

//#include "disp_rel_fluid.h"
#include "plasmaSlabField.h"


int main( int argc, char * argv []) {

  int i = 1;
  int NIonSpecies= atoi(argv[i++]);
  double alpha = atof(argv[i++]);             
  double dimless_wave_freq = atof(argv[i++]);
  double B0 = atof(argv[i++]); 
  double theta_sample = atof(argv[i++]);   
  int species_index = atoi(argv[i++]);    // index of speacies "a"

  // general timestepping variables
  int tstep_fac = atoi(argv[i++]);
//  int print_fac = atoi(argv[i++]);
   

  double X0 = atof(argv[i++]);
  double Y0 = atof(argv[i++]);
  double rho0 = atof(argv[i++]);
  double theta0 = atof(argv[i++]);

  // cartesian coordinates
  double x[3];
  double v[3];

  std::vector<Ion> ions;
  for(int j=0; j<NIonSpecies; j++) {
    Ion new_ion;
    new_ion.charge = atoi(argv[i++]);
    new_ion.mass = atof(argv[i++]);
    new_ion.density = atof(argv[i++]);

    ions.emplace_back(new_ion);   
  }
   
  PlasmaSlabField field = PlasmaSlabField(alpha, dimless_wave_freq, B0, species_index, ions); 

  x[0] = rho0 * cos(theta0)+ X0;
  x[1] = Y0 - rho0 * sin(theta0);
  x[2] = 0;

  for (int j = 1; j < i; j++){
    printf("%f ",atof(argv[j]));
  }
  printf("\n");

  BZUPusher pusher  = BZUPusher();

  Particle particle = Particle(0, std::string("dimless"), 1, 1, pusher, field); // ion

  particle.setX(x[0], x[1], x[2]);


  double position[3] = {0};
  double velocity[3] = {0};
  double A[3] = {0};

  double old_theta, rho, Y, new_theta;
  field.getVectorPotential(0, x, A);      

  v[0] = -rho0 * sin(theta0) - A[0];
  v[1] = X0 - A[1];
  v[2] = 0;

  particle.setV(v[0], v[1], v[2]);
  old_theta = theta0;

  double dt = M_PI / (100*tstep_fac);
  i=0;
  new_theta =   atan2(- (v[0]+A[0]),x[0] - (v[1]+A[1]));
  rho = sqrt(pow(v[0]+A[0],2) + pow(x[0] - (v[1]+A[1]),2));
  Y  = std::fmod(x[1] + rho * sin(new_theta), 2*M_PI);
  fprintf(stderr,"Y0 = %22.15e\n", Y);
  printf("%10u %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e\n", 0, i*dt,
             x[0], x[1], x[2], v[0], v[1], v[2],
              (v[0]+A[0]), (v[1]+A[1]), (v[2]+A[2]),
              rho, new_theta, Y);
  for (unsigned long int i = 0; i*dt<2*1e4; i++) {
  //  if (i % (tstep_fac * print_fac) == 0) {
      particle.getX(position);
      particle.getV(velocity);
      field.getVectorPotential(i*dt, position, A);      

      new_theta =   atan2(- (velocity[0]+A[0]),position[0] - (velocity[1]+A[1]));
      rho = sqrt(pow(velocity[0]+A[0],2) + pow(position[0] - (velocity[1]+A[1]),2));
      Y  = std::fmod(position[1]-dimless_wave_freq * i*dt + rho * sin(new_theta), 2*M_PI);
        
      if (old_theta < theta_sample && new_theta >= theta_sample) {
        printf("%10u %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e\n", 0, i*dt,
             position[0], position[1], position[2], velocity[0], velocity[1], velocity[2],
              (velocity[0]+A[0]), (velocity[1]+A[1]), (velocity[2]+A[2]),
              rho, new_theta, Y);
      }
  //}
    particle.push( dt);
    old_theta = new_theta;
  }
  return 0;
}


