#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

//const maximum number of types
constexpr int max_type = 2;

//Unknown purpose?
constexpr int groupbit = 1;

//Class declaration
class MyClass {
  double rayleigh_time;
  double r_min;
  double hertz_time;

  // Method to compute square of velocity comp
   template <typename T> double vectorMag3DSquared(const T &v) const {
    return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
  }
   
public:
    double fun(const int);
};

//Explicit specialization
template <> inline double MyClass::vectorMag3DSquared(const std::array<double,3> &v) const {
    return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
}

double MyClass::fun(const int nlocal) {
    //nlocal = 10...10^7

    // PREAMBLE: Each particle in the system has its own density, radius, and
    // type (material) The type alternates between particles, and the counting
    // for types starts at 1!

    //A vector v of size nlocal that contains 3 double values in all elements
    //Possibly to store the 3 vel components of all particles
    std::vector<std::array<double, 3>> v(
        nlocal, {1., 0.3 * log10(nlocal), 0.5 * log10(nlocal)});

    //A vector 'density' of size nlocal that is uniformly initialized to 1000.   
    std::vector<double> density(nlocal, 1000.);

    //A vector 'r' of size nlocal that is uniformly initialized to 1e-3.    
    //Possibly denoting radius of particles
    std::vector<double> r(nlocal, 1e-3);
   
    //A vector 'mask' of size nlocal
    std::vector<int> mask(nlocal);

    //const maximum number of types
    //const int max_type = 2;

    //A vector 'type' that contains just 2 integers
    //Possibly denoting types of particles
    const std::array<int, max_type> type({1, 2});
    
    //Vector 'Y' storing 2 values, possibly young's modulus of 2 'types' of particles
    const std::array<double, max_type> Y({1e6, 2e6});
    
    //vector 'nu' storing 2 similar values, possibly viscosity of the particles
    const std::array<double, max_type> nu({0.2, 0.2});
    
    //A vector of vector of size 2 that has 2 values each; possibly effective young's modulus
    const std::array<std::array<double, max_type>, max_type> Yeff{{{1.25e6, 1.75e6}, {1.25e6, 1.75e6}}};

    //Variables to store the coefficients used within the Rayleigh time compute loop
    std::array<double, max_type>coeff1_nu({0.});
    std::array<double, max_type>coeff2_nu({0.});
    std::array<double, max_type>coeff_shearMod({0.});

    //Store '1' for half the points and '3' for the remaining
    int half_nlocal = nlocal / 2;
    for (int i = 0; i < half_nlocal; ++i) {
      mask[i] = 1;
    }
    for (int i = half_nlocal; i < nlocal; ++i) {
      mask[i] = 3;
    }
    // END PREAMBLE

    // check rayleigh time and vmax of particles
    //Set rayleigh time to max possible value 
    rayleigh_time = std::numeric_limits<double>::max();

    //Set value of r_min to max possible double value
    r_min = std::numeric_limits<double>::max();
    
    //Initializers 
    double rayleigh_time_i;
    double rad = 0.;
    double shear_mod =0.;
    double nuForEachParticle =0.;
    int typeIndex=-2;
    
    //Initialize hertz time
    double hertz_time_min = 1000000.;

    //Declare local var for the loop
    double hertz_time_i, meff, reff;

    //Initializer of max vel of the mesh, possibly mesh based velocity
    double vmax_sqr = 0.;
    double vmag_sqr;
    //const double vmax_sqr_mesh = 0.;
    const double sqrt_vmax_sqr_mesh = 0.;
    //Declare sqr of vmax_sqr
    double coeff_sqrt_vmax_sqr = 0.;
    double coeff_sumOfvMaxAndvMaxMesh = 0.;
    double v_rel_max_simulation = 0.;
    //Coefficients for meff
    const double coeff_meff = M_PI/3.0;
    //Eff
    double Eeff = 0.;

  
    //Clock start
    auto start = std::chrono::high_resolution_clock::now();

    //Find rmin
    auto minIterator = std::min_element(r.begin(), r.end());
    r_min = *minIterator;

    //Loop over the types of particles, precompute and store:
    //coeff involving nu
    //Shear Mod
    for (int i = 0; i < max_type; ++i) {
      coeff1_nu[i] = 1./(2. * (nu[i] + 1.));
      coeff2_nu[i]= 1./(0.1631 * nu[i] + 0.8766);
      coeff_shearMod[i] = 1.0/(Y[i]*coeff1_nu[i]);
    }

    //Loop over nlocal particles
    for (int i = 0; i < nlocal; i++) {
      //type access index 
      typeIndex = i%type.size();

      //Compute Rayleigh time 
      rayleigh_time_i = M_PI * r[i] * sqrt(density[i] * coeff_shearMod[typeIndex]) * coeff2_nu[typeIndex];
      
      //This if condition works for both types of particles
      //?? Looks redundant for the given type mask and groupbit values
      if (mask[i] & groupbit) {
        rayleigh_time = std::min(rayleigh_time, rayleigh_time_i);
      }
    }

    // check estimation for hertz time
    // this is not exact...
    // loop over all material combinations
    //  loop all particles
    //     test collision of particle with itself
    
    //For the particle pair {1,1}
    Eeff = Yeff[0][0];
    int ti = 1; int tj=1;

    //Loop over number of particles
    for (int i = 0; i < nlocal; i++) {
      
      //type access index 
      typeIndex = i%type.size();

      // decide vmax - either particle-particle or particle-mesh contact
      vmag_sqr = vectorMag3DSquared(v[i]);
      vmax_sqr = std::max(vmag_sqr, vmax_sqr);
      coeff_sqrt_vmax_sqr = 2.* sqrt(vmax_sqr);
      coeff_sumOfvMaxAndvMaxMesh = 0.5*coeff_sqrt_vmax_sqr + sqrt_vmax_sqr_mesh;

      //Compute rel max velocity in the simulation between particle and mesh vel
      v_rel_max_simulation = std::max(coeff_sqrt_vmax_sqr, coeff_sumOfvMaxAndvMaxMesh );
      
      //This branch works for all types
      //?? Redundant:Can we remove it?
      if (mask[i] & groupbit) {

        if (type[typeIndex] != ti || type[typeIndex] != tj)
          continue;
        
        //Compute effective mass
        meff = 4. * r[i] * r[i] * r[i] * coeff_meff * density[i];
        
        //Compute effective radius
        reff = r[i] * 0.5;

        //Compute hertz time
        hertz_time_i = meff * meff / (reff * Eeff * Eeff * v_rel_max_simulation);
        //Find min hertz time among all particles 
        hertz_time_min = std::min(hertz_time_i, hertz_time_min);  
      }
    }
    //For the particle pair {1,2}
    Eeff = Yeff[0][1];
    ti = 1; tj=2;
     
    //Loop over number of particles
    for (int i = 0; i < nlocal; i++) {
      
      //type access index 
      typeIndex = i%type.size();

      // decide vmax - either particle-particle or particle-mesh contact
      vmag_sqr = vectorMag3DSquared(v[i]);
      vmax_sqr = std::max(vmag_sqr, vmax_sqr);
      coeff_sqrt_vmax_sqr = 2.* sqrt(vmax_sqr);
      coeff_sumOfvMaxAndvMaxMesh = 0.5*coeff_sqrt_vmax_sqr + sqrt_vmax_sqr_mesh;

      //Compute rel max velocity in the simulation between particle and mesh vel
      v_rel_max_simulation = std::max(coeff_sqrt_vmax_sqr, coeff_sumOfvMaxAndvMaxMesh );
      
      //This branch works for all types
      //??Redundant: Can we remove it?
      if (mask[i] & groupbit) {

        if (type[typeIndex] != ti || type[typeIndex] != tj)
          continue;
        
        //Compute effective mass
        meff = 4. * r[i] * r[i] * r[i] * coeff_meff * density[i];
        
        //Compute effective radius
        reff = r[i] * 0.5;

        //Compute hertz time 
        hertz_time_i = meff * meff / (reff * Eeff * Eeff * v_rel_max_simulation);
        //Find min hertz time among all particles 
        hertz_time_min = std::min(hertz_time_i, hertz_time_min);  
      }
    }
    //For the particle pair {2,1}
    Eeff = Yeff[1][0];
    ti = 2; tj=1;

     //Loop over particles
    for (int i = 0; i < nlocal; i++) {
    
        //type access index 
        typeIndex = i%type.size();

        // decide vmax - either particle-particle or particle-mesh contact
        vmag_sqr = vectorMag3DSquared(v[i]);
        vmax_sqr = std::max(vmag_sqr, vmax_sqr);
        coeff_sqrt_vmax_sqr = 2.* sqrt(vmax_sqr);
        coeff_sumOfvMaxAndvMaxMesh = 0.5*coeff_sqrt_vmax_sqr + sqrt_vmax_sqr_mesh;

        //Compute rel max velocity in the simulation between particle and mesh vel
        v_rel_max_simulation = std::max(coeff_sqrt_vmax_sqr, coeff_sumOfvMaxAndvMaxMesh );
        
        //This branch works for all types
        //??Can we remove it?
        if (mask[i] & groupbit) {

          if (type[typeIndex] != ti || type[typeIndex] != tj)
            continue;
          
          //Compute effective mass
          meff = 4. * r[i] * r[i] * r[i] * coeff_meff * density[i];
          
          //Compute effective radius
          reff = r[i] * 0.5;

          //Compute hertz time
          hertz_time_i = meff * meff / (reff * Eeff * Eeff * v_rel_max_simulation);
          //Find min hertz time among all particles 
          hertz_time_min = std::min(hertz_time_i, hertz_time_min);  
        }
      }
    
    //For the particle pair {2,2}
    Eeff = Yeff[1][1];
    ti = 2; tj=2;
    
    //Loop over number of particles
    for (int i = 0; i < nlocal; i++) {
      
      //type access index 
      typeIndex = i%type.size();

      // decide vmax - either particle-particle or particle-mesh contact
      vmag_sqr = vectorMag3DSquared(v[i]);
      vmax_sqr = std::max(vmag_sqr, vmax_sqr);
      coeff_sqrt_vmax_sqr = 2.* sqrt(vmax_sqr);
      coeff_sumOfvMaxAndvMaxMesh = 0.5*coeff_sqrt_vmax_sqr + sqrt_vmax_sqr_mesh;

      //Compute rel max velocity in the simulation between particle and mesh vel
      v_rel_max_simulation = std::max(coeff_sqrt_vmax_sqr, coeff_sumOfvMaxAndvMaxMesh );
      
      //This branch works for all types
      //??Can we remove it?
      if (mask[i] & groupbit) {

        if (type[typeIndex] != ti || type[typeIndex] != tj)
          continue;
        
        //Compute effective mass
        meff = 4. * r[i] * r[i] * r[i] * coeff_meff * density[i];
        
        //Compute effective radius
        reff = r[i] * 0.5;

        //Compute hertz time 
        hertz_time_i = meff * meff / (reff * Eeff * Eeff * v_rel_max_simulation);
        //Find min hertz time among all particles 
        hertz_time_min = std::min(hertz_time_i, hertz_time_min);  
      }
    }
    
    //Store the final minimal hertz time as the simulation deltaT
    hertz_time = 2.87*pow(hertz_time_min, 0.2);

    //End timer and compute elapsed time
    auto end = std::chrono::high_resolution_clock::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    double duration_count = static_cast<double>(duration.count()) / 1000;

    std::cout << std::setw(10) << nlocal << std::setw(15) << hertz_time
              << std::setw(15) << rayleigh_time << std::setw(15)
              << duration_count << "\n";
    return duration_count;
  }


int main() {
  //Class instantiation
  MyClass c;

  //??This start seems to be unused.
  auto start = std::chrono::high_resolution_clock::now();
  //Console output
  std::cout << std::setw(10) << "nlocal" << std::setw(15) << "hertz_time"
            << std::setw(15) << "rayleigh_time" << std::setw(15)
            << "run_time [ms]" << "\n";
  //Initialize duration
  int duration_sum = 0;
  //Vary the number of points in the powers of 10 and call the fun() with it as arguments
  for (int i = 1; i < 8; ++i) {
    duration_sum += c.fun(pow(10, i));
  }
  //Print the total elapsed time
  std::cout << std::setw(40) << "" << std::setw(15) << duration_sum << "\n";
  return 0;
}
