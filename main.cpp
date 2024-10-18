#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

/*

Have a look at the following function:

https://github.com/CFDEMproject/LIGGGHTS-PUBLIC/blob/master/src/fix_check_timestep_gran.cpp#L233

Task 1: Please explain what this function calculates, why is it needed?

This function was extracted into the code below. We have removed some
code related to meshes and adapted some code to make the code compile and output
something.

Task 2: Test the code and fix any obvious bugs you find. Explain what tools you
used to find the bugs.

Task 3: There are two for loops over all particles in the function
"fun". Do you see any way of optimizing these? If so, implement the
optimization such that the code runs as fast as possible (while maintaining
correct results). Please attach the output (timings) you get for your compiled
code, and also indicate what compilation command / flags you used

Bonus task 4: There is a comment saying "// this is not exact...". Can
you identify what this comment refers to? I.e. what is not exact?

*/

class MyClass {
  double rayleigh_time;
  double r_min;
  double hertz_time;

  // Method to compute square of velocity comp
  //??Not const-qualified for arguments or function signature
  //??Possibly inline to prevent repeated function pointer indirection  
  template <typename T> double vectorMag3DSquared(T &v) {
    return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
  }

public:
  auto fun(const int nlocal) {
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

    //A vector 'type' that contains just 2 integers
    //Possibly denoting types of particles
    std::vector<int> type({1, 2});
    
    //A vector 'mask' of size nlocal
    std::vector<int> mask(nlocal);
    
    //Store '1' for half the points and '3' for the remaining
    int half_nlocal = nlocal / 2;
    for (int i = 0; i < half_nlocal; ++i) {
      mask[i] = 1;
    }
    for (int i = half_nlocal; i < nlocal; ++i) {
      mask[i] = 3;
    }
    
    //const maximum number of types
    const int max_type = 2;
    
    //Vector 'Y' storing 2 values, possibly young's modulus of 2 'types' of particles
    std::vector<double> Y({1e6, 2e6});
    //vector 'nu' storing 2 similar values, possibly viscosity of the particles
    std::vector<double> nu(max_type, 0.2);
    //A vector of vector of size 2 that has 2 values each; possibly effective young's modulus
    std::vector<std::vector<double>> Yeff(max_type, {1.25e6, 1.75e6});

    //?? What is this used for?
    int groupbit = 1;

    //Initializer of max vel of the mesh, possibly mesh based velocity
    const double vmax_sqr_mesh = 0.;
    // END PREAMBLE

    // check rayleigh time and vmax of particles
    //Set rayleigh time to max possible value 
    rayleigh_time = std::numeric_limits<double>::max();

    //Set value of r_min to max possible double value
    r_min = std::numeric_limits<double>::max();
    
    //Initializer
    double vmax_sqr = 0.;
    double vmag_sqr;
    double rayleigh_time_i;
    
    //Clock start
    auto start = std::chrono::high_resolution_clock::now();

    //Loop over nlocal particles
    //??Think about thread parallelism here; each particle to a thread with a parallel red?
    for (int i = 0; i < nlocal; i++) {
      //Read and store radius of the ith particle
      //?? declaration of double can be shifted outside?
      double rad = r[i];

      // Compute shear modulus of the ith particle
      //?? Shift declaration of the shear_mod var to outside the loop
      //++ Why is type[i] not showing seg fault for type[2] when i=2 : RESOLVED
      //?? Is Y[type[i]] and nu[type[]] optimal in memory access?
      double shear_mod = Y[type[i] - 1] / (2. * (nu[type[i] - 1] + 1.));

      //Compute rayleigh time of the ith particle
      //?? is nu[type[i]] optimal in memory access?
      rayleigh_time_i = M_PI * rad * sqrt(density[i] / shear_mod) /
                        (0.1631 * nu[type[i] - 1] + 0.8766);
      
      //This if condition works for both types of particles
      //?? Can be removed since it is redundant
      if (mask[i] & groupbit) {
        //Tries to find out if the computed time is smaller than the stored time and replaces it if true
        //?? Can we do better than an if?
        if (rayleigh_time_i < rayleigh_time)
          rayleigh_time = rayleigh_time_i;
        
        //Tries to find out if the computed rad is smaller than the obtained rad and replaces it if true
        //?? can we do better than an if?
        //Since rad is not computed, move the whole computation outside the loop and use std::min
        if (rad < r_min)
          r_min = rad;
      }
    }

    // check estimation for hertz time
    // this is not exact...
    // loop over all material combinations
    //  loop all particles
    //     test collision of particle with itself
    
    //Initialize hertz time
    double hertz_time_min = 1000000.;

    //Declare local var for the loop
    double hertz_time_i, meff, reff;
    
    //?? Can we do better than 4 iterations here? Maybe a lookup table?
    //First loop over the type of particles {1,2}   
    for (int ti = 1; ti < max_type + 1; ti++) {
      // Second loop over the type of particles {1,2}
      for (int tj = 1; tj < max_type + 1; tj++) {

        //Store the effective shear_mod of the particle pairs [1][1], [1][2], [2][1], [2][2]
        //??Seg fault since Yeff[0] Yeff[1] only exists by declaration
        //?? move variable declarations outisde the loop
        const double Eeff = Yeff[ti-1][tj-1];

        //Loop over number of particles
        for (int i = 0; i < nlocal; i++) {
          // decide vmax - either particle-particle or particle-mesh contact
          vmag_sqr = vectorMag3DSquared(v[i]);
          if (vmag_sqr > vmax_sqr)
            vmax_sqr = vmag_sqr;
          
          //Compute rel max velocity in the simulation between particle and mesh vel
          //?? Move variable declaration outside
          //??sqrt() could be expensive within the loop. 
          double v_rel_max_simulation = std::max(
              2. * sqrt(vmax_sqr), sqrt(vmax_sqr) + sqrt(vmax_sqr_mesh));
          
          //This branch works for all types
          //??Can we remove it?
          if (mask[i] & groupbit) {

            if (type[i] != ti || type[i] != tj)
              continue;
            
            //Compute effective mass
            //??Avoid divisions  
            meff = 4. * r[i] * r[i] * r[i] * M_PI / 3. * density[i];
            //Compute effective radius
            //??Avoid divisions
            reff = r[i] / 2.;
            
            //Compute hertz time
            //??Avoid pow() and divisions
            hertz_time_i =
                2.87 *
                pow(meff * meff / (reff * Eeff * Eeff * v_rel_max_simulation),
                    0.2);
            //Find min hertz time among all particles
            //??Can we do better than 'if' branch?
            if (hertz_time_i < hertz_time_min)
              hertz_time_min = hertz_time_i;
          }
        }
      }
    }
    //Store the final minimal hertz time as the simulation deltaT
    hertz_time = hertz_time_min;

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
};

int main() {
  MyClass c;
  auto start = std::chrono::high_resolution_clock::now();
  std::cout << std::setw(10) << "nlocal" << std::setw(15) << "hertz_time"
            << std::setw(15) << "rayleigh_time" << std::setw(15)
            << "run_time [ms]" << "\n";
  int duration_sum = 0;
  for (int i = 1; i < 8; ++i) {
    duration_sum += c.fun(pow(10, i));
  }
  std::cout << std::setw(40) << "" << std::setw(15) << duration_sum << "\n";
  return 0;
}
