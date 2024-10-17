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

  template <typename T> double vectorMag3DSquared(T &v) {
    return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
  }

public:
  auto fun(const int nlocal) {

    // PREAMBLE: Each particle in the system has its own density, radius, and
    // type (material) The type alternates between particles, and the counting
    // for types starts at 1!
    std::vector<std::array<double, 3>> v(
        nlocal, {1., 0.3 * log10(nlocal), 0.5 * log10(nlocal)});
    std::vector<double> density(nlocal, 1000.);
    std::vector<double> r(nlocal, 1e-3);
    std::vector<int> type({1, 2});

    std::vector<int> mask(nlocal);
    int half_nlocal = nlocal / 2;
    for (int i = 0; i < half_nlocal; ++i) {
      mask[i] = 1;
    }
    for (int i = half_nlocal; i < nlocal; ++i) {
      mask[i] = 3;
    }

    const int max_type = 2;
    std::vector<double> Y({1e6, 2e6});
    std::vector<double> nu(max_type, 0.2);
    std::vector<std::vector<double>> Yeff(max_type, {1.25e6, 1.75e6});
    int groupbit = 1;
    const double vmax_sqr_mesh = 0.;
    // END PREAMBLE

    // check rayleigh time and vmax of particles
    rayleigh_time = std::numeric_limits<double>::max();
    r_min = std::numeric_limits<double>::max();
    double vmax_sqr = 0.;
    double vmag_sqr;
    double rayleigh_time_i;

    auto start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < nlocal; i++) {
      double rad = r[i];

      double shear_mod = Y[type[i] - 1] / (2. * (nu[type[i] - 1] + 1.));
      rayleigh_time_i = M_PI * rad * sqrt(density[i] / shear_mod) /
                        (0.1631 * nu[type[i] - 1] + 0.8766);
      if (mask[i] & groupbit) {
        if (rayleigh_time_i < rayleigh_time)
          rayleigh_time = rayleigh_time_i;

        if (rad < r_min)
          r_min = rad;
      }
    }

    // check estimation for hertz time
    // this is not exact...
    // loop over all material combinations
    //  loop all particles
    //     test collision of particle with itself
    double hertz_time_min = 1000000.;
    double hertz_time_i, meff, reff;

    for (int ti = 1; ti < max_type + 1; ti++) {
      for (int tj = 1; tj < max_type + 1; tj++) {
        const double Eeff = Yeff[ti][tj];

        for (int i = 0; i < nlocal; i++) {
          // decide vmax - either particle-particle or particle-mesh contact
          vmag_sqr = vectorMag3DSquared(v[i]);
          if (vmag_sqr > vmax_sqr)
            vmax_sqr = vmag_sqr;
          double v_rel_max_simulation = std::max(
              2. * sqrt(vmax_sqr), sqrt(vmax_sqr) + sqrt(vmax_sqr_mesh));
          if (mask[i] & groupbit) {
            if (type[i] != ti || type[i] != tj)
              continue;
            meff = 4. * r[i] * r[i] * r[i] * M_PI / 3. * density[i];
            reff = r[i] / 2.;
            hertz_time_i =
                2.87 *
                pow(meff * meff / (reff * Eeff * Eeff * v_rel_max_simulation),
                    0.2);
            if (hertz_time_i < hertz_time_min)
              hertz_time_min = hertz_time_i;
          }
        }
      }
    }

    hertz_time = hertz_time_min;

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
