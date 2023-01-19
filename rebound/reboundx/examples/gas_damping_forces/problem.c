/**
 * Gas Damping of Eccentricity and Inclination
 * 
 * This example shows how to add eccentricity and inclination damping to bodies 
 * in the presence of gas.  The calculations of the damping timescale follow
 * Eq. 16 in Dawson et al. 2016.  You can adjust the variable "d" [depletion factor]
 * to either increase or decrease the amount of gas present.
 * 
 * d=1 corresponds to the minimum mass solar nebula, and d>1 corrsponds to a more
 * depleted nebula (less gas).
 *
 * 
 * code unit: AU, 𝑀⊙, yr/2𝜋, G = 1
 * The central body is 1 solar mass.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include "rebound.h"
#include "reboundx.h"

double tmax;
void heartbeat(struct reb_simulation* sim);

int main(int argc, char* argv[]){
    
    // initialize your simulation
    struct reb_simulation* sim = reb_create_simulation();

    // set up constants
    sim->integrator = REB_INTEGRATOR_MERCURIUS;  // choose your integrator
    sim->dt = 0.0137*2.*M_PI;    // choose your initial timestep [this corresponds to 5 days]
    sim->heartbeat = heartbeat;

    // add a host star
    struct reb_particle star = {0};
    star.m      = 1.;
    star.r      = 0.005;
    star.hash   = reb_hash("CentralStar");
    reb_add(sim, star);

    // choose your initial conditions for your planets
    double a1 = 1.; // Semi-major axis (au)
    double a2 = 2.;
    double e1 = 0.05; // Eccentricity
    double e2 = 0.2;
    double inc1 = 0.5; // Inc 
    double inc2 = 0.;
    double Omega = 0.; // Argument of peri 
    double omega = 0.; // Longitude of ascending node 
    double f = 0.; // true anomaly 
    double m1 = 3.e-6; // Mass in solar
    double m2 = 1.e-6;

    
    // add planets
    struct reb_particle p1 = reb_tools_orbit_to_particle(sim->G, star, m1, a1, e1, inc1, Omega, omega, f);
    struct reb_particle p2 = reb_tools_orbit_to_particle(sim->G, star, m2, a2, e2, inc2, Omega, omega, f);
    reb_add(sim, p1); 
    reb_add(sim, p2);
    reb_move_to_com(sim); // move simulation to center of mass frame
    
    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_force* gd = rebx_load_force(rebx, "gas_damping_forces");
    rebx_add_force(rebx, gd);

    rebx_set_param_double(rebx, &sim->particles[1].ap, "gdf_damp_coeff", 1.); // damping coefficient set to 10
    rebx_set_param_double(rebx, &sim->particles[2].ap, "gdf_damp_coeff", 1.);

    double tmax = 1.e4*2*M_PI;  // set integration time for 10,000 years
    reb_integrate(sim, tmax);
    rebx_free(rebx);
    reb_free_simulation(sim);
}

void heartbeat(struct reb_simulation* sim){
    if(reb_output_check(sim, 5.e2*2*M_PI)){
        reb_move_to_hel(sim);
        reb_output_orbits(sim, "orbits.txt");
        reb_output_timing(sim, tmax);
    }
}
