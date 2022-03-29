/**
 * Planetary evolution
 * code unit: AU, 𝑀⊙, yr/2𝜋, G = 1
 * fordj.txt: Columns are a, e, i, peri, node, M, mass 
 * with all angles in degrees, all elements in heliocentric coordinates, and masses relative to solar
 * The central body is 1 solar mass.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include "rebound.h"
#include "reboundx.h"

#define PLANET_MAX 1000

void heartbeat(struct reb_simulation* sim);
double tmax;

int main(int argc, char* argv[]){

    // check if the problem.c file is starting within Huang2017 or SolarSystem
    // then change initial filename based on which
    char path[200];
    getcwd(path, 200);
    char *giantPlanetConfig;
    char* filename = "default.txt";
    giantPlanetConfig = strstr(path, "SolarSystem"); 

    if (giantPlanetConfig) {
        filename = "initialPlanetProperties_SolarSystem.txt";
    }
    else {
        filename = "initialPlanetProperties_Huang2017.txt";
    }    

    // check if existing instances of simulation exist
    // then restart simulation using the current_orbits.txt file
    FILE* initialFile = fopen(filename, "r");
    FILE* currentOrbits = fopen("current_orbits.txt", "r");
    if (currentOrbits == NULL) {
        
        /* import from initial   a         e        i      peri     node       M        mass  */
        double a[PLANET_MAX]; // Semi-major axis (au)
        double e[PLANET_MAX]; // Eccentricity
        double inc[PLANET_MAX]; // Inc (deg)
        double omega[PLANET_MAX]; // Argument of peri (deg)
        double Omega[PLANET_MAX]; // Longitude of ascending node (deg)
        double M[PLANET_MAX]; // Mean longitude (deg)
        double mass[PLANET_MAX]; // Mass in solar

        memset(a, 0, sizeof(a));
        memset(e, 0, sizeof(e));
        memset(inc, 0, sizeof(inc));
        memset(omega, 0, sizeof(omega));
        memset(Omega, 0, sizeof(Omega));
        memset(M, 0, sizeof(M));
        memset(mass, 0, sizeof(mass));
        
        int n_embryos = 0;
        while (!feof(initialFile) && (n_embryos < PLANET_MAX)){
            // a     e    i   peri   node   M    mass
            fscanf(initialFile, "%lf %lf %lf %lf %lf %lf %lf", 
                &(a[n_embryos]), &(e[n_embryos]), 
                &(inc[n_embryos]), &(omega[n_embryos]), &(Omega[n_embryos]), &(M[n_embryos]), 
                &(mass[n_embryos]));
            n_embryos++;
        }
        printf("Number of planetary embryos: %d\n", n_embryos);
        fclose(initialFile);

        /* start the rebound simulation here */
        struct reb_simulation* sim = reb_create_simulation();
        
        if (!strcmp(filename, "initialPlanetProperties_Huang2017.txt")) {
            sim->integrator = REB_INTEGRATOR_IAS15;  // Use the same integrator as Huang 2017
        }
        else {
            sim->integrator = REB_INTEGRATOR_MERCURIUS;
        }
        sim->ri_mercurius.hillfac = 3;
        sim->dt = 0.00137*2.*M_PI;
        
        
        sim->collision            = REB_COLLISION_DIRECT;
        sim->collision_resolve    = reb_collision_resolve_merge;       // Choose merger collision routine.
        
        sim->heartbeat = heartbeat;

        // add the host stars
        struct reb_particle star = {0};
        star.m      = 1.;
        star.r      = 0.005;
        star.hash   = reb_hash("CentralStar");
        reb_add(sim, star);
        
        // add planets
        for (int i = 0; i < n_embryos; i++){
            double this_f = reb_tools_M_to_f(e[i], M[i]*(M_PI/180.));
            double this_inc = inc[i]*(M_PI/180.);
            double this_Omega = Omega[i]*(M_PI/180.);
            double this_omega = omega[i]*(M_PI/180.);
            struct reb_particle p = reb_tools_orbit_to_particle(sim->G, star, mass[i], a[i], e[i], this_inc, 
                                                                this_Omega, this_omega, this_f);
            double density = 1.68372e6; // 1 g cm^-3 in Msun AU^-3
            p.r = pow(3.*mass[i]/(4.*M_PI*density),1./3.); // calculate collision r
            p.lastcollision = 0;
            // assign planet number to each body
            p.hash = i;
            reb_add(sim, p);
        }
        
        reb_move_to_com(sim);

        struct rebx_extras* rebx = rebx_attach(sim);
        struct rebx_force* gd = rebx_load_force(rebx, "gas_damping_forces");
        rebx_add_force(rebx, gd);

        for (int i = 0; i < sim->N; i++){
            rebx_set_param_double(rebx, &sim->particles[i+1].ap, "gdf_damp_coeff", 100.); // damping coefficient set to 100
        }

        system("rm -v orbits.txt");
        system("rm -v collisions.txt");

        tmax = 1.e6*2*M_PI;
        reb_integrate(sim, tmax);
        rebx_free(rebx);
        reb_free_simulation(sim);

        exit(1);
    } 
    else {
        printf("Restarting integration...");
        /* import from currentOrbits   t     a       e      i     node     peri      M     P     f    mass    planetNum */
        double t[PLANET_MAX]; // current simulation time (yr/2𝜋)
        double a[PLANET_MAX]; // Semi-major axis (au)
        double e[PLANET_MAX]; // Eccentricity
        double inc[PLANET_MAX]; // Inc (radians)
        double Omega[PLANET_MAX]; // Longitude of ascending node (radians) 
        double omega[PLANET_MAX]; // Argument of peri (radians)
        double M[PLANET_MAX]; // Mean longitude (radians)
        double P[PLANET_MAX]; // orbital period (yr/2𝜋)
        double f[PLANET_MAX]; // true anomaly (radians)
        double mass[PLANET_MAX]; // Mass in solar
        int planetNum[PLANET_MAX]; // corresponding planet number

        memset(t, 0, sizeof(t));
        memset(a, 0, sizeof(a));
        memset(e, 0, sizeof(e));
        memset(inc, 0, sizeof(inc));
        memset(Omega, 0, sizeof(Omega));
        memset(omega, 0, sizeof(omega));
        memset(M, 0, sizeof(M));
        memset(P, 0, sizeof(P));
        memset(f, 0, sizeof(f));
        memset(mass, 0, sizeof(mass));
        memset(planetNum, 0, sizeof(planetNum));

        int n_embryos = 0;
        while (!feof(currentOrbits) && (n_embryos < PLANET_MAX)){
            // t    a     e    i   node     peri     M     P     f   mass     planetNum
            fscanf(currentOrbits, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d", 
                &(t[n_embryos]), &(a[n_embryos]), &(e[n_embryos]), 
                &(inc[n_embryos]), &(Omega[n_embryos]), &(omega[n_embryos]), &(M[n_embryos]), 
                &(P[n_embryos]), &(f[n_embryos]), &(mass[n_embryos]), &(planetNum[n_embryos]));
            n_embryos++;
        }
        printf("Number of planetary embryos: %d\n", n_embryos);
        fclose(currentOrbits);

        /* start the rebound simulation here */
        struct reb_simulation* sim = reb_create_simulation();
        
        if (!strcmp(filename, "initialPlanetProperties_Huang2017.txt")) {
            sim->integrator = REB_INTEGRATOR_IAS15;  // Use the same integrator as Huang 2017
        }
        else {
            sim->integrator = REB_INTEGRATOR_MERCURIUS;
        }

        sim->ri_mercurius.hillfac = 3;
        sim->dt = 0.00137*2.*M_PI;
        sim->t = t[0];
        
        sim->collision            = REB_COLLISION_DIRECT;
        sim->collision_resolve    = reb_collision_resolve_merge;       // Choose merger collision routine.
        
        sim->heartbeat = heartbeat;

        // add the host stars
        struct reb_particle star = {0};
        star.m      = 1.;
        star.r      = 0.005;
        star.hash   = reb_hash("CentralStar");
        reb_add(sim, star);
        
        // add planets
        for (int i = 0; i < n_embryos; i++){
            struct reb_particle p = reb_tools_orbit_to_particle(sim->G, star, mass[i], a[i], e[i], inc[i], Omega[i], omega[i], f[i]);
            double density = 1.68372e6; // 1 g cm^-3 in Msun AU^-3
            p.r = pow(3.*mass[i]/(4.*M_PI*density),1./3.); // calculate collision r
            p.lastcollision = 0;
            // assign planet number to each body
            p.hash = planetNum[i];
            reb_add(sim, p);
        }
        
        reb_move_to_com(sim);

        struct rebx_extras* rebx = rebx_attach(sim);
        struct rebx_force* gd = rebx_load_force(rebx, "gas_damping_forces");
        rebx_add_force(rebx, gd);

        for (int i = 0; i < sim->N; i++){
            rebx_set_param_double(rebx, &sim->particles[i+1].ap, "damp_coeff", 100.); // damping coefficient set to 100
        }

        tmax = 1.e6*2*M_PI;
        reb_integrate(sim, tmax);
        rebx_free(rebx);
        reb_free_simulation(sim);

        exit(1);
    }

    return 0;
}

void heartbeat(struct reb_simulation* sim){
    if(reb_output_check(sim, 2.74*2*M_PI)){
        reb_move_to_hel(sim);
        reb_output_orbits(sim, "orbits.txt");
        system("rm -v current_orbits.txt");
        reb_output_orbits(sim, "current_orbits.txt");
        system("truncate -s -1 current_orbits.txt");
    }
}
