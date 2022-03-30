/**
 * @file    gas_damping.c
 * @brief   Update orbits with prescribed timescales by directly changing orbital elements after each timestep, based on Dawson et al. 2016.
 * @author  Modified by Phoebe Sandhaus <pjs5535@psu.edu>; Original file by Dan Tamayo <tamayo.daniel@gmail.com>
 * 
 * @section     LICENSE
 * Copyright (c) 2015 Dan Tamayo, Hanno Rein
 *
 * This file is part of reboundx.
 *
 * reboundx is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * reboundx is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 * The section after the dollar signs gets built into the documentation by a script.  All lines must start with space * space like below.
 * Tables always must be preceded and followed by a blank line.  See http://docutils.sourceforge.net/docs/user/rst/quickstart.html for a primer on rst.
 * $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 *
 * $Orbit Modifications$       // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script).
 *
 * ======================= ===============================================
 * Authors                 D. Tamayo
 * Implementation Paper    `Tamayo, Rein, Shi and Hernandez, 2019 <https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.2885T/abstract>`_. 
 * Based on                `Lee & Peale 2002 <http://labs.adsabs.harvard.edu/adsabs/abs/2002ApJ...567..596L/>`_. 
 * C Example               :ref:`c_example_modify_orbits`
 * Python Example          `Migration.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/Migration.ipynb>`_,
 *                         `EccAndIncDamping.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/EccAndIncDamping.ipynb>`_.
 * ======================= ===============================================
 * 
 * This updates particles' positions and velocities between timesteps to achieve the desired changes to the osculating orbital elements (exponential growth/decay for a, e, inc, linear progression/regression for Omega/omega.
 * This nicely isolates changes to particular osculating elements, making it easier to interpret the resulting dynamics.  
 * One can also adjust the coupling parameter `p` between eccentricity and semimajor axis evolution, as well as whether the damping is done on Jacobi, barycentric or heliocentric elements.
 * Since this method changes osculating (i.e., two-body) elements, it can give unphysical results in highly perturbed systems.
 * 
 * **Effect Parameters**
 *
 * If p is not set, it defaults to 0.  If coordinates not set, defaults to using Jacobi coordinates.
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * coordinates (enum)           No          Type of elements to use for modification (Jacobi, barycentric or particle).
 *                                          See the examples for usage.
 * ============================ =========== ==================================================================
 *
 * **Particle Parameters**
 *
 * One can pick and choose which particles have which parameters set.  
 * For each particle, any unset parameter is ignored.
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * gdf_damp_coeff (double)          Yes         Damping coefficient d in Equation 16 from Dawson et al. 2016
 * ============================ =========== ==================================================================
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

static struct reb_vec3d rebx_calculate_gas_damping_forces(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* planet, struct reb_particle* star){
    struct rebx_extras* const rebx = sim->extras;
    struct reb_orbit o = reb_tools_particle_to_orbit(sim->G, *planet, *star);

    const double* const gdf_damp_coeff = rebx_get_param(rebx, planet->ap, "gdf_damp_coeff");

    // initialize damping timescales
    double tau_a = INFINITY;

    // initialize positions and velocities
    const double dvx = planet->vx - star->vx;
    const double dvy = planet->vy - star->vy;
    const double dvz = planet->vz - star->vz;
    const double dx = planet->x-star->x;
    const double dy = planet->y-star->y;
    const double dz = planet->z-star->z;
    const double r2 = dx*dx + dy*dy + dz*dz;

    // initial semimajor axis, eccentricity, and inclination
    const double a0 = o.a;
    const double e0 = o.e;
    const double inc0 = o.inc;
    const double starMass = star->m;
    const double planetMass = planet->m;

    // eccentricity and inclination timescales from Dawson+16 Eqn 16
    const double cs_coeff = 0.272125; // 1.29 km s^-1 in yr AU^-1
    double coeff;

    double v = sqrt(pow(e0, 2.)+pow(inc0, 2.))*pow(starMass, 1./2.)*pow(a0, -1./2.);
    double cs = cs_coeff/(2.*M_PI)*pow(a0, -1./4.);
    double v_over_cs = v/cs;

    if (v <= cs){
        coeff = 1.;
    }
    else {
        if (inc0 < cs/v) {
            coeff = pow(v_over_cs, 3.);
        }
        else {
            coeff = pow(v_over_cs, 4.);
        }
    }

    double tau_e = -0.003*(*gdf_damp_coeff)*pow(a0, 2.)*(starMass/planetMass)*2.*M_PI*coeff;
    double tau_inc = tau_e/2.;

    struct reb_vec3d a = {0};

    a.x =  dvx/(2.*tau_a);
    a.y =  dvy/(2.*tau_a);
    a.z =  dvz/(2.*tau_a);

    if (tau_e < INFINITY || tau_inc < INFINITY){
        const double vdotr = dx*dvx + dy*dvy + dz*dvz;
        const double prefac = 2*vdotr/r2/tau_e;
        a.x += prefac*dx;
        a.y += prefac*dy;
        a.z += prefac*dz + 2.*dvz/tau_inc;
    }
    return a;
}

void rebx_gas_damping_forces(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
    int* ptr = rebx_get_param(sim->extras, force->ap, "coordinates");
    enum REBX_COORDINATES coordinates = REBX_COORDINATES_JACOBI; // Default
    if (ptr != NULL){
        coordinates = *ptr;
    }
    const int back_reactions_inclusive = 1;
    const char* reference_name = "primary";
    rebx_com_force(sim, force, coordinates, back_reactions_inclusive, reference_name, rebx_calculate_gas_damping_forces, particles, N);
}