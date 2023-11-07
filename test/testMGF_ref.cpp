// Author: Shashwat Sharma

// Copyright 2021 Shashwat Sharma and Piero Triverio

// This file is part of Strata.

// Strata is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// Strata is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with Strata.  If not, see <https://www.gnu.org/licenses/>.

// -----------------------------------------------------------

// The basic usage of Strata is demonstrated in this example.
// The Michalski-Zheng formulation-C MGF is computed using
// numerical integration for the example in Fig. 3 of
// Ling, Jin, IEEE Microw. Guided Wave Lett., 2000.

// -----------------------------------------------------------

#include <cmath>
#include <complex>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdexcept>

#include "MGF.hpp"

// THIS TEST FILE IS ONLY USED FOR GENERATING REFERENCE VALUES FOR NEW INTERPOLATION SCHEME OF MGF
int main(int argc, char** argv)
{

    std::cout << "===========================" << std::endl;
    std::cout << "TestMGF()" << std::endl;
    std::cout << "===========================" << std::endl;

    std::string tech_file, out_file1, out_file2;

    if (argc < 2)
    {
        throw std::runtime_error("[ERROR] Layer file was not provided.");
    }
    else if (argc < 3)
    {
        tech_file = argv[1];
        out_file1 = "MGFdata_ref_rho3_real.txt";
        out_file2 = "MGFdata_ref_rho3_imag.txt";
    }
    else
    {
        tech_file = argv[1];
        out_file1 = argv[2];
        out_file2 = argv[3];
    }


    // ====== Stackup ======

    // Create a layer management object and parse the layer file
    LayerManager lm;
    lm.ProcessTechFile(tech_file);

    // Set the analysis frequency and wave number
    double f = 30.0e9;
    double omega = 2.0*M_PI*f;

    // Some useful constants are provided via the Strata namespace
    double k0 = omega*std::sqrt(strata::eps0*strata::mu0);
    double lambda0 = 2.0*M_PI/k0;

    // Precompute frequency-dependent layer data
    lm.ProcessLayers(f);

    // Print layer data to the terminal for verification
    lm.PrintLayerData(NULL, true);


    // ====== Set up the source and observation points ======

    // Set the source and observation (obs) points
    // Uncomment / comment for the appropriate example

    // ------ Ling, Jin, 2000 ------

    // For this example, we'll sweep the observation point along the x axis from 10^{-4} wavelengths to 10 wavelengths away from the source point

    double x_src = 0.0, y_src = 0.0;
    double x_obs = 0.05e-3, y_obs = 0.0;

    double dis_threshold = 1e-8;
    double z_min = lm.layers.back().zmin + dis_threshold;
    double z_max = lm.layers.front().zmax - dis_threshold;

    //double z_interp_min = lm.layers.back().zmin;
    //double z_interp_max = lm.layers.front().zmax;

    double z_interp_min = z_min;
    double z_interp_max = z_max;

    // ====== Initialize the MGF class ======

    // In this example, we'll compute the MGF using straightforward numerical integration

    // This class stores all the settings we want to use
    MGF_settings s;


    // ------ Generate the points -----

    std::vector<double> z_nodes;
    lm.SetZnodes_interp(f, z_nodes, s.N_lambda); // Number of points for interpolation along z

    // This class is the "engine" which will compute the MGF
    MGF mgf;

    // Tell the class that we want to use the numerical integration method
    s.method = MGF_INTEGRATE;

    // Initialize the class with the given stackup and chosen settings
    mgf.Initialize(f, lm, s);


    // ====== Compute the MGF ======

    // Create an output file where the MGF will be exported for post-processing
    std::ofstream outfile1(out_file1);
    std::ofstream outfile2(out_file2);

    // For post-processing, we'll store the frequency and positions along the z axis in the header
    outfile1 << "Frequency: " << f << " Hz" << std::endl;
    outfile2 << "Frequency: " << f << " Hz" << std::endl;

    // The lateral separation between source and observation points (rho) and all the MGF components is tabulated for each observation point
    outfile1 << "\nz_prime z Gxx Gxy Gxz Gyx Gyy Gyz Gzx Gzy Gzz Gphi" << std::endl;
    outfile2 << "\nz_prime z Gxx Gxy Gxz Gyx Gyy Gyz Gzx Gzy Gzz Gphi" << std::endl;

    for (int ii = 0; ii < z_nodes.size(); ii++)
    {
        for (int jj = 0; jj < z_nodes.size(); jj++)
        {
            // The MGF class needs to know which layers we're working with, in order to perform some precomputations which may save time later on.
            // The working source and observation layers can be found with the FindLayer() method of the layer manager.
            // In a realistic MoM setting, if we are looping through the points on the mesh of an object, and we know in which layer that object resides, we can just set the source and observation layers once for each pair of source and observation objects.
            int i = lm.FindLayer(z_nodes[ii]);
            int m = lm.FindLayer(z_nodes[jj]);
            mgf.SetLayers(i, m); // Source first, observation second


            // In the x and y directions, the MGF only depends on the separation between source and observation points, rather than the actual coordinates
            double x_diff = x_obs - x_src;
            double y_diff = y_obs - y_src;

            // The 3x3 dyadic MGF has 9 components which will be stored in a C++ standard array, in row major order.
            // In addition, there is a scalar component which is just a complex number.
            std::array<std::complex<double>, 9> G_dyadic;
            std::complex<double> G_phi;

            // Compute the MGF
            mgf.ComputeMGF(x_diff, y_diff, z_nodes[jj], z_nodes[ii], G_dyadic, G_phi);

            // The dyadic MGF components are now stored in G_dyadic, while the scalar MGF is stored in G_phi.


            // ====== Optional modifications to the output ======

            // With reference to Michalski, Zheng, TAP 1990, 38 (3), equations (50)--(53), the MGF we computed above ** does include ** the cos and sin pre-factors. However, in literature, the MGF is often reported without these pre-factors. Therefore, for the purpose of this example, and to make direct comparisons to data from literature, those prefactors are cancelled out below. This section of code would not be needed in an actual MoM-type application.

            std::complex<double> zeta = std::atan2(y_diff, x_diff);
            std::complex<double> cos_term = std::cos(zeta);
            std::complex<double> sin_term = std::sin(zeta);

            if (std::abs(cos_term) > 0.0)
            {
                G_dyadic[2] /= cos_term;
                G_dyadic[6] /= cos_term;
            }
            if (std::abs(sin_term) > 0.0)
            {
                G_dyadic[5] /= sin_term;
                G_dyadic[7] /= sin_term;
            }


            // ====== Export data to the output text file ======

            // Lateral separation
            outfile1 << z_nodes[ii] << " " << z_nodes[jj] << " " <<
                     std::real(G_dyadic[0]) << " " << std::real(G_dyadic[1]) << " " << std::real(G_dyadic[2]) << " " <<
                     std::real(G_dyadic[3]) << " " << std::real(G_dyadic[4]) << " " << std::real(G_dyadic[5]) << " " <<
                     std::real(G_dyadic[6]) << " " << std::real(G_dyadic[7]) << " " << std::real(G_dyadic[8]) << " "
                     << std::real(G_phi) << std::endl;

            outfile2 << z_nodes[ii] << " " << z_nodes[jj] << " " <<
                     std::imag(G_dyadic[0]) << " " << std::imag(G_dyadic[1]) << " " << std::imag(G_dyadic[2]) << " " <<
                     std::imag(G_dyadic[3]) << " " << std::imag(G_dyadic[4]) << " " << std::imag(G_dyadic[5]) << " " <<
                     std::imag(G_dyadic[6]) << " " << std::imag(G_dyadic[7]) << " " << std::imag(G_dyadic[8]) << " "
                     << std::imag(G_phi) << std::endl;

        }
    }
    return 0;

}






