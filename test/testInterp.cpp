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

// The Michalski-Zheng formulation-C MGF is computed using
// an interpolation table for the example in Fig. 3 of
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


int main(int argc, char** argv)
{

    std::cout << "===========================" << std::endl;
    std::cout << "TestInterp()" << std::endl;
    std::cout << "===========================" << std::endl;

    std::string tech_file, out_file;

    if (argc < 2)
    {
        throw std::runtime_error("[ERROR] Layer file was not provided.");
    }
    else if (argc < 3)
    {
        tech_file = argv[1];
        out_file = "MGFdata_interp.txt";
    }
    else
    {
        tech_file = argv[1];
        out_file = argv[2];
    }


    // ====== Stackup ======

    // Create a layer management object and parse the layer file
    LayerManager lm;
    lm.ProcessTechFile(tech_file);

    // Set the analysis frequency and wave number
    double f = 1.0e9;
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
    // For this example, we'll sweep the observation point along the x axis from 10^{-4} wavelengths to 10 wavelengths away from the source point

    double x_src = 0.0, y_src = 0.0, z_src = 0.2e-3;
    double y_obs = 0.1e-3, z_obs = 0.0e-3;

    int Nx = 500; // Number of points in the sweep
    double x_obs_min = std::abs(1.6e-4*lambda0);
    double x_obs_max = std::abs(1.6*lambda0);

    double z_min = 0.1e-3;
    double z_max = 0.7e-3;

    // We can use the Matlab-like linspace or logspace functions to create linearly- or logarithmically-spaced vectors points, provided via the Strata namespace
    std::vector<double> x_vec;
    strata::logspace(std::log10(x_obs_min), std::log10(x_obs_max), Nx, x_vec);


    // ====== Create the interpolation grid ======

    // The interpolation grid consists of nodes along the z and rho axes.
    // When computing the MGF for arbitrary source and observation coordinates:
    // - along the z axis, the nearest z-node is chosen (0-order interpolation), and,
    // - along the rho axis, Lagrange polynomial interpolation of a given order is performed.

    // In this example, only two points are needed along z
    //std::vector<double> z_nodes = {z_src,z_obs};

    // Along rho, we'll pick 5 samples per wavelength based on the layer with maximum permittivity (12.5)
    //double lambda = lambda0/std::sqrt(9.8);
    //double electrical_size = (x_obs_max - x_obs_min)/lambda;
    //int N_rho = 5.0*electrical_size;

    // Strata comes with Matlab-like linspace and logspace functions for convenience
    //std::vector<double> rho_nodes;
    //strata::linspace(std::sqrt(std::pow(x_obs_min, 2) + std::pow(y_obs, 2)), std::sqrt(std::pow(x_obs_max, 2) + std::pow(y_obs, 2)), N_rho, rho_nodes);

    // Provide the layer manager with the z and rho nodes

    // z-nodes and their corresponding layer indices
    std::vector<double> z_nodes = {
        7.8e-05, 8.134e-05, 8.468e-05, 8.802e-05,
        3.9962e-05, 4.2888e-05, 4.5814e-05, 4.874e-05,
        5.1666e-05, 5.4592e-05, 5.7518e-05, 6.0444e-05,
        6.337e-05, 6.6296e-05, 6.9222e-05, 7.2148e-05,
        7.5074e-05, 7.8e-05
    };

    std::vector<int> z_layers = {
        0, 0, 0, 0,   // First 4 z-nodes are in layer 0
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 // Remaining 14 z-nodes are in layer 1
    };

    // rho-nodes as provided
    std::vector<double> rho_nodes = {
        0.0,
        2.4902439487712573e-05,
        4.9804878975425146e-05,
        7.470731846313773e-05,
        9.960975795085029e-05,
        0.00012451219743856286,
        0.00014941463692627545,
        0.00017431707641398802,
        0.00019921951590170058,
        0.00022412195538941315,
        0.0002490243948771257,
        0.0002739268343648383,
        0.0002988292738525509,
        0.00032373171334026344,
        0.00034863415282797603,
        0.00037353659231568857,
        0.00039843903180340116,
        0.00042334147129111376,
        0.0004482439107788263,
        0.0004731463502665389,
        0.0004980487897542514,
        0.0005229512292419641,
        0.0005478536687296766,
        0.0005727561082173892,
        0.0005976585477051018,
        0.0006225609871928143,
        0.0006474634266805269,
        0.0006723658661682394,
        0.0006972683056559521,
        0.0007221707451436646,
        0.0007470731846313771,
        0.0007719756241190898,
        0.0007968780636068023,
        0.0008217805030945149,
        0.0008466829425822275,
        0.0008715853820699401,
        0.0008964878215576526,
        0.0009213902610453652,
        0.0009462927005330778,
        0.0009711951400207903
    };


    lm.ClearNodes_z();
    std::vector<double> z_node{1};
    for (int ii = 0; ii < z_nodes.size(); ii++) {
        z_node[0] = z_nodes[ii];
        lm.InsertNodes_z(z_node, z_layers[ii]);
    }
    lm.InsertNodes_z(z_nodes);

    lm.ClearNodes_rho();
    lm.InsertNodes_rho(rho_nodes);


    // ====== Initialize the MGF class ======

    // This class stores all the settings we want to use
    MGF_settings s;

    // This class is the "engine" which will compute the MGF
    MGF mgf;

    // Tell the class that we want to use the interpolation method.
    s.method = MGF_INTERPOLATE;

    // For source and observation points close to each other, the interpolation can be inaccurate due to the singularity of the MGF.
    // So, instruct the MGF engine to extract the quasistatic part in the spectral domain, and then add it back in the spatial domain analytically.
    s.extract_quasistatic = true;

    // Polynomial interpolation order
    s.order = 3;

    // Initialize the class with the given stackup and chosen settings.
    // The interpolation table will be generated at this point, so the initialization step could take a while.
    mgf.Initialize(f, lm, s);


    // ====== Compute the MGF ======

    // Create an output file where the MGF will be exported for post-processing
    std::ofstream outfile(out_file);

    // For post-processing, we'll store the frequency and positions along the z axis in the header
    outfile << "Frequency: " << f << " Hz" << std::endl;
    outfile << "z_src: " << z_src << " m" << std::endl;
    outfile << "z_obs: " << z_obs << " m" << std::endl;

    // The lateral separation between source and observation points (rho) and all the MGF components is tabulated for each observation point
    outfile << "\nrho Gxx Gxy Gxz Gyx Gyy Gyz Gzx Gzy Gzz Gphi" << std::endl;

    // The MGF class needs to know which layers we're working with, in order to perform some precomputations which may save time later on.
    // The working source and observation layers can be found with the FindLayer() method of the layer manager.
    // In a realistic MoM setting, if we are looping through the points on the mesh of an object, and we know in which layer that object resides, we can just set the source and observation layers once for each pair of source and observation objects.
    int i = lm.FindLayer(z_src);
    int m = lm.FindLayer(z_obs);
    mgf.SetLayers(i, m); // Source first, observation second

    for (int ii = 0; ii < Nx; ii++)
    {

        double x_obs = x_vec[ii];

        // In the x and y directions, the MGF only depends on the separation between source and observation points, rather than the actual coordinates
        double x_diff = x_obs - x_src;
        double y_diff = y_obs - y_src;

        // The 3x3 dyadic MGF has 9 components which will be stored in a C++ standard array, in row major order.
        // In addition, there is a scalar component which is just a complex number.
        std::array<std::complex<double>, 9> G_dyadic;
        std::complex<double> G_phi;

        // Compute the MGF
        mgf.ComputeMGF(x_diff, y_diff, z_obs, z_src, G_dyadic, G_phi);

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
        double rho = std::sqrt( std::pow(x_diff, 2) + std::pow(y_diff, 2) );
        outfile << rho << " " <<
                std::abs(G_dyadic[0]) << " " << std::abs(G_dyadic[1]) << " " << std::abs(G_dyadic[2]) << " " <<
                std::abs(G_dyadic[3]) << " " << std::abs(G_dyadic[4]) << " " << std::abs(G_dyadic[5]) << " " <<
                std::abs(G_dyadic[6]) << " " << std::abs(G_dyadic[7]) << " " << std::abs(G_dyadic[8]) << " " << std::abs(G_phi) << std::endl;

    }

    return 0;

}