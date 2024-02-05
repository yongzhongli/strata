

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
#include <numeric>

#include "MGF.hpp"


int main(int argc, char** argv)
{
    bool print_log = true;
    std::cout << "===========================" << std::endl;
    std::cout << "TestInterp()" << std::endl;
    std::cout << "===========================" << std::endl;

    std::string tech_file, out_file1, out_file2;

    if (argc < 2)
    {
        throw std::runtime_error("[ERROR] Layer file was not provided.");
    }
    else if (argc < 3)
    {
        tech_file = argv[1];
        out_file1 = "plate_interp_real.txt";
        //out_file1 = "multiple_layers_interp_real_noAdaptive.txt";
        //out_file2 = "multiple_layers_interp_imag_noAdaptive.txt";
        out_file2 = "plate_interp_imag.txt";
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
    double f = 3e8;
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

    double x_src = 0.0, y_src = 0.0;
    double x_obs = 5.0e-3, y_obs = 0.0;

    int Nx = 500; // Number of points in the sweep
    int Nz_interpolated = 100; // Number of points to be interpolated

    // For this example, we'll sweep the observation point along the x axis from 10^{-4} wavelengths to 1 wavelengths away from the source point
    double x_obs_min = std::abs(1e-4*lambda0); // define the range of rho
    double x_obs_max = std::abs(1.0*lambda0);

    std::cout << "x_obs_min: " << x_obs_min << "(m)" << std::endl;
    std::cout << "x_obs_max: " << x_obs_max << "(m)" << std::endl;

    double lambda = lambda0/std::sqrt(12.5);

    double dis_threshold = 1e-8;
    double z_min = lm.layers.back().zmin + dis_threshold;
    double z_max = lm.layers.front().zmax - dis_threshold;

    double z_interp_min = z_min;
    double z_interp_max = z_max;

    // We can use the Matlab-like linspace or logspace functions to create linearly- or logarithmically-spaced vectors points, provided via the Strata namespace
    std::vector<double> x_vec;
    std::vector<double> z_src_vec;
    std::vector<double> z_obs_vec;
    strata::linspace(x_obs_min, x_obs_max, Nx, x_vec);
    strata::linspace(z_interp_min, z_interp_max, Nz_interpolated, z_src_vec);
    strata::linspace(z_interp_min, z_interp_max, Nz_interpolated, z_obs_vec);


    // ====== Create the interpolation grid ======

    // The interpolation grid consists of nodes along the z and rho axes.
    // When computing the MGF for arbitrary source and observation coordinates:
    // - along the z axis, Lagrange polynomial interpolation of a given order is performed,
    // - along the rho axis, Lagrange polynomial interpolation of a given order is performed.

    // In this example, Nz nodes are used to create the interpolation table
    // ====== Initialize the MGF class ======

    // This class stores all the settings we want to use
    MGF_settings s;

    std::vector<std::vector<double>> z_nodes_2D;
    std::vector<double> z_nodes_1D;

    lm.SetZnodes_interpGrid(f, z_nodes_2D, s.N_lambda); // Number of interpolation points
    s.interpolate_z = true;
    // Along rho, we'll pick 5 samples per wavelength based on the layer with maximum permittivity (12.5)
    double electrical_size_x = (x_obs_max - x_obs_min)/lambda;
    int N_rho = 10.0*electrical_size_x;


    // Strata comes with Matlab-like linspace and logspace functions for convenience
    std::vector<double> rho_nodes; strata::linspace(std::sqrt(std::pow(x_obs_min, 2) + std::pow(y_obs, 2)), std::sqrt(std::pow(x_obs_max, 2) + std::pow(y_obs, 2)), std::max(N_rho,10), rho_nodes);


    // Provide the layer manager with the z and rho nodes
    lm.ClearNodes_z();
    // Iterate over each sub-vector and insert its elements into z_nodes_1D
    for (const auto& subVec : z_nodes_2D)
    {
        z_nodes_1D.insert(z_nodes_1D.end(), subVec.begin(), subVec.end());
    }
    lm.InsertNodes_z(z_nodes_1D);

    lm.ClearNodes_rho();
    lm.InsertNodes_rho(rho_nodes);

    // This class is the "engine" which will compute the MGF
    MGF mgf;

    // Tell the class that we want to use the interpolation method.
    s.method = MGF_INTERPOLATE;
    s.filename = "./interpolation_table.txt";
    s.export_table = false;

    // For source and observation points close to each other, the interpolation can be inaccurate due to the singularity of the MGF.
    // So, instruct the MGF engine to extract the quasistatic part in the spectral domain, and then add it back in the spatial domain analytically.
    s.extract_quasistatic = true;

    // Polynomial interpolation order
    s.order = 2;
    s.order_z = 2;

    // Initialize the class with the given stackup and chosen settings.
    // The interpolation table will be generated at this point, so the initialization step could take a while.
    std::chrono::steady_clock::time_point begin1 = std::chrono::steady_clock::now();
    mgf.Initialize(f, lm, s);
    std::chrono::steady_clock::time_point end1 = std::chrono::steady_clock::now();
    static int Interpolation_initialization = std::chrono::duration_cast<std::chrono::milliseconds>(end1 - begin1).count();

    // ====== Compute the MGF ======

    // Create an output file where the MGF will be exported for post-processing
    std::ofstream outfile1(out_file1);
    std::ofstream outfile2(out_file2);

    // For post-processing, we'll store the frequency and positions along the z axis in the header
    outfile1 << "Frequency: " << f << " Hz" << std::endl;
    outfile2 << "Frequency: " << f << " Hz" << std::endl;

    // The lateral separation between the source and observation points (rho) and all the MGF components is tabulated for each observation point
    outfile1 << "\nz_prime z Gxx Gxy Gxz Gyx Gyy Gyz Gzx Gzy Gzz Gphi" << std::endl;
    outfile2 << "\nz_prime z Gxx Gxy Gxz Gyx Gyy Gyz Gzx Gzy Gzz Gphi" << std::endl;

    std::cout << "Number of z points: " << z_nodes_1D.size() << std::endl;
    std::cout << "Number of rho points: " << rho_nodes.size() << std::endl;

    std::chrono::steady_clock::time_point begin3 = std::chrono::steady_clock::now();
    mgf.adaptiveInterpolation(mgf);
    std::chrono::steady_clock::time_point end3 = std::chrono::steady_clock::now();
    static int AdaptiveInterpolation = std::chrono::duration_cast<std::chrono::milliseconds>(end3 - begin3).count();

    std::chrono::steady_clock::time_point begin2 = std::chrono::steady_clock::now();
    for (int ii = 0; ii < Nz_interpolated; ii++) {
        for (int jj = 0; jj < Nz_interpolated; jj++) {

            // The MGF class needs to know which layers we're working with, in order to perform some precomputations which may save time later on.
            // The working source and observation layers can be found with the FindLayer() method of the layer manager.
            // In a realistic MoM setting, if we are looping through the points on the mesh of an object, and we know in which layer that object resides, we can just set the source and observation layers once for each pair of source and observation objects.
            int i = lm.FindLayer(z_src_vec[ii]);
            int m = lm.FindLayer(z_obs_vec[jj]);
            mgf.SetLayers(i, m); // Source first, observation second


            // In the x and y directions, the MGF only depends on the separation between source and observation points, rather than the actual coordinates
            double x_diff = x_obs - x_src;
            double y_diff = y_obs - y_src;

            // The 3x3 dyadic MGF has 9 components which will be stored in a C++ standard array, in row major order.
            // In addition, there is a scalar component which is just a complex number.
            std::array<std::complex<double>, 9> G_dyadic;
            std::complex<double> G_phi;

            // Compute the MGF
            mgf.ComputeMGF(x_diff, y_diff, z_obs_vec[jj], z_src_vec[ii], G_dyadic, G_phi);

            // The dyadic MGF components are now stored in G_dyadic, while the scalar MGF is stored in G_phi.


            // ====== Export data to the output text file ======

            // Lateral separation
            outfile1 << z_src_vec[ii] << " " << z_obs_vec[jj] << " " <<
                    std::real(G_dyadic[0]) << " " << std::real(G_dyadic[1]) << " " << std::real(G_dyadic[2]) << " " <<
                    std::real(G_dyadic[3]) << " " << std::real(G_dyadic[4]) << " " << std::real(G_dyadic[5]) << " " <<
                    std::real(G_dyadic[6]) << " " << std::real(G_dyadic[7]) << " " << std::real(G_dyadic[8]) << " "
                    << std::real(G_phi) << std::endl;

            outfile2 << z_src_vec[ii] << " " << z_obs_vec[jj] << " " <<
                     std::imag(G_dyadic[0]) << " " << std::imag(G_dyadic[1]) << " " << std::imag(G_dyadic[2]) << " " <<
                     std::imag(G_dyadic[3]) << " " << std::imag(G_dyadic[4]) << " " << std::imag(G_dyadic[5]) << " " <<
                     std::imag(G_dyadic[6]) << " " << std::imag(G_dyadic[7]) << " " << std::imag(G_dyadic[8]) << " "
                     << std::imag(G_phi) << std::endl;

        }
    }
    std::chrono::steady_clock::time_point end2 = std::chrono::steady_clock::now();
    static int MGF_computation = std::chrono::duration_cast<std::chrono::milliseconds>(end2 - begin2).count();

    if (print_log)
    {
        for (int i = 0; i < mgf.lm.z_nodes.size(); i++)
            std::cout << "Number of z nodes for " << i << "th layer is: " << mgf.lm.z_nodes[i].size() << std::endl;
        std::cout << "Nz_interpolated: " << Nz_interpolated << std::endl;
        std::cout << "z_interp_min: " << z_interp_min << ", z_interp_max: " << z_interp_max << std::endl;
        std::cout << "Initialize the interpolation table: " << Interpolation_initialization << "(ms)" << std::endl;
        std::cout << "Update the interpolation table: " << AdaptiveInterpolation << "(ms)" << std::endl;
        std::cout << "Compute the MGF: " << MGF_computation << "(ms)" << std::endl;
    }


    return 0;

}






