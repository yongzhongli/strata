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
// numerical integration for the example in Fig. 3 of
// Ling, Jin, IEEE Microw. Guided Wave Lett., 2000. 
// When computing MoM integrals with numerical quadrature
// for overlapping or nearby mesh elements, one typically
// extracts the singularity of the Green's function by 
// expanding exp(-jkr)/r in a Taylor series. The 1/r term
// of the series is removed, and the remainder is used
// for the quadrature. The integral over the source element
// for the 1/r term is computed separately using an 
// analytical or semi-analytical approach (e.g., Gibson 2014).
// The same thing can be done with the MGF by extracting the
// quasistatic contribution in spectral domain and adding it
// back in spatial domain. When adding it back in spatial
// domain, each quasistatic term which has a 1/r behaviour
// is expanded in a Taylor series, and the 1/r term is
// removed, in anticipation of the user computing 1/r type
// integrals separately. The only difference is that in the
// quasistatic MGF, there are several terms with a 1/r
// behaviour, and some have multiplicative factors which
// depend on the layer materials etc. So, to add these 1/r
// terms during the triangle integration, the user needs to
// also know all these multiplicative factors. This example
// shows how to do this singularity extraction and then add
// back the 1/r terms along with their corresponding factors.

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
	std::cout << "TestSingularity()" << std::endl;
	std::cout << "===========================" << std::endl;

    std::string tech_file, out_file;
	
	if (argc < 2)
	{
		throw std::runtime_error("[ERROR] Layer file was not provided.");
	}
	else if (argc < 3)
	{
		tech_file = argv[1];
		out_file = "MGFdata_singularity.txt";
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
	// For this example, we'll sweep the observation point along the x axis from 10^{-4} wavelengths to 10 wavelengths away from the source point
	
	double x_src = 0.0, y_src = 0.0, z_src = 0.4e-3;
	double y_obs = 0.0, z_obs = 1.4e-3;
	
	int Nx = 500; // Number of points in the sweep
	double x_obs_min = std::abs(1.6e-4*lambda0);
	double x_obs_max = std::abs(1.6e1*lambda0);
		
	// We can use the Matlab-like linspace or logspace functions to create linearly- or logarithmically-spaced vectors points, provided via the Strata namespace
	std::vector<double> x_vec;
	strata::logspace(std::log10(x_obs_min), std::log10(x_obs_max), Nx, x_vec);


	// ====== Initialize the MGF class ======
	
	// This class stores all the settings we want to use
	MGF_settings s;
	
	// This class is the "engine" which will compute the MGF
	MGF mgf;

	// Tell the class that we want to use the numerical integration method
	s.method = MGF_INTEGRATE;

	// The singularity extraction requires either the quasistatic MGF or the homogeneous Green's function to be extracted in spectral domain and added back in spatial domain
	s.extract_quasistatic = true;
	s.extract_singularities = true;

	// Initialize the class with the given stackup and chosen settings
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


		// ====== Add singular terms ======

		// Add back the 1/r singular terms with the appropriate multiplicative factors.
		// In a real MoM context with a Galerkin scheme, r could be 0, so one would first integrate 1/r over the mesh element and then add the result.

		double rho = std::sqrt(std::pow(x_diff, 2) + std::pow(y_diff, 2));
		double r = std::sqrt(std::pow(rho, 2) + std::pow(z_obs - z_src, 2));

		for (int jj = 0; jj < 9; jj++)
			G_dyadic[jj] += mgf.GetSingularityFactor(jj)/r;
		G_phi += mgf.GetSingularityFactor(-1)/r;


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
		outfile << rho << " " << 
				std::abs(G_dyadic[0]) << " " << std::abs(G_dyadic[1]) << " " << std::abs(G_dyadic[2]) << " " << 
				std::abs(G_dyadic[3]) << " " << std::abs(G_dyadic[4]) << " " << std::abs(G_dyadic[5]) << " " << 
				std::abs(G_dyadic[6]) << " " << std::abs(G_dyadic[7]) << " " << std::abs(G_dyadic[8]) << " " << std::abs(G_phi) << std::endl;
		
	}
		
	return 0;

}






