/***************************************************************************
 *
 * Authors:    Jose Luis Vilas, 					  jlvilas@cnb.csic.es
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include "resolution_pdb_bfactor.h"
#include <core/bilib/kernel.h>
#include "data/pdb.h"
#include <numeric>
#include <algorithm>
#include <iostream>


void ProgResBFactor::readParams()
{
	fn_pdb = getParam("--atmodel");
	fn_locres = getParam("--vol");
	sampling = getDoubleParam("--sampling");
	fnOut = getParam("-o");
}


void ProgResBFactor::defineParams()
{
	addUsageLine("The matching between a b-factor of an atomic model and "
			"the local resolution of a cryoEM map is analyzed.");
	addParamsLine("  --atmodel <vol_file=\"\">   		: Atomic model (pdb)");
	addParamsLine("  --vol <vol_file=\"\">				: Local resolution map");
	addParamsLine("  [--sampling <sampling=1>]			: Sampling Rate (A)");
	addParamsLine("  [--mean]			                : The resolution an bfactor per residue are averaged instead of computed the median");
	addParamsLine("  -o <output=\"amap.mrc\">			: Output of the algorithm");
}


void ProgResBFactor::analyzePDB()
{
	//Open the pdb file
	std::ifstream f2parse;
	f2parse.open(fn_pdb.c_str());

	double maxx=-1e-38, maxy, maxz, minx = 1e-38, miny, minz;
	maxy = maxx;
	maxz = maxx;
	miny = minx;
	minz = minx;

	numberOfAtoms = 0;

	int last_resi = 0;

	while (!f2parse.eof())
	{
		std::string line;
		getline(f2parse, line);

		// The type of record (line) is defined in the first 6 characters of the pdb
		std::string typeOfline = line.substr(0,4);
		//std::cout << typeOfline << std::endl;

		if ( (typeOfline == "ATOM") || (typeOfline == "HETA"))
		{
			numberOfAtoms++;
			double x = textToFloat(line.substr(30,8));
			double y = textToFloat(line.substr(38,8));
			double z = textToFloat(line.substr(46,8));

			// The furthest points along each axis are found to set the boxsize with the sampling
			if (x<minx)
				minx = x;
			if (y<miny)
			    miny = y;
			if (z<minz)
				minz = z;
			if (x>maxx)
				maxx = x;
			if (y>maxy)
				maxy = y;
			if (z>maxz)
				maxz = z;

			// Getting coordinates
			at_pos.x.push_back(x);
			at_pos.y.push_back(y);
			at_pos.z.push_back(z);

			int resi = (int) textToFloat(line.substr(23,5));

//			if (resi > last_resi)
//				last_resi = resi;

			at_pos.residue.push_back(resi);

			// Type of Atom
			std::string at = line.substr(13,2);

			// Getting the bfactor =8pi^2*u
			double bfactorRad = sqrt(textToFloat(line.substr(60,6))/(8*PI*PI));
			at_pos.b.push_back(bfactorRad);

			double rad = atomCovalentRadius(line.substr(13,2));

			at_pos.atomCovRad.push_back(rad);
		}
	}
}


template <typename T>
std::vector<size_t> ProgResBFactor::sort_indexes(const std::vector<T> &v)
{
  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values
  stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}


void ProgResBFactor::sweepByResidue(MultidimArray<int> &mask)
{
	Image<double> imgResVol;
	imgResVol.read(fn_locres);
	MultidimArray<double> resvol;
	resvol = imgResVol();

	size_t xdim, ydim, zdim, ndim;
	resvol.getDimensions(xdim, ydim, zdim, ndim);

	if (zdim == 1)
		zdim = ndim;

	std::vector<size_t> idx_residue;
	idx_residue = sort_indexes(at_pos.residue);

	std::cout << ".-.-.-.-.-.-..-.-.-.-.-.-.-.-..-.-" << std::endl;

	std::vector<double> resolution_per_residue(0);

	std::vector<double> resolution_to_estimate(0), resNumberList(0), bfactor_per_residue(0);

	// Setting the first residue before the loop with all residue
	size_t r = 0;
	size_t idx = idx_residue[r];
//	for (size_t i = 0; i<idx_residue.size(); ++i)
//	{
//		size_t idx = idx_residue[i];
//		std::cout << idx << std::endl;
//	}

	// Selecting the residue
	int resi = at_pos.residue[idx];
	int last_resi = resi;

	// Getting the atom position
	int k = round(at_pos.x[idx]/sampling) + floor(zdim/2);
	int i = round(at_pos.y[idx]/sampling) + floor(ydim/2);
	int j = round(at_pos.z[idx]/sampling) + floor(xdim/2);

	std::cout << i << " " << j << " " << k << std::endl;

	// Covalent Radius of the atom
	double covRad = at_pos.atomCovRad[idx];

	// Thermal displacement
	double bfactorRad = at_pos.b[idx];

	bfactor_per_residue.push_back(bfactorRad);

	// Total Displacement in voxels
	int totRad = round( (covRad + bfactorRad)/sampling );

	std::cout << totRad << std::endl;

	int dim;// = round(totRad/2 + 1);
	dim = totRad*totRad;
	for (size_t kk = 0; kk<totRad; ++kk)
	{
		size_t kk2 = kk * kk;
		for (size_t jj = 0; jj<totRad; ++jj)
		{
			size_t jj2kk2 = jj * jj + kk2;
			for (size_t ii = 0; ii<totRad; ++ii)
			{
				size_t dist2 = (ii)*(ii) + (jj)*(jj) + (kk)*(kk);



				if (dist2 <= dim)
				{
					std::cout << "entro" <<std::endl;
					A3D_ELEM(mask, k-kk, i-ii, j-jj) = 1;
					A3D_ELEM(mask, k-kk, i-ii, j+jj) = 1;
					A3D_ELEM(mask, k-kk, i+ii, j-jj) = 1;
					A3D_ELEM(mask, k-kk, i+ii, j+jj) = 1;
					A3D_ELEM(mask, k+kk, i-ii, j-jj) = 1;
					A3D_ELEM(mask, k+kk, i-ii, j+jj) = 1;
					A3D_ELEM(mask, k+kk, i+ii, j-jj) = 1;
					A3D_ELEM(mask, k+kk, i+ii, j+jj) = 1;

					resolution_to_estimate.push_back(A3D_ELEM(resvol, k-kk, i-ii, j-jj));
					resolution_to_estimate.push_back(A3D_ELEM(resvol, k-kk, i-ii, j+jj));
					resolution_to_estimate.push_back(A3D_ELEM(resvol, k-kk, i+ii, j-jj));
					resolution_to_estimate.push_back(A3D_ELEM(resvol, k-kk, i+ii, j+jj));
					resolution_to_estimate.push_back(A3D_ELEM(resvol, k+kk, i-ii, j-jj));
					resolution_to_estimate.push_back(A3D_ELEM(resvol, k+kk, i-ii, j+jj));
					resolution_to_estimate.push_back(A3D_ELEM(resvol, k+kk, i+ii, j-jj));
					resolution_to_estimate.push_back(A3D_ELEM(resvol, k+kk, i+ii, j+jj));
				}
			}
		}
	}


	MetaData md;
	size_t objId;


	std::cout << "numberOfAtoms" << numberOfAtoms << std::endl;


	for (size_t r=1; r<numberOfAtoms; ++r)
	{
		//std::cout << idx_residue[r] << std::endl;
		idx = idx_residue[r];

		// Selecting the residue
		resi = at_pos.residue[r];

		// Getting the atom position
		k = round(at_pos.x[idx]/sampling) + floor(zdim/2);
		i = round(at_pos.y[idx]/sampling) + floor(ydim/2);
		j = round(at_pos.z[idx]/sampling) + floor(xdim/2);
		//std::cout << k << " " << i << " " << j << "resi " << resi << std::endl;

		// Covalent Radius of the atom
		covRad = at_pos.atomCovRad[idx];

		// Thermal displacement
		bfactorRad = at_pos.b[idx];

		// Total Displacement in voxels
		totRad = round( (covRad + bfactorRad)/sampling );

		dim = totRad*totRad;


		if (resi != last_resi)
		{

			//std::cout << "New residue "<< resi << std::endl;
			std::sort(resolution_to_estimate.begin(), resolution_to_estimate.end());
			std::sort(bfactor_per_residue.begin(), bfactor_per_residue.end());

			double res_resi, bfactor_resi;
			//std::cout << resolution_to_estimate.size() << " " << size_t(resolution_to_estimate.size()*0.5) << std::endl;
			res_resi = resolution_to_estimate[size_t(resolution_to_estimate.size()*0.5)];
			bfactor_resi = bfactor_per_residue[size_t(bfactor_per_residue.size()*0.5)];

			resolution_per_residue.push_back(res_resi);
			resNumberList.push_back(res_resi);

			objId = md.addObject();
			md.setValue(MDL_BFACTOR, bfactor_resi, objId);
			md.setValue(MDL_RESIDUE, last_resi, objId);
			md.setValue(MDL_RESOLUTION_LOCAL_RESIDUE, last_resi, objId);
			std::cout << last_resi << "  " << res_resi << ";" << std::endl;

			last_resi = resi;

			resolution_to_estimate.clear();
			bfactor_per_residue.clear();
		}


//		std::cout << "--------------------------------------" << std::endl;

		for (size_t kk = 0; kk<totRad; ++kk)
		{
			size_t kk2 = kk * kk;
			for (size_t jj = 0; jj<totRad; ++jj)
			{
				size_t jj2kk2 = jj * jj + kk2;
				for (size_t ii = 0; ii<totRad; ++ii)
				{
					size_t dist2 = (ii)*(ii) + (jj)*(jj) + (kk)*(kk);
					if (dist2 <= dim)
					{
						A3D_ELEM(mask, k-kk, i-ii, j-jj) = 1;
						A3D_ELEM(mask, k-kk, i-ii, j+jj) = 1;
						A3D_ELEM(mask, k-kk, i+ii, j-jj) = 1;
						A3D_ELEM(mask, k-kk, i+ii, j+jj) = 1;
						A3D_ELEM(mask, k+kk, i-ii, j-jj) = 1;
						A3D_ELEM(mask, k+kk, i-ii, j+jj) = 1;
						A3D_ELEM(mask, k+kk, i+ii, j-jj) = 1;
						A3D_ELEM(mask, k+kk, i+ii, j+jj) = 1;

						resolution_to_estimate.push_back(A3D_ELEM(resvol, k-kk, i-ii, j-jj));
						resolution_to_estimate.push_back(A3D_ELEM(resvol, k-kk, i-ii, j+jj));
						resolution_to_estimate.push_back(A3D_ELEM(resvol, k-kk, i+ii, j-jj));
						resolution_to_estimate.push_back(A3D_ELEM(resvol, k-kk, i+ii, j+jj));
						resolution_to_estimate.push_back(A3D_ELEM(resvol, k+kk, i-ii, j-jj));
						resolution_to_estimate.push_back(A3D_ELEM(resvol, k+kk, i-ii, j+jj));
						resolution_to_estimate.push_back(A3D_ELEM(resvol, k+kk, i+ii, j-jj));
						resolution_to_estimate.push_back(A3D_ELEM(resvol, k+kk, i+ii, j+jj));
					}
				}
			}
		}

		bfactor_per_residue.push_back(bfactorRad);

	}

	md.write("bfactor_resolution.xmd");

	Image<int> imMask;
	imMask() = mask;
	imMask.write("mascara.mrc");

}




template <typename T>
void ProgResBFactor::maskFromPDBData(struct pdbdata &coord, MultidimArray<T> &mask)
{
	MultidimArray<double> resVol;
	Image<double> resVolImg;
	resVolImg.read(fn_locres);
	resVol = resVolImg();

	mask.resizeNoCopy(resVol);
	mask.initZeros();

	size_t xdim, ydim, zdim, ndim;
	resVol.getDimensions(xdim, ydim, zdim, ndim);

	std::cout << xdim << " " << ydim << " " << ndim << std::endl;

//	double Xorig, Yorig, Zorig;
//
//	Xorig = xdim/2.0;
//	Yorig = ydim/2.0;
//	Zorig = zdim/2.0;
//
//	std::cout << Xorig << "  " << Yorig << "  " << Zorig << std::endl;
//
//	for (size_t at = 0; at<numberOfAtoms; ++at)
//	{
//		int k = round(coord.x[at]/sampling);
//		int i = round(coord.y[at]/sampling);
//		int j = round(coord.z[at]/sampling);
//
//		std::cout << k << " " << i << " " << j << std::endl;
//
//		A3D_ELEM(mask, k, i, j) = 1;
//	}
}


void ProgResBFactor::run()
{
	MultidimArray<int> mask;

	analyzePDB();

	std::cout << "The pdb was parsed" << std::endl;

	maskFromPDBData(at_pos, mask);
	std::cout << "The pdb was parsed" << std::endl;

	sweepByResidue(mask);
	
//	Image<int> imgsave;
//	imgsave() = mask;
//	imgsave.write("mascara.mrc");

}
