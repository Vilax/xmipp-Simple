/***************************************************************************
 *
 * Authors:    Jose Luis Vilas (joseluis.vilas-prieto@yale.edu)
 *                             or (jlvilas@cnb.csic.es)
 *              Hemant. D. Tagare (hemant.tagare@yale.edu)
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

#ifndef _PROG_RES_DIR
#define _PROG_RES_DIR

#include <core/xmipp_program.h>
#include <core/xmipp_fftw.h>
#include <core/metadata_extension.h>
#include <data/monogenic_signal.h>


class ProgFSO : public XmippProgram
{
public:
        // Filenames
        FileName fnhalf1, fnhalf2, fnmask, fn_root, fn_3dfsc, fn_fscmd_folder, fn_ani, fnOut;
    
        // Double Params
        double sampling, ang_con, thrs;

		// Int params
		int Nthreads;
		size_t xvoldim, yvoldim, zvoldim;

		//long params
		//number of Fourier elements
		long Ncomps;
        
        // Bool params
        bool test, do_3dfsc_filter;

        //Matrix2d for the projection angles
        Matrix2D<double> angles;

		//Frequency vectors and frequency map
		Matrix1D<double> freq_fourier_x;
		Matrix1D<double> freq_fourier_y;
		Matrix1D<double> freq_fourier_z;
		MultidimArray<double> fx, fy, fz, threeD_FSC, normalizationMap;
		MultidimArray< double > freqMap;

		//Half maps
		MultidimArray< std::complex< double > > FT1, FT1_vec;
		MultidimArray< std::complex< double > > FT2, FT2_vec;

		//Access indices
		MultidimArray<long> freqElems, cumpos, freqidx, arr2indx;


public:
        // Defining the params and help of the algorithm
        void defineParams();

        // It reads the input parameters
        void readParams();

        // The output is a multidim array that define INVERSE of the value of the frequency in Fourier Space
        // To do that, it makes use of myfftV to detemine the size of the output map (myfftV and
        // the output will have the same size), and the vectors freq_fourier_x, freq_fourier_y, 
        // and freq_fourier_z that defines the frequencies along the 3 axis. The output will be
        // sqrt(freq_fourier_x^2+freq_fourier_y^2+freq_fourier_x^2)
        void defineFrequencies(const MultidimArray< std::complex<double> > &myfftV,
    		                                    const MultidimArray<double> &inputVol);

        // Estimates the directional FSC between two half maps FT1 and FT2 (in Fourier Space)
        // requires the sampling rate, and the frequency vectors, 
        void fscDir(MultidimArray< std::complex< double > > & FT1,
            	 MultidimArray< std::complex< double > > & FT2,
                 double sampling_rate,
				 Matrix1D<double> &freq_fourier_x,
				 Matrix1D<double> &freq_fourier_y,
				 Matrix1D<double> &freq_fourier_z,
				 MultidimArray< double >& freqMap,
                 MultidimArray< double >& freq,
                 MultidimArray< double >& frc,
    			 double maxFreq, int m1sizeX, int m1sizeY, int m1sizeZ,
				 double rot, double tilt, double ang_con, double &dres, double &thrs);

        // Estimates the global FSC between two half maps FT1 and FT2 (in Fourier Space)
        void fscGlobal(double sampling_rate,
                 MultidimArray< double >& freq,
                 MultidimArray< double >& frc,
    			 double maxFreq, int m1sizeX, int m1sizeY, int m1sizeZ, MetaData &mdRes,
				 double &fscFreq, double &thrs, double &resInterp);

        // Defines a map (sphere) in Fourier space after ffsshift (complete Fourier space) for which
        // each radius defines the frequency of that point in Fourier Space
        void createfrequencySphere(MultidimArray<double> &sphere,
    		Matrix1D<double> &freq_fourier_x,
			 Matrix1D<double> &freq_fourier_y,
			 Matrix1D<double> &freq_fourier_z);

        void crossValues(Matrix2D<double> &indexesFourier, double &rot, double &tilt, double &angCon,
			 MultidimArray<std::complex<double>> &f1, MultidimArray<std::complex<double>> &f2,
			 std::complex<double> &f1_mean, std::complex<double> &f2_mean);

        void weights(double freq, Matrix2D<double> &indexesFourier, Matrix2D<int> &indexesFourier2, double &rot, double &tilt, double &angCon,
			 MultidimArray<std::complex<double>> &f1, MultidimArray<std::complex<double>> &f2,
			 MultidimArray<std::complex<double>> &FT1, MultidimArray<std::complex<double>> &FT2,
			 Matrix1D<double> &freq_fourier_x,
			 Matrix1D<double> &freq_fourier_y,
			 Matrix1D<double> &freq_fourier_z,
			 double &cross);


        void fscShell(MultidimArray< std::complex< double > > & FT1,
    		 MultidimArray< std::complex< double > > & FT2,
			 Matrix1D<double> &freq_fourier_x,
			 Matrix1D<double> &freq_fourier_y,
			 Matrix1D<double> &freq_fourier_z,
			 MultidimArray< double >& freqMap,
			 int m1sizeX, Matrix2D<double> &indexesFourier, Matrix2D<int> &indexesFourier2, double &cutoff,
			 MultidimArray<std::complex<double>> &f1, MultidimArray<std::complex<double>> &f2);
        
        // Defines a Matrix2D with coordinates Rot and tilt achieving a uniform coverage of the
        // projection sphere. Bool alot = True, implies a dense converage
        void generateDirections(Matrix2D<double> &angles, bool alot);

        void interpolationCoarse(MultidimArray< double > fsc,
    		const Matrix2D<double> &angles,
			Matrix1D<double> &freq_fourier_x,
			Matrix1D<double> &freq_fourier_y,
			Matrix1D<double> &freq_fourier_z,
    		MultidimArray<double> &threeD_FSC,
			MultidimArray<double> &counterMap,
			MultidimArray< double >& freqMap,
			MultidimArray< double >& freq,
			double maxFreq, int m1sizeX, int m1sizeY, int m1sizeZ,
			double rot, double tilt, double ang_con);

        void anistropyParameter(const MultidimArray<double> FSC,
    		MultidimArray<double> &directionAnisotropy, size_t dirnumber,
			MultidimArray<double> &aniParam, double thrs);
    
        void anistropyParameterSimple(const MultidimArray<double> FSC,
			MultidimArray<double> &aniParam, double thrs);
        
        void prepareData(MultidimArray<double> &half1, MultidimArray<double> &half2, bool test);

        void saveFSCToMetadata(MetaData &mdRes,
    		const MultidimArray<double> &freq,
			const MultidimArray<double> &FSC, FileName &fnmd);

        void saveAnisotropyToMetadata(MetaData &mdAnisotropy,
    		const MultidimArray<double> &freq,
			const MultidimArray<double> &anisotropy);

        void directionalFilter(MultidimArray<std::complex<double>> &FThalf1,
    		MultidimArray<double> &threeDfsc, MultidimArray<double> &filteredMap, 
            int m1sizeX, int m1sizeY, int m1sizeZ);

		void directionalFilter_reading(MultidimArray<std::complex<double>> &FThalf1,
    		MultidimArray<double> &threeDfsc, MultidimArray<double> &filteredMap, 
            int m1sizeX, int m1sizeY, int m1sizeZ);
        
        void resolutionDistribution(MultidimArray<double> &resDirFSC, FileName &fn);

        void getCompleteFourier(MultidimArray<double> &V, MultidimArray<double> &newV,
    		int m1sizeX, int m1sizeY, int m1sizeZ);

        void createFullFourier(MultidimArray<double> &fourierHalf, FileName &fnMap,
    		int m1sizeX, int m1sizeY, int m1sizeZ);

        void run();

		void run_fast();

		void run_old();

        void estimateSSNR(MultidimArray<double> &half1, MultidimArray<double> &half2,
		                int m1sizeX, int m1sizeY, int m1sizeZ);

        void directionalSSNR(MultidimArray< std::complex< double > > & FT1,
					 MultidimArray< std::complex< double > > & FT2, double sampling_rate,
					 Matrix1D<double> &freq_fourier_x,
					 Matrix1D<double> &freq_fourier_y,
					 Matrix1D<double> &freq_fourier_z,
					 MultidimArray< double >& freqMap, MultidimArray< double >& sig,
					 MultidimArray< double >& noi,
					 double maxFreq, int m1sizeX, int m1sizeY, int m1sizeZ,
					 double rot, double tilt, double ang_con, size_t dire);

        void getErrorCurves(int &m1sizeX, int &m1sizeY, int &m1sizeZ,
	                        Matrix1D<double> &freq_fourier_x,
	                        Matrix1D<double> &freq_fourier_y,
	                        Matrix1D<double> &freq_fourier_z, 
                            MultidimArray<double> &freqMap, size_t Nrealization, double thrs);
		
		void noiseStatisticsInMask(MultidimArray<double> &map, MultidimArray<double> &mask,
 									double &mean, double &stdev);

		void createNoisyMap(MultidimArray<double> &map, double mean, double stddev);

		void createNoisyFringePattern(MultidimArray<double> &map, MultidimArray<double> &noise, 
									  MultidimArray<double> &mask,
										double sqrtpowernoise, double wavelength);

		void findBestConeAngle(Matrix2D<int> &fscShell, double resolutionfsc);

		// void shellValue(Matrix2D<double> &indexesFourier, double &rot, double &tilt,
		// 	MultidimArray<std::complex<double>> &f1, MultidimArray<std::complex<double>> &f2,
		// 	std::complex<double> &f_coeff_1, std::complex<double>  &f_coeff_2);

		void shellValue(double freq, double &rot, double &tilt,
			MultidimArray<std::complex<double>> &FT1, MultidimArray<std::complex<double>> &FT2,
			Matrix1D<double> &freq_fourier_x,
			Matrix1D<double> &freq_fourier_y,
			Matrix1D<double> &freq_fourier_z,
			std::complex<double> &f_coeff_1, std::complex<double>  &f_coeff_2);

		void findIndexinVector(double freq, double x_dir, size_t &idx,
								Matrix1D<double> &freq_fourier);

		void arrangeFSC_and_fscGlobal(double sampling_rate, 
				double &fscFreq, double &thrs, double &resInterp, MultidimArray<double> &freq);

		void fscDir_fast(MultidimArray<double> &fsc, double rot, double tilt,
				         MetaData &mdRes, MultidimArray<double> &threeD_FSC, 
						 MultidimArray<double> &normalizationMap,
						 double &fscFreq, double &thrs, double &resol, size_t dirnumber);

};

#endif
