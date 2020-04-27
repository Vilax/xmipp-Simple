/***************************************************************************
 *
 * Authors:     Jose Luis Vilas (jlvilas@cnb.csic.es)
 *              Hemant. D. Tagare
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

#include <core/xmipp_program.h>
#include <core/xmipp_fftw.h>
#include <core/metadata_extension.h>
#include <data/monogenic.h>
#include <ctime>

class ProgResolutionDirectionalFsc : public XmippProgram
{
public:

    FileName    fnhalf1, fnhalf2, fnmask, fn_root, fn_out;
    double      sampling, ang_con;
    bool        test;

    Matrix2D<double> angles;
    FileName    fn_sel;

    void defineParams()
    {
        addUsageLine("Calculate the directional FSC resolution of two half maps along a set of directions");
        addUsageLine("+ ");
        addUsageLine("+* Fourier Shell Correlation (FSC)", true);
        addUsageLine("+This program may then be used to calculate FSC between these two volumes. The resulting plots");
        addUsageLine("+are commonly used to assess the high-resolution limit of the initial reconstruction.");
        addUsageLine(" ");
        addUsageLine("The program writes out filename.frc files, for each input volume or image, or selfilename.frc, the");
        addUsageLine("set_of_images mode. These ACSII files contain the DPR, FRC and SSNR as a function of resolution (in 1/Angstrom).");
        addUsageLine(" The .frc files also contain a column for the FRC expected for pure noise.");
        addSeeAlsoLine("resolution_ssnr");

        addParamsLine("   --half1 <input_file>      : Half map1");
        addParamsLine("   --half2 <input_file>      : Half map2");
        addParamsLine("   [--mask <input_file=\"\">]     : (Optional) Smooth mask to remove noise");
        addParamsLine("   [-o <output_file=\"\">]   : Output file name.");
        addParamsLine("   [--vol <Ts=20>]   : semi angle of the cone in degrees");
        addParamsLine("   [--sampling_rate <Ts=1>]  : Pixel size (Angstrom)");
        addParamsLine("   [--test]                  : It executes an unitary test");
        addExampleLine("Resolution of subset2.vol volume with respect to subset1.vol reference volume using 5.6 pixel size (in Angstrom):", false);
        addExampleLine("xmipp_resolution_fsc --ref subset1.vol  -i subset2.vol --sampling_rate 5.6 ");
        addExampleLine("Resolution of a set of images using 5.6 pixel size (in Angstrom):", false);
        addExampleLine("xmipp_resolution_fsc --set_of_images selfile.sel --sampling_rate 5.6");
    }

    void readParams()
    {
        fnhalf1 = getParam("--half1");
        fnhalf2 = getParam("--half2");
        fnmask = getParam("--mask");

        sampling = getDoubleParam("--sampling_rate");
//        ang_con  = getDoubleParam("--vol");

        fn_out = getParam("-o");
        test = checkParam("--test");
    }

    void defineFrequencies(MultidimArray<double> &vol,
    		MultidimArray<std::complex<double>> &fftvol,
    		MultidimArray<double> &freq, Matrix1D<double> &freq_fourier_vec)
    {
    	// Frequency volume
    	double uz, uy, ux, uz2, u2, uz2y2;
    	long n=0;

    	freq.resizeNoCopy(fftvol);

    	for(size_t k=0; k<ZSIZE(fftvol); ++k)
    	{
    		FFT_IDX2DIGFREQ(k,ZSIZE(vol),uz);
    		uz2=uz*uz;

    		for(size_t i=0; i<YSIZE(fftvol); ++i)
    		{
    			FFT_IDX2DIGFREQ(i,YSIZE(vol),uy);
    			uz2y2=uz2+uy*uy;

    			for(size_t j=0; j<XSIZE(fftvol); ++j)
    			{
    				FFT_IDX2DIGFREQ(j,XSIZE(vol), ux);
    				u2=uz2y2+ux*ux;
    				DIRECT_MULTIDIM_ELEM(freq,n) = sqrt(u2);
    				++n;
    			}
    		}
    	}

    	double u;
    	size_t dimfft, dimreal;

    	if (XSIZE(fftvol) <= YSIZE(fftvol))
    		if (XSIZE(fftvol) <= ZSIZE(fftvol)){
    			dimfft = XSIZE(fftvol);
    			dimreal = XSIZE(vol);
    		}

    	if (YSIZE(fftvol) <= XSIZE(fftvol))
    		if (YSIZE(fftvol) <= ZSIZE(fftvol)){
    			dimfft = YSIZE(fftvol);
    			dimreal = YSIZE(vol);
    		}

    	if (ZSIZE(fftvol) <= XSIZE(fftvol))
    		if (ZSIZE(fftvol) <= YSIZE(fftvol)){
    			dimfft = ZSIZE(fftvol);
    			dimreal = ZSIZE(vol);
    		}

    	freq_fourier_vec.initZeros(dimfft);

    	VEC_ELEM(freq_fourier_vec,0) = 1e-38;
    	for(size_t k=1; k<dimfft; ++k){
    		FFT_IDX2DIGFREQ(k,dimreal, u);
    		VEC_ELEM(freq_fourier_vec,k) = u;
    	}



    }
    MultidimArray<double> defineFrequencies(const MultidimArray< std::complex<double> > &myfftV,
    		const MultidimArray<double> &inputVol,
    		Matrix1D<double> &freq_fourier_x,
    		Matrix1D<double> &freq_fourier_y,
    		Matrix1D<double> &freq_fourier_z)
    {
    	double u;

    	freq_fourier_z.initZeros(ZSIZE(myfftV));
    	freq_fourier_x.initZeros(XSIZE(myfftV));
    	freq_fourier_y.initZeros(YSIZE(myfftV));

    	VEC_ELEM(freq_fourier_z,0) = 1e-38;
    	for(size_t k=1; k<ZSIZE(myfftV); ++k){
    		FFT_IDX2DIGFREQ(k,ZSIZE(inputVol), u);
    		VEC_ELEM(freq_fourier_z,k) = u;
    	}

    	VEC_ELEM(freq_fourier_y,0) = 1e-38;
    	for(size_t k=1; k<YSIZE(myfftV); ++k){
    		FFT_IDX2DIGFREQ(k,YSIZE(inputVol), u);
    		VEC_ELEM(freq_fourier_y,k) = u;
    	}

    	VEC_ELEM(freq_fourier_x,0) = 1e-38;
    	for(size_t k=1; k<XSIZE(myfftV); ++k){
    		FFT_IDX2DIGFREQ(k,XSIZE(inputVol), u);
    		VEC_ELEM(freq_fourier_x,k) = u;
    	}


    	MultidimArray<double> iu;

    	iu.initZeros(myfftV);

    	double uz, uy, ux, uz2, u2, uz2y2;
    	long n=0;
    	//  TODO: reasign uz = uz*uz to save memory
    	//  TODO: Take ZSIZE(myfftV) out of the loop
    	//	TODO: Use freq_fourier_x instead of calling FFT_IDX2DIGFREQ

    	for(size_t k=0; k<ZSIZE(myfftV); ++k)
    	{
    		FFT_IDX2DIGFREQ(k,ZSIZE(inputVol),uz);
    		uz2 = uz*uz;
    		for(size_t i=0; i<YSIZE(myfftV); ++i)
    		{
    			FFT_IDX2DIGFREQ(i,YSIZE(inputVol),uy);
    			uz2y2 = uz2 + uy*uy;

    			for(size_t j=0; j<XSIZE(myfftV); ++j)
    			{
    				FFT_IDX2DIGFREQ(j,XSIZE(inputVol), ux);
    				u2 = uz2y2 + ux*ux;
//   					DIRECT_MULTIDIM_ELEM(iu,n) = sqrt(u2);
   					if ((k != 0) || (i != 0) || (j != 0))
   					{
   						DIRECT_MULTIDIM_ELEM(iu,n) = 1/sqrt(u2);
   					}
   					else
   					{
   						DIRECT_MULTIDIM_ELEM(iu,n) = 1e38;
   					}
    				++n;
    			}
    		}
    	}


    	return iu;
    }


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
				 double rot, double tilt, double ang_con)
    {
        MultidimArray< int > radial_count(m1sizeX/2+1);
        MultidimArray<double> num, den1, den2;
//        MultidimArray<float> num, den1, den2;
//        Matrix1D<double> f(3);

        num.initZeros(radial_count);
        den1.initZeros(radial_count);
        den2.initZeros(radial_count);

        freq.initZeros(radial_count);
        frc.initZeros(radial_count);

        int ZdimFT1=(int)ZSIZE(FT1);
        int YdimFT1=(int)YSIZE(FT1);
        int XdimFT1=(int)XSIZE(FT1);

        double maxFreq_2 =0.;
        maxFreq_2 = maxFreq;

    	double x_dir, y_dir, z_dir, uz, uy, ux, cosAngle;
    	x_dir = sin(tilt*PI/180)*cos(rot*PI/180);
    	y_dir = sin(tilt*PI/180)*sin(rot*PI/180);
    	z_dir = cos(tilt*PI/180);
    	cosAngle = cos(ang_con);

        long n = 0;
        for (int k=0; k<ZdimFT1; k++)
        {
            double uz = VEC_ELEM(freq_fourier_z,k);
            uz *= z_dir;
            for (int i=0; i<YdimFT1; i++)
            {
            	double uy = VEC_ELEM(freq_fourier_y,i);
                uy *= y_dir;
                for (int j=0; j<XdimFT1; j++)
                {
                	double ux = VEC_ELEM(freq_fourier_x,j);
                    ux *= x_dir;
                    double iun = DIRECT_MULTIDIM_ELEM(freqMap,n);
                    double f = 1/iun;
                    iun *= (ux + uy + uz);

                    double cosine =fabs(iun);
                    ++n;

					if (cosine>=cosAngle)
						{
							if (f>maxFreq_2)
								continue;

							int idx = (int) round(f * m1sizeX);
							std::complex<double> &z1 = dAkij(FT1, k, i, j);
							std::complex<double> &z2 = dAkij(FT2, k, i, j);
							double absz1 = abs(z1);
							double absz2 = abs(z2);
							dAi(num,idx) += real(conj(z1) * z2);
							dAi(den1,idx) += absz1*absz1;
							dAi(den2,idx) += absz2*absz2;
						}
                }
            }
        }

        FOR_ALL_ELEMENTS_IN_ARRAY1D(freq)
        {
            dAi(freq,i) = (float) i / (m1sizeX * sampling_rate);
            dAi(frc,i) = dAi(num,i)/sqrt(dAi(den1,i)*dAi(den2,i));
        }
    }


    void generateDirections(Matrix2D<double> &angles, bool alot)
    {
    	if (alot == true)
    	{
    		angles.initZeros(2,384);
    		MAT_ELEM(angles,0,0) =85.2198;   MAT_ELEM(angles,1,0) =45;
    		MAT_ELEM(angles,0,1) =80.4059;   MAT_ELEM(angles,1,1) =50.625;
    		MAT_ELEM(angles,0,2) =80.4059;   MAT_ELEM(angles,1,2) =39.375;
    		MAT_ELEM(angles,0,3) =75.5225;   MAT_ELEM(angles,1,3) =45;
    		MAT_ELEM(angles,0,4) =75.5225;   MAT_ELEM(angles,1,4) =56.25;
    		MAT_ELEM(angles,0,5) =70.5288;   MAT_ELEM(angles,1,5) =61.875;
    		MAT_ELEM(angles,0,6) =70.5288;   MAT_ELEM(angles,1,6) =50.625;
    		MAT_ELEM(angles,0,7) =65.3757;   MAT_ELEM(angles,1,7) =56.25;
    		MAT_ELEM(angles,0,8) =75.5225;   MAT_ELEM(angles,1,8) =33.75;
    		MAT_ELEM(angles,0,9) =70.5288;   MAT_ELEM(angles,1,9) =39.375;
    		MAT_ELEM(angles,0,10) =70.5288;   MAT_ELEM(angles,1,10) =28.125;
    		MAT_ELEM(angles,0,11) =65.3757;   MAT_ELEM(angles,1,11) =33.75;
    		MAT_ELEM(angles,0,12) =65.3757;   MAT_ELEM(angles,1,12) =45;
    		MAT_ELEM(angles,0,13) =60;   MAT_ELEM(angles,1,13) =50.625;
    		MAT_ELEM(angles,0,14) =60;   MAT_ELEM(angles,1,14) =39.375;
    		MAT_ELEM(angles,0,15) =54.3147;   MAT_ELEM(angles,1,15) =45;
    		MAT_ELEM(angles,0,16) =65.3757;   MAT_ELEM(angles,1,16) =67.5;
    		MAT_ELEM(angles,0,17) =60;   MAT_ELEM(angles,1,17) =73.125;
    		MAT_ELEM(angles,0,18) =60;   MAT_ELEM(angles,1,18) =61.875;
    		MAT_ELEM(angles,0,19) =54.3147;   MAT_ELEM(angles,1,19) =67.5;
    		MAT_ELEM(angles,0,20) =54.3147;   MAT_ELEM(angles,1,20) =78.75;
    		MAT_ELEM(angles,0,21) =48.1897;   MAT_ELEM(angles,1,21) =84.375;
    		MAT_ELEM(angles,0,22) =48.1897;   MAT_ELEM(angles,1,22) =73.125;
    		MAT_ELEM(angles,0,23) =41.8588;   MAT_ELEM(angles,1,23) =83.5714;
    		MAT_ELEM(angles,0,24) =54.3147;   MAT_ELEM(angles,1,24) =56.25;
    		MAT_ELEM(angles,0,25) =48.1897;   MAT_ELEM(angles,1,25) =61.875;
    		MAT_ELEM(angles,0,26) =48.1897;   MAT_ELEM(angles,1,26) =50.625;
    		MAT_ELEM(angles,0,27) =41.8588;   MAT_ELEM(angles,1,27) =57.8571;
    		MAT_ELEM(angles,0,28) =41.8588;   MAT_ELEM(angles,1,28) =70.7143;
    		MAT_ELEM(angles,0,29) =35.6591;   MAT_ELEM(angles,1,29) =82.5;
    		MAT_ELEM(angles,0,30) =35.6591;   MAT_ELEM(angles,1,30) =67.5;
    		MAT_ELEM(angles,0,31) =29.5656;   MAT_ELEM(angles,1,31) =81;
    		MAT_ELEM(angles,0,32) =65.3757;   MAT_ELEM(angles,1,32) =22.5;
    		MAT_ELEM(angles,0,33) =60;   MAT_ELEM(angles,1,33) =28.125;
    		MAT_ELEM(angles,0,34) =60;   MAT_ELEM(angles,1,34) =16.875;
    		MAT_ELEM(angles,0,35) =54.3147;   MAT_ELEM(angles,1,35) =22.5;
    		MAT_ELEM(angles,0,36) =54.3147;   MAT_ELEM(angles,1,36) =33.75;
    		MAT_ELEM(angles,0,37) =48.1897;   MAT_ELEM(angles,1,37) =39.375;
    		MAT_ELEM(angles,0,38) =48.1897;   MAT_ELEM(angles,1,38) =28.125;
    		MAT_ELEM(angles,0,39) =41.8588;   MAT_ELEM(angles,1,39) =32.1429;
    		MAT_ELEM(angles,0,40) =54.3147;   MAT_ELEM(angles,1,40) =11.25;
    		MAT_ELEM(angles,0,41) =48.1897;   MAT_ELEM(angles,1,41) =16.875;
    		MAT_ELEM(angles,0,42) =48.1897;   MAT_ELEM(angles,1,42) =5.625;
    		MAT_ELEM(angles,0,43) =41.8588;   MAT_ELEM(angles,1,43) =6.4286;
    		MAT_ELEM(angles,0,44) =41.8588;   MAT_ELEM(angles,1,44) =19.2857;
    		MAT_ELEM(angles,0,45) =35.6591;   MAT_ELEM(angles,1,45) =22.5;
    		MAT_ELEM(angles,0,46) =35.6591;   MAT_ELEM(angles,1,46) =7.5;
    		MAT_ELEM(angles,0,47) =29.5656;   MAT_ELEM(angles,1,47) =9;
    		MAT_ELEM(angles,0,48) =41.8588;   MAT_ELEM(angles,1,48) =45;
    		MAT_ELEM(angles,0,49) =35.6591;   MAT_ELEM(angles,1,49) =52.5;
    		MAT_ELEM(angles,0,50) =35.6591;   MAT_ELEM(angles,1,50) =37.5;
    		MAT_ELEM(angles,0,51) =29.5656;   MAT_ELEM(angles,1,51) =45;
    		MAT_ELEM(angles,0,52) =29.5656;   MAT_ELEM(angles,1,52) =63;
    		MAT_ELEM(angles,0,53) =23.5565;   MAT_ELEM(angles,1,53) =78.75;
    		MAT_ELEM(angles,0,54) =23.5565;   MAT_ELEM(angles,1,54) =56.25;
    		MAT_ELEM(angles,0,55) =17.6124;   MAT_ELEM(angles,1,55) =75;
    		MAT_ELEM(angles,0,56) =29.5656;   MAT_ELEM(angles,1,56) =27;
    		MAT_ELEM(angles,0,57) =23.5565;   MAT_ELEM(angles,1,57) =33.75;
    		MAT_ELEM(angles,0,58) =23.5565;   MAT_ELEM(angles,1,58) =11.25;
    		MAT_ELEM(angles,0,59) =17.6124;   MAT_ELEM(angles,1,59) =15;
    		MAT_ELEM(angles,0,60) =17.6124;   MAT_ELEM(angles,1,60) =45;
    		MAT_ELEM(angles,0,61) =11.7159;   MAT_ELEM(angles,1,61) =67.5;
    		MAT_ELEM(angles,0,62) =11.7159;   MAT_ELEM(angles,1,62) =22.5;
    		MAT_ELEM(angles,0,63) =5.8503;   MAT_ELEM(angles,1,63) =45;
    		MAT_ELEM(angles,0,64) =85.2198;   MAT_ELEM(angles,1,64) =-45;
    		MAT_ELEM(angles,0,65) =80.4059;   MAT_ELEM(angles,1,65) =-39.375;
    		MAT_ELEM(angles,0,66) =80.4059;   MAT_ELEM(angles,1,66) =-50.625;
    		MAT_ELEM(angles,0,67) =75.5225;   MAT_ELEM(angles,1,67) =-45;
    		MAT_ELEM(angles,0,68) =75.5225;   MAT_ELEM(angles,1,68) =-33.75;
    		MAT_ELEM(angles,0,69) =70.5288;   MAT_ELEM(angles,1,69) =-28.125;
    		MAT_ELEM(angles,0,70) =70.5288;   MAT_ELEM(angles,1,70) =-39.375;
    		MAT_ELEM(angles,0,71) =65.3757;   MAT_ELEM(angles,1,71) =-33.75;
    		MAT_ELEM(angles,0,72) =75.5225;   MAT_ELEM(angles,1,72) =-56.25;
    		MAT_ELEM(angles,0,73) =70.5288;   MAT_ELEM(angles,1,73) =-50.625;
    		MAT_ELEM(angles,0,74) =70.5288;   MAT_ELEM(angles,1,74) =-61.875;
    		MAT_ELEM(angles,0,75) =65.3757;   MAT_ELEM(angles,1,75) =-56.25;
    		MAT_ELEM(angles,0,76) =65.3757;   MAT_ELEM(angles,1,76) =-45;
    		MAT_ELEM(angles,0,77) =60;   MAT_ELEM(angles,1,77) =-39.375;
    		MAT_ELEM(angles,0,78) =60;   MAT_ELEM(angles,1,78) =-50.625;
    		MAT_ELEM(angles,0,79) =54.3147;   MAT_ELEM(angles,1,79) =-45;
    		MAT_ELEM(angles,0,80) =65.3757;   MAT_ELEM(angles,1,80) =-22.5;
    		MAT_ELEM(angles,0,81) =60;   MAT_ELEM(angles,1,81) =-16.875;
    		MAT_ELEM(angles,0,82) =60;   MAT_ELEM(angles,1,82) =-28.125;
    		MAT_ELEM(angles,0,83) =54.3147;   MAT_ELEM(angles,1,83) =-22.5;
    		MAT_ELEM(angles,0,84) =54.3147;   MAT_ELEM(angles,1,84) =-11.25;
    		MAT_ELEM(angles,0,85) =48.1897;   MAT_ELEM(angles,1,85) =-5.625;
    		MAT_ELEM(angles,0,86) =48.1897;   MAT_ELEM(angles,1,86) =-16.875;
    		MAT_ELEM(angles,0,87) =41.8588;   MAT_ELEM(angles,1,87) =-6.4286;
    		MAT_ELEM(angles,0,88) =54.3147;   MAT_ELEM(angles,1,88) =-33.75;
    		MAT_ELEM(angles,0,89) =48.1897;   MAT_ELEM(angles,1,89) =-28.125;
    		MAT_ELEM(angles,0,90) =48.1897;   MAT_ELEM(angles,1,90) =-39.375;
    		MAT_ELEM(angles,0,91) =41.8588;   MAT_ELEM(angles,1,91) =-32.1429;
    		MAT_ELEM(angles,0,92) =41.8588;   MAT_ELEM(angles,1,92) =-19.2857;
    		MAT_ELEM(angles,0,93) =35.6591;   MAT_ELEM(angles,1,93) =-7.5;
    		MAT_ELEM(angles,0,94) =35.6591;   MAT_ELEM(angles,1,94) =-22.5;
    		MAT_ELEM(angles,0,95) =29.5656;   MAT_ELEM(angles,1,95) =-9;
    		MAT_ELEM(angles,0,96) =65.3757;   MAT_ELEM(angles,1,96) =-67.5;
    		MAT_ELEM(angles,0,97) =60;   MAT_ELEM(angles,1,97) =-61.875;
    		MAT_ELEM(angles,0,98) =60;   MAT_ELEM(angles,1,98) =-73.125;
    		MAT_ELEM(angles,0,99) =54.3147;   MAT_ELEM(angles,1,99) =-67.5;
    		MAT_ELEM(angles,0,100) =54.3147;   MAT_ELEM(angles,1,100) =-56.25;
    		MAT_ELEM(angles,0,101) =48.1897;   MAT_ELEM(angles,1,101) =-50.625;
    		MAT_ELEM(angles,0,102) =48.1897;   MAT_ELEM(angles,1,102) =-61.875;
    		MAT_ELEM(angles,0,103) =41.8588;   MAT_ELEM(angles,1,103) =-57.8571;
    		MAT_ELEM(angles,0,104) =54.3147;   MAT_ELEM(angles,1,104) =-78.75;
    		MAT_ELEM(angles,0,105) =48.1897;   MAT_ELEM(angles,1,105) =-73.125;
    		MAT_ELEM(angles,0,106) =48.1897;   MAT_ELEM(angles,1,106) =-84.375;
    		MAT_ELEM(angles,0,107) =41.8588;   MAT_ELEM(angles,1,107) =-83.5714;
    		MAT_ELEM(angles,0,108) =41.8588;   MAT_ELEM(angles,1,108) =-70.7143;
    		MAT_ELEM(angles,0,109) =35.6591;   MAT_ELEM(angles,1,109) =-67.5;
    		MAT_ELEM(angles,0,110) =35.6591;   MAT_ELEM(angles,1,110) =-82.5;
    		MAT_ELEM(angles,0,111) =29.5656;   MAT_ELEM(angles,1,111) =-81;
    		MAT_ELEM(angles,0,112) =41.8588;   MAT_ELEM(angles,1,112) =-45;
    		MAT_ELEM(angles,0,113) =35.6591;   MAT_ELEM(angles,1,113) =-37.5;
    		MAT_ELEM(angles,0,114) =35.6591;   MAT_ELEM(angles,1,114) =-52.5;
    		MAT_ELEM(angles,0,115) =29.5656;   MAT_ELEM(angles,1,115) =-45;
    		MAT_ELEM(angles,0,116) =29.5656;   MAT_ELEM(angles,1,116) =-27;
    		MAT_ELEM(angles,0,117) =23.5565;   MAT_ELEM(angles,1,117) =-11.25;
    		MAT_ELEM(angles,0,118) =23.5565;   MAT_ELEM(angles,1,118) =-33.75;
    		MAT_ELEM(angles,0,119) =17.6124;   MAT_ELEM(angles,1,119) =-15;
    		MAT_ELEM(angles,0,120) =29.5656;   MAT_ELEM(angles,1,120) =-63;
    		MAT_ELEM(angles,0,121) =23.5565;   MAT_ELEM(angles,1,121) =-56.25;
    		MAT_ELEM(angles,0,122) =23.5565;   MAT_ELEM(angles,1,122) =-78.75;
    		MAT_ELEM(angles,0,123) =17.6124;   MAT_ELEM(angles,1,123) =-75;
    		MAT_ELEM(angles,0,124) =17.6124;   MAT_ELEM(angles,1,124) =-45;
    		MAT_ELEM(angles,0,125) =11.7159;   MAT_ELEM(angles,1,125) =-22.5;
    		MAT_ELEM(angles,0,126) =11.7159;   MAT_ELEM(angles,1,126) =-67.5;
    		MAT_ELEM(angles,0,127) =5.8503;   MAT_ELEM(angles,1,127) =-45;
    		MAT_ELEM(angles,0,128) =125.6853;   MAT_ELEM(angles,1,128) =-1.4033e-14;
    		MAT_ELEM(angles,0,129) =120;   MAT_ELEM(angles,1,129) =5.625;
    		MAT_ELEM(angles,0,130) =120;   MAT_ELEM(angles,1,130) =-5.625;
    		MAT_ELEM(angles,0,131) =114.6243;   MAT_ELEM(angles,1,131) =-1.4033e-14;
    		MAT_ELEM(angles,0,132) =114.6243;   MAT_ELEM(angles,1,132) =11.25;
    		MAT_ELEM(angles,0,133) =109.4712;   MAT_ELEM(angles,1,133) =16.875;
    		MAT_ELEM(angles,0,134) =109.4712;   MAT_ELEM(angles,1,134) =5.625;
    		MAT_ELEM(angles,0,135) =104.4775;   MAT_ELEM(angles,1,135) =11.25;
    		MAT_ELEM(angles,0,136) =114.6243;   MAT_ELEM(angles,1,136) =-11.25;
    		MAT_ELEM(angles,0,137) =109.4712;   MAT_ELEM(angles,1,137) =-5.625;
    		MAT_ELEM(angles,0,138) =109.4712;   MAT_ELEM(angles,1,138) =-16.875;
    		MAT_ELEM(angles,0,139) =104.4775;   MAT_ELEM(angles,1,139) =-11.25;
    		MAT_ELEM(angles,0,140) =104.4775;   MAT_ELEM(angles,1,140) =-1.4033e-14;
    		MAT_ELEM(angles,0,141) =99.5941;   MAT_ELEM(angles,1,141) =5.625;
    		MAT_ELEM(angles,0,142) =99.5941;   MAT_ELEM(angles,1,142) =-5.625;
    		MAT_ELEM(angles,0,143) =94.7802;   MAT_ELEM(angles,1,143) =-1.4033e-14;
    		MAT_ELEM(angles,0,144) =104.4775;   MAT_ELEM(angles,1,144) =22.5;
    		MAT_ELEM(angles,0,145) =99.5941;   MAT_ELEM(angles,1,145) =28.125;
    		MAT_ELEM(angles,0,146) =99.5941;   MAT_ELEM(angles,1,146) =16.875;
    		MAT_ELEM(angles,0,147) =94.7802;   MAT_ELEM(angles,1,147) =22.5;
    		MAT_ELEM(angles,0,148) =94.7802;   MAT_ELEM(angles,1,148) =33.75;
    		MAT_ELEM(angles,0,149) =90;   MAT_ELEM(angles,1,149) =39.375;
    		MAT_ELEM(angles,0,150) =90;   MAT_ELEM(angles,1,150) =28.125;
    		MAT_ELEM(angles,0,151) =85.2198;   MAT_ELEM(angles,1,151) =33.75;
    		MAT_ELEM(angles,0,152) =94.7802;   MAT_ELEM(angles,1,152) =11.25;
    		MAT_ELEM(angles,0,153) =90;   MAT_ELEM(angles,1,153) =16.875;
    		MAT_ELEM(angles,0,154) =90;   MAT_ELEM(angles,1,154) =5.625;
    		MAT_ELEM(angles,0,155) =85.2198;   MAT_ELEM(angles,1,155) =11.25;
    		MAT_ELEM(angles,0,156) =85.2198;   MAT_ELEM(angles,1,156) =22.5;
    		MAT_ELEM(angles,0,157) =80.4059;   MAT_ELEM(angles,1,157) =28.125;
    		MAT_ELEM(angles,0,158) =80.4059;   MAT_ELEM(angles,1,158) =16.875;
    		MAT_ELEM(angles,0,159) =75.5225;   MAT_ELEM(angles,1,159) =22.5;
    		MAT_ELEM(angles,0,160) =104.4775;   MAT_ELEM(angles,1,160) =-22.5;
    		MAT_ELEM(angles,0,161) =99.5941;   MAT_ELEM(angles,1,161) =-16.875;
    		MAT_ELEM(angles,0,162) =99.5941;   MAT_ELEM(angles,1,162) =-28.125;
    		MAT_ELEM(angles,0,163) =94.7802;   MAT_ELEM(angles,1,163) =-22.5;
    		MAT_ELEM(angles,0,164) =94.7802;   MAT_ELEM(angles,1,164) =-11.25;
    		MAT_ELEM(angles,0,165) =90;   MAT_ELEM(angles,1,165) =-5.625;
    		MAT_ELEM(angles,0,166) =90;   MAT_ELEM(angles,1,166) =-16.875;
    		MAT_ELEM(angles,0,167) =85.2198;   MAT_ELEM(angles,1,167) =-11.25;
    		MAT_ELEM(angles,0,168) =94.7802;   MAT_ELEM(angles,1,168) =-33.75;
    		MAT_ELEM(angles,0,169) =90;   MAT_ELEM(angles,1,169) =-28.125;
    		MAT_ELEM(angles,0,170) =90;   MAT_ELEM(angles,1,170) =-39.375;
    		MAT_ELEM(angles,0,171) =85.2198;   MAT_ELEM(angles,1,171) =-33.75;
    		MAT_ELEM(angles,0,172) =85.2198;   MAT_ELEM(angles,1,172) =-22.5;
    		MAT_ELEM(angles,0,173) =80.4059;   MAT_ELEM(angles,1,173) =-16.875;
    		MAT_ELEM(angles,0,174) =80.4059;   MAT_ELEM(angles,1,174) =-28.125;
    		MAT_ELEM(angles,0,175) =75.5225;   MAT_ELEM(angles,1,175) =-22.5;
    		MAT_ELEM(angles,0,176) =85.2198;   MAT_ELEM(angles,1,176) =-1.4033e-14;
    		MAT_ELEM(angles,0,177) =80.4059;   MAT_ELEM(angles,1,177) =5.625;
    		MAT_ELEM(angles,0,178) =80.4059;   MAT_ELEM(angles,1,178) =-5.625;
    		MAT_ELEM(angles,0,179) =75.5225;   MAT_ELEM(angles,1,179) =-1.4033e-14;
    		MAT_ELEM(angles,0,180) =75.5225;   MAT_ELEM(angles,1,180) =11.25;
    		MAT_ELEM(angles,0,181) =70.5288;   MAT_ELEM(angles,1,181) =16.875;
    		MAT_ELEM(angles,0,182) =70.5288;   MAT_ELEM(angles,1,182) =5.625;
    		MAT_ELEM(angles,0,183) =65.3757;   MAT_ELEM(angles,1,183) =11.25;
    		MAT_ELEM(angles,0,184) =75.5225;   MAT_ELEM(angles,1,184) =-11.25;
    		MAT_ELEM(angles,0,185) =70.5288;   MAT_ELEM(angles,1,185) =-5.625;
    		MAT_ELEM(angles,0,186) =70.5288;   MAT_ELEM(angles,1,186) =-16.875;
    		MAT_ELEM(angles,0,187) =65.3757;   MAT_ELEM(angles,1,187) =-11.25;
    		MAT_ELEM(angles,0,188) =65.3757;   MAT_ELEM(angles,1,188) =-1.4033e-14;
    		MAT_ELEM(angles,0,189) =60;   MAT_ELEM(angles,1,189) =5.625;
    		MAT_ELEM(angles,0,190) =60;   MAT_ELEM(angles,1,190) =-5.625;
    		MAT_ELEM(angles,0,191) =54.3147;   MAT_ELEM(angles,1,191) =-1.4033e-14;
    		MAT_ELEM(angles,0,192) =125.6853;   MAT_ELEM(angles,1,192) =90;
    		MAT_ELEM(angles,0,193) =120;   MAT_ELEM(angles,1,193) =84.375;
    		MAT_ELEM(angles,0,194) =114.6243;   MAT_ELEM(angles,1,194) =90;
    		MAT_ELEM(angles,0,195) =114.6243;   MAT_ELEM(angles,1,195) =78.75;
    		MAT_ELEM(angles,0,196) =109.4712;   MAT_ELEM(angles,1,196) =84.375;
    		MAT_ELEM(angles,0,197) =109.4712;   MAT_ELEM(angles,1,197) =73.125;
    		MAT_ELEM(angles,0,198) =104.4775;   MAT_ELEM(angles,1,198) =78.75;
    		MAT_ELEM(angles,0,199) =104.4775;   MAT_ELEM(angles,1,199) =90;
    		MAT_ELEM(angles,0,200) =99.5941;   MAT_ELEM(angles,1,200) =84.375;
    		MAT_ELEM(angles,0,201) =94.7802;   MAT_ELEM(angles,1,201) =90;
    		MAT_ELEM(angles,0,202) =104.4775;   MAT_ELEM(angles,1,202) =67.5;
    		MAT_ELEM(angles,0,203) =99.5941;   MAT_ELEM(angles,1,203) =73.125;
    		MAT_ELEM(angles,0,204) =99.5941;   MAT_ELEM(angles,1,204) =61.875;
    		MAT_ELEM(angles,0,205) =94.7802;   MAT_ELEM(angles,1,205) =67.5;
    		MAT_ELEM(angles,0,206) =94.7802;   MAT_ELEM(angles,1,206) =78.75;
    		MAT_ELEM(angles,0,207) =90;   MAT_ELEM(angles,1,207) =84.375;
    		MAT_ELEM(angles,0,208) =90;   MAT_ELEM(angles,1,208) =73.125;
    		MAT_ELEM(angles,0,209) =85.2198;   MAT_ELEM(angles,1,209) =78.75;
    		MAT_ELEM(angles,0,210) =94.7802;   MAT_ELEM(angles,1,210) =56.25;
    		MAT_ELEM(angles,0,211) =90;   MAT_ELEM(angles,1,211) =61.875;
    		MAT_ELEM(angles,0,212) =90;   MAT_ELEM(angles,1,212) =50.625;
    		MAT_ELEM(angles,0,213) =85.2198;   MAT_ELEM(angles,1,213) =56.25;
    		MAT_ELEM(angles,0,214) =85.2198;   MAT_ELEM(angles,1,214) =67.5;
    		MAT_ELEM(angles,0,215) =80.4059;   MAT_ELEM(angles,1,215) =73.125;
    		MAT_ELEM(angles,0,216) =80.4059;   MAT_ELEM(angles,1,216) =61.875;
    		MAT_ELEM(angles,0,217) =75.5225;   MAT_ELEM(angles,1,217) =67.5;
    		MAT_ELEM(angles,0,218) =85.2198;   MAT_ELEM(angles,1,218) =90;
    		MAT_ELEM(angles,0,219) =80.4059;   MAT_ELEM(angles,1,219) =84.375;
    		MAT_ELEM(angles,0,220) =75.5225;   MAT_ELEM(angles,1,220) =90;
    		MAT_ELEM(angles,0,221) =75.5225;   MAT_ELEM(angles,1,221) =78.75;
    		MAT_ELEM(angles,0,222) =70.5288;   MAT_ELEM(angles,1,222) =84.375;
    		MAT_ELEM(angles,0,223) =70.5288;   MAT_ELEM(angles,1,223) =73.125;
    		MAT_ELEM(angles,0,224) =65.3757;   MAT_ELEM(angles,1,224) =78.75;
    		MAT_ELEM(angles,0,225) =65.3757;   MAT_ELEM(angles,1,225) =90;
    		MAT_ELEM(angles,0,226) =60;   MAT_ELEM(angles,1,226) =84.375;
    		MAT_ELEM(angles,0,227) =54.3147;   MAT_ELEM(angles,1,227) =90;
    		MAT_ELEM(angles,0,228) =120;   MAT_ELEM(angles,1,228) =-84.375;
    		MAT_ELEM(angles,0,229) =114.6243;   MAT_ELEM(angles,1,229) =-78.75;
    		MAT_ELEM(angles,0,230) =109.4712;   MAT_ELEM(angles,1,230) =-73.125;
    		MAT_ELEM(angles,0,231) =109.4712;   MAT_ELEM(angles,1,231) =-84.375;
    		MAT_ELEM(angles,0,232) =104.4775;   MAT_ELEM(angles,1,232) =-78.75;
    		MAT_ELEM(angles,0,233) =99.5941;   MAT_ELEM(angles,1,233) =-84.375;
    		MAT_ELEM(angles,0,234) =104.4775;   MAT_ELEM(angles,1,234) =-67.5;
    		MAT_ELEM(angles,0,235) =99.5941;   MAT_ELEM(angles,1,235) =-61.875;
    		MAT_ELEM(angles,0,236) =99.5941;   MAT_ELEM(angles,1,236) =-73.125;
    		MAT_ELEM(angles,0,237) =94.7802;   MAT_ELEM(angles,1,237) =-67.5;
    		MAT_ELEM(angles,0,238) =94.7802;   MAT_ELEM(angles,1,238) =-56.25;
    		MAT_ELEM(angles,0,239) =90;   MAT_ELEM(angles,1,239) =-50.625;
    		MAT_ELEM(angles,0,240) =90;   MAT_ELEM(angles,1,240) =-61.875;
    		MAT_ELEM(angles,0,241) =85.2198;   MAT_ELEM(angles,1,241) =-56.25;
    		MAT_ELEM(angles,0,242) =94.7802;   MAT_ELEM(angles,1,242) =-78.75;
    		MAT_ELEM(angles,0,243) =90;   MAT_ELEM(angles,1,243) =-73.125;
    		MAT_ELEM(angles,0,244) =90;   MAT_ELEM(angles,1,244) =-84.375;
    		MAT_ELEM(angles,0,245) =85.2198;   MAT_ELEM(angles,1,245) =-78.75;
    		MAT_ELEM(angles,0,246) =85.2198;   MAT_ELEM(angles,1,246) =-67.5;
    		MAT_ELEM(angles,0,247) =80.4059;   MAT_ELEM(angles,1,247) =-61.875;
    		MAT_ELEM(angles,0,248) =80.4059;   MAT_ELEM(angles,1,248) =-73.125;
    		MAT_ELEM(angles,0,249) =75.5225;   MAT_ELEM(angles,1,249) =-67.5;
    		MAT_ELEM(angles,0,250) =80.4059;   MAT_ELEM(angles,1,250) =-84.375;
    		MAT_ELEM(angles,0,251) =75.5225;   MAT_ELEM(angles,1,251) =-78.75;
    		MAT_ELEM(angles,0,252) =70.5288;   MAT_ELEM(angles,1,252) =-73.125;
    		MAT_ELEM(angles,0,253) =70.5288;   MAT_ELEM(angles,1,253) =-84.375;
    		MAT_ELEM(angles,0,254) =65.3757;   MAT_ELEM(angles,1,254) =-78.75;
    		MAT_ELEM(angles,0,255) =60;   MAT_ELEM(angles,1,255) =-84.375;
    		MAT_ELEM(angles,0,256) =174.1497;   MAT_ELEM(angles,1,256) =45;
    		MAT_ELEM(angles,0,257) =168.2841;   MAT_ELEM(angles,1,257) =67.5;
    		MAT_ELEM(angles,0,258) =168.2841;   MAT_ELEM(angles,1,258) =22.5;
    		MAT_ELEM(angles,0,259) =162.3876;   MAT_ELEM(angles,1,259) =45;
    		MAT_ELEM(angles,0,260) =162.3876;   MAT_ELEM(angles,1,260) =75;
    		MAT_ELEM(angles,0,261) =156.4435;   MAT_ELEM(angles,1,261) =78.75;
    		MAT_ELEM(angles,0,262) =156.4435;   MAT_ELEM(angles,1,262) =56.25;
    		MAT_ELEM(angles,0,263) =150.4344;   MAT_ELEM(angles,1,263) =63;
    		MAT_ELEM(angles,0,264) =162.3876;   MAT_ELEM(angles,1,264) =15;
    		MAT_ELEM(angles,0,265) =156.4435;   MAT_ELEM(angles,1,265) =33.75;
    		MAT_ELEM(angles,0,266) =156.4435;   MAT_ELEM(angles,1,266) =11.25;
    		MAT_ELEM(angles,0,267) =150.4344;   MAT_ELEM(angles,1,267) =27;
    		MAT_ELEM(angles,0,268) =150.4344;   MAT_ELEM(angles,1,268) =45;
    		MAT_ELEM(angles,0,269) =144.3409;   MAT_ELEM(angles,1,269) =52.5;
    		MAT_ELEM(angles,0,270) =144.3409;   MAT_ELEM(angles,1,270) =37.5;
    		MAT_ELEM(angles,0,271) =138.1412;   MAT_ELEM(angles,1,271) =45;
    		MAT_ELEM(angles,0,272) =150.4344;   MAT_ELEM(angles,1,272) =81;
    		MAT_ELEM(angles,0,273) =144.3409;   MAT_ELEM(angles,1,273) =82.5;
    		MAT_ELEM(angles,0,274) =144.3409;   MAT_ELEM(angles,1,274) =67.5;
    		MAT_ELEM(angles,0,275) =138.1412;   MAT_ELEM(angles,1,275) =70.7143;
    		MAT_ELEM(angles,0,276) =138.1412;   MAT_ELEM(angles,1,276) =83.5714;
    		MAT_ELEM(angles,0,277) =131.8103;   MAT_ELEM(angles,1,277) =84.375;
    		MAT_ELEM(angles,0,278) =131.8103;   MAT_ELEM(angles,1,278) =73.125;
    		MAT_ELEM(angles,0,279) =125.6853;   MAT_ELEM(angles,1,279) =78.75;
    		MAT_ELEM(angles,0,280) =138.1412;   MAT_ELEM(angles,1,280) =57.8571;
    		MAT_ELEM(angles,0,281) =131.8103;   MAT_ELEM(angles,1,281) =61.875;
    		MAT_ELEM(angles,0,282) =131.8103;   MAT_ELEM(angles,1,282) =50.625;
    		MAT_ELEM(angles,0,283) =125.6853;   MAT_ELEM(angles,1,283) =56.25;
    		MAT_ELEM(angles,0,284) =125.6853;   MAT_ELEM(angles,1,284) =67.5;
    		MAT_ELEM(angles,0,285) =120;   MAT_ELEM(angles,1,285) =73.125;
    		MAT_ELEM(angles,0,286) =120;   MAT_ELEM(angles,1,286) =61.875;
    		MAT_ELEM(angles,0,287) =114.6243;   MAT_ELEM(angles,1,287) =67.5;
    		MAT_ELEM(angles,0,288) =150.4344;   MAT_ELEM(angles,1,288) =9;
    		MAT_ELEM(angles,0,289) =144.3409;   MAT_ELEM(angles,1,289) =22.5;
    		MAT_ELEM(angles,0,290) =144.3409;   MAT_ELEM(angles,1,290) =7.5;
    		MAT_ELEM(angles,0,291) =138.1412;   MAT_ELEM(angles,1,291) =19.2857;
    		MAT_ELEM(angles,0,292) =138.1412;   MAT_ELEM(angles,1,292) =32.1429;
    		MAT_ELEM(angles,0,293) =131.8103;   MAT_ELEM(angles,1,293) =39.375;
    		MAT_ELEM(angles,0,294) =131.8103;   MAT_ELEM(angles,1,294) =28.125;
    		MAT_ELEM(angles,0,295) =125.6853;   MAT_ELEM(angles,1,295) =33.75;
    		MAT_ELEM(angles,0,296) =138.1412;   MAT_ELEM(angles,1,296) =6.4286;
    		MAT_ELEM(angles,0,297) =131.8103;   MAT_ELEM(angles,1,297) =16.875;
    		MAT_ELEM(angles,0,298) =131.8103;   MAT_ELEM(angles,1,298) =5.625;
    		MAT_ELEM(angles,0,299) =125.6853;   MAT_ELEM(angles,1,299) =11.25;
    		MAT_ELEM(angles,0,300) =125.6853;   MAT_ELEM(angles,1,300) =22.5;
    		MAT_ELEM(angles,0,301) =120;   MAT_ELEM(angles,1,301) =28.125;
    		MAT_ELEM(angles,0,302) =120;   MAT_ELEM(angles,1,302) =16.875;
    		MAT_ELEM(angles,0,303) =114.6243;   MAT_ELEM(angles,1,303) =22.5;
    		MAT_ELEM(angles,0,304) =125.6853;   MAT_ELEM(angles,1,304) =45;
    		MAT_ELEM(angles,0,305) =120;   MAT_ELEM(angles,1,305) =50.625;
    		MAT_ELEM(angles,0,306) =120;   MAT_ELEM(angles,1,306) =39.375;
    		MAT_ELEM(angles,0,307) =114.6243;   MAT_ELEM(angles,1,307) =45;
    		MAT_ELEM(angles,0,308) =114.6243;   MAT_ELEM(angles,1,308) =56.25;
    		MAT_ELEM(angles,0,309) =109.4712;   MAT_ELEM(angles,1,309) =61.875;
    		MAT_ELEM(angles,0,310) =109.4712;   MAT_ELEM(angles,1,310) =50.625;
    		MAT_ELEM(angles,0,311) =104.4775;   MAT_ELEM(angles,1,311) =56.25;
    		MAT_ELEM(angles,0,312) =114.6243;   MAT_ELEM(angles,1,312) =33.75;
    		MAT_ELEM(angles,0,313) =109.4712;   MAT_ELEM(angles,1,313) =39.375;
    		MAT_ELEM(angles,0,314) =109.4712;   MAT_ELEM(angles,1,314) =28.125;
    		MAT_ELEM(angles,0,315) =104.4775;   MAT_ELEM(angles,1,315) =33.75;
    		MAT_ELEM(angles,0,316) =104.4775;   MAT_ELEM(angles,1,316) =45;
    		MAT_ELEM(angles,0,317) =99.5941;   MAT_ELEM(angles,1,317) =50.625;
    		MAT_ELEM(angles,0,318) =99.5941;   MAT_ELEM(angles,1,318) =39.375;
    		MAT_ELEM(angles,0,319) =94.7802;   MAT_ELEM(angles,1,319) =45;
    		MAT_ELEM(angles,0,320) =174.1497;   MAT_ELEM(angles,1,320) =-45;
    		MAT_ELEM(angles,0,321) =168.2841;   MAT_ELEM(angles,1,321) =-22.5;
    		MAT_ELEM(angles,0,322) =168.2841;   MAT_ELEM(angles,1,322) =-67.5;
    		MAT_ELEM(angles,0,323) =162.3876;   MAT_ELEM(angles,1,323) =-45;
    		MAT_ELEM(angles,0,324) =162.3876;   MAT_ELEM(angles,1,324) =-15;
    		MAT_ELEM(angles,0,325) =156.4435;   MAT_ELEM(angles,1,325) =-11.25;
    		MAT_ELEM(angles,0,326) =156.4435;   MAT_ELEM(angles,1,326) =-33.75;
    		MAT_ELEM(angles,0,327) =150.4344;   MAT_ELEM(angles,1,327) =-27;
    		MAT_ELEM(angles,0,328) =162.3876;   MAT_ELEM(angles,1,328) =-75;
    		MAT_ELEM(angles,0,329) =156.4435;   MAT_ELEM(angles,1,329) =-56.25;
    		MAT_ELEM(angles,0,330) =156.4435;   MAT_ELEM(angles,1,330) =-78.75;
    		MAT_ELEM(angles,0,331) =150.4344;   MAT_ELEM(angles,1,331) =-63;
    		MAT_ELEM(angles,0,332) =150.4344;   MAT_ELEM(angles,1,332) =-45;
    		MAT_ELEM(angles,0,333) =144.3409;   MAT_ELEM(angles,1,333) =-37.5;
    		MAT_ELEM(angles,0,334) =144.3409;   MAT_ELEM(angles,1,334) =-52.5;
    		MAT_ELEM(angles,0,335) =138.1412;   MAT_ELEM(angles,1,335) =-45;
    		MAT_ELEM(angles,0,336) =150.4344;   MAT_ELEM(angles,1,336) =-9;
    		MAT_ELEM(angles,0,337) =144.3409;   MAT_ELEM(angles,1,337) =-7.5;
    		MAT_ELEM(angles,0,338) =144.3409;   MAT_ELEM(angles,1,338) =-22.5;
    		MAT_ELEM(angles,0,339) =138.1412;   MAT_ELEM(angles,1,339) =-19.2857;
    		MAT_ELEM(angles,0,340) =138.1412;   MAT_ELEM(angles,1,340) =-6.4286;
    		MAT_ELEM(angles,0,341) =131.8103;   MAT_ELEM(angles,1,341) =-5.625;
    		MAT_ELEM(angles,0,342) =131.8103;   MAT_ELEM(angles,1,342) =-16.875;
    		MAT_ELEM(angles,0,343) =125.6853;   MAT_ELEM(angles,1,343) =-11.25;
    		MAT_ELEM(angles,0,344) =138.1412;   MAT_ELEM(angles,1,344) =-32.1429;
    		MAT_ELEM(angles,0,345) =131.8103;   MAT_ELEM(angles,1,345) =-28.125;
    		MAT_ELEM(angles,0,346) =131.8103;   MAT_ELEM(angles,1,346) =-39.375;
    		MAT_ELEM(angles,0,347) =125.6853;   MAT_ELEM(angles,1,347) =-33.75;
    		MAT_ELEM(angles,0,348) =125.6853;   MAT_ELEM(angles,1,348) =-22.5;
    		MAT_ELEM(angles,0,349) =120;   MAT_ELEM(angles,1,349) =-16.875;
    		MAT_ELEM(angles,0,350) =120;   MAT_ELEM(angles,1,350) =-28.125;
    		MAT_ELEM(angles,0,351) =114.6243;   MAT_ELEM(angles,1,351) =-22.5;
    		MAT_ELEM(angles,0,352) =150.4344;   MAT_ELEM(angles,1,352) =-81;
    		MAT_ELEM(angles,0,353) =144.3409;   MAT_ELEM(angles,1,353) =-67.5;
    		MAT_ELEM(angles,0,354) =144.3409;   MAT_ELEM(angles,1,354) =-82.5;
    		MAT_ELEM(angles,0,355) =138.1412;   MAT_ELEM(angles,1,355) =-70.7143;
    		MAT_ELEM(angles,0,356) =138.1412;   MAT_ELEM(angles,1,356) =-57.8571;
    		MAT_ELEM(angles,0,357) =131.8103;   MAT_ELEM(angles,1,357) =-50.625;
    		MAT_ELEM(angles,0,358) =131.8103;   MAT_ELEM(angles,1,358) =-61.875;
    		MAT_ELEM(angles,0,359) =125.6853;   MAT_ELEM(angles,1,359) =-56.25;
    		MAT_ELEM(angles,0,360) =138.1412;   MAT_ELEM(angles,1,360) =-83.5714;
    		MAT_ELEM(angles,0,361) =131.8103;   MAT_ELEM(angles,1,361) =-73.125;
    		MAT_ELEM(angles,0,362) =131.8103;   MAT_ELEM(angles,1,362) =-84.375;
    		MAT_ELEM(angles,0,363) =125.6853;   MAT_ELEM(angles,1,363) =-78.75;
    		MAT_ELEM(angles,0,364) =125.6853;   MAT_ELEM(angles,1,364) =-67.5;
    		MAT_ELEM(angles,0,365) =120;   MAT_ELEM(angles,1,365) =-61.875;
    		MAT_ELEM(angles,0,366) =120;   MAT_ELEM(angles,1,366) =-73.125;
    		MAT_ELEM(angles,0,367) =114.6243;   MAT_ELEM(angles,1,367) =-67.5;
    		MAT_ELEM(angles,0,368) =125.6853;   MAT_ELEM(angles,1,368) =-45;
    		MAT_ELEM(angles,0,369) =120;   MAT_ELEM(angles,1,369) =-39.375;
    		MAT_ELEM(angles,0,370) =120;   MAT_ELEM(angles,1,370) =-50.625;
    		MAT_ELEM(angles,0,371) =114.6243;   MAT_ELEM(angles,1,371) =-45;
    		MAT_ELEM(angles,0,372) =114.6243;   MAT_ELEM(angles,1,372) =-33.75;
    		MAT_ELEM(angles,0,373) =109.4712;   MAT_ELEM(angles,1,373) =-28.125;
    		MAT_ELEM(angles,0,374) =109.4712;   MAT_ELEM(angles,1,374) =-39.375;
    		MAT_ELEM(angles,0,375) =104.4775;   MAT_ELEM(angles,1,375) =-33.75;
    		MAT_ELEM(angles,0,376) =114.6243;   MAT_ELEM(angles,1,376) =-56.25;
    		MAT_ELEM(angles,0,377) =109.4712;   MAT_ELEM(angles,1,377) =-50.625;
    		MAT_ELEM(angles,0,378) =109.4712;   MAT_ELEM(angles,1,378) =-61.875;
    		MAT_ELEM(angles,0,379) =104.4775;   MAT_ELEM(angles,1,379) =-56.25;
    		MAT_ELEM(angles,0,380) =104.4775;   MAT_ELEM(angles,1,380) =-45;
    		MAT_ELEM(angles,0,381) =99.5941;   MAT_ELEM(angles,1,381) =-39.375;
    		MAT_ELEM(angles,0,382) =99.5941;   MAT_ELEM(angles,1,382) =-50.625;
    		MAT_ELEM(angles,0,383) =94.7802;   MAT_ELEM(angles,1,383) =-45;
    	}
    	else
    	{
    	//TODO: use coordinates on the sphere instead of angles
    	angles.initZeros(2,81);
    	MAT_ELEM(angles, 0, 0) = 0.000000;	 	 MAT_ELEM(angles, 1, 0) = 0.000000;
    	MAT_ELEM(angles, 0, 1) = 36.000000;	 	 MAT_ELEM(angles, 1, 1) = 15.858741;
    	MAT_ELEM(angles, 0, 2) = 36.000000;	 	 MAT_ELEM(angles, 1, 2) = 31.717482;
    	MAT_ELEM(angles, 0, 3) = 36.000000;	 	 MAT_ELEM(angles, 1, 3) = 47.576224;
    	MAT_ELEM(angles, 0, 4) = 36.000000;	 	 MAT_ELEM(angles, 1, 4) = 63.434965;
    	MAT_ELEM(angles, 0, 5) = 62.494295;	 	 MAT_ELEM(angles, 1, 5) = -76.558393;
    	MAT_ELEM(angles, 0, 6) = 54.000000;	 	 MAT_ELEM(angles, 1, 6) = 90.000000;
    	MAT_ELEM(angles, 0, 7) = 45.505705;	 	 MAT_ELEM(angles, 1, 7) = 76.558393;
    	MAT_ELEM(angles, 0, 8) = 108.000000;	 MAT_ELEM(angles, 1, 8) = 15.858741;
    	MAT_ELEM(angles, 0, 9) = 108.000000;	 MAT_ELEM(angles, 1, 9) = 31.717482;
    	MAT_ELEM(angles, 0, 10) = 108.000000;	 MAT_ELEM(angles, 1, 10) = 47.576224;
    	MAT_ELEM(angles, 0, 11) = 108.000000;	 MAT_ELEM(angles, 1, 11) = 63.434965;
    	MAT_ELEM(angles, 0, 12) = 134.494295;	 MAT_ELEM(angles, 1, 12) = -76.558393;
    	MAT_ELEM(angles, 0, 13) = 126.000000;	 MAT_ELEM(angles, 1, 13) = 90.000000;
    	MAT_ELEM(angles, 0, 14) = 117.505705;	 MAT_ELEM(angles, 1, 14) = 76.558393;
    	MAT_ELEM(angles, 0, 15) = 144.000000;	 MAT_ELEM(angles, 1, 15) = -15.858741;
    	MAT_ELEM(angles, 0, 16) = 144.000000;	 MAT_ELEM(angles, 1, 16) = -31.717482;
    	MAT_ELEM(angles, 0, 17) = 144.000000;	 MAT_ELEM(angles, 1, 17) = -47.576224;
    	MAT_ELEM(angles, 0, 18) = 144.000000;	 MAT_ELEM(angles, 1, 18) = -63.434965;
    	MAT_ELEM(angles, 0, 19) = 170.494295;	 MAT_ELEM(angles, 1, 19) = 76.558393;
    	MAT_ELEM(angles, 0, 20) = 162.000000;	 MAT_ELEM(angles, 1, 20) = 90.000000;
    	MAT_ELEM(angles, 0, 21) = 153.505705;	 MAT_ELEM(angles, 1, 21) = -76.558393;
    	MAT_ELEM(angles, 0, 22) = 72.000000;	 MAT_ELEM(angles, 1, 22) = -15.858741;
    	MAT_ELEM(angles, 0, 23) = 72.000000;	 MAT_ELEM(angles, 1, 23) = -31.717482;
    	MAT_ELEM(angles, 0, 24) = 72.000000;	 MAT_ELEM(angles, 1, 24) = -47.576224;
    	MAT_ELEM(angles, 0, 25) = 72.000000;	 MAT_ELEM(angles, 1, 25) = -63.434965;
    	MAT_ELEM(angles, 0, 26) = 98.494295;	 MAT_ELEM(angles, 1, 26) = 76.558393;
    	MAT_ELEM(angles, 0, 27) = 90.000000;	 MAT_ELEM(angles, 1, 27) = 90.000000;
    	MAT_ELEM(angles, 0, 28) = 81.505705;	 MAT_ELEM(angles, 1, 28) = -76.558393;
    	MAT_ELEM(angles, 0, 29) = 0.000000;	 	 MAT_ELEM(angles, 1, 29) = -15.858741;
    	MAT_ELEM(angles, 0, 30) = 0.000000;	 	 MAT_ELEM(angles, 1, 30) = -31.717482;
    	MAT_ELEM(angles, 0, 31) = 0.000000;	 	 MAT_ELEM(angles, 1, 31) = -47.576224;
    	MAT_ELEM(angles, 0, 32) = 0.000000;	 	 MAT_ELEM(angles, 1, 32) = -63.434965;
    	MAT_ELEM(angles, 0, 33) = 26.494295;	 MAT_ELEM(angles, 1, 33) = 76.558393;
    	MAT_ELEM(angles, 0, 34) = 18.000000;	 MAT_ELEM(angles, 1, 34) = 90.000000;
    	MAT_ELEM(angles, 0, 35) = 9.505705;	 	 MAT_ELEM(angles, 1, 35) = -76.558393;
    	MAT_ELEM(angles, 0, 36) = 12.811021;	 MAT_ELEM(angles, 1, 36) = 42.234673;
    	MAT_ELEM(angles, 0, 37) = 18.466996;	 MAT_ELEM(angles, 1, 37) = 59.620797;
    	MAT_ELEM(angles, 0, 38) = 0.000000;	 	 MAT_ELEM(angles, 1, 38) = 90.000000;
    	MAT_ELEM(angles, 0, 39) = 8.867209;	 	 MAT_ELEM(angles, 1, 39) = 75.219088;
    	MAT_ELEM(angles, 0, 40) = 72.000000;	 MAT_ELEM(angles, 1, 40) = 26.565058;
    	MAT_ELEM(angles, 0, 41) = 59.188979;	 MAT_ELEM(angles, 1, 41) = 42.234673;
    	MAT_ELEM(angles, 0, 42) = 84.811021;	 MAT_ELEM(angles, 1, 42) = 42.234673;
    	MAT_ELEM(angles, 0, 43) = 53.533003;	 MAT_ELEM(angles, 1, 43) = 59.620797;
    	MAT_ELEM(angles, 0, 44) = 72.000000;	 MAT_ELEM(angles, 1, 44) = 58.282544;
    	MAT_ELEM(angles, 0, 45) = 90.466996;	 MAT_ELEM(angles, 1, 45) = 59.620797;
    	MAT_ELEM(angles, 0, 46) = 72.000000;	 MAT_ELEM(angles, 1, 46) = 90.000000;
    	MAT_ELEM(angles, 0, 47) = 63.132791;	 MAT_ELEM(angles, 1, 47) = 75.219088;
    	MAT_ELEM(angles, 0, 48) = 80.867209;	 MAT_ELEM(angles, 1, 48) = 75.219088;
    	MAT_ELEM(angles, 0, 49) = 144.000000;	 MAT_ELEM(angles, 1, 49) = 26.565058;
    	MAT_ELEM(angles, 0, 50) = 131.188979;	 MAT_ELEM(angles, 1, 50) = 42.234673;
    	MAT_ELEM(angles, 0, 51) = 156.811021;	 MAT_ELEM(angles, 1, 51) = 42.234673;
    	MAT_ELEM(angles, 0, 52) = 125.533003;	 MAT_ELEM(angles, 1, 52) = 59.620797;
    	MAT_ELEM(angles, 0, 53) = 144.000000;	 MAT_ELEM(angles, 1, 53) = 58.282544;
    	MAT_ELEM(angles, 0, 54) = 162.466996;	 MAT_ELEM(angles, 1, 54) = 59.620797;
    	MAT_ELEM(angles, 0, 55) = 144.000000;	 MAT_ELEM(angles, 1, 55) = 90.000000;
    	MAT_ELEM(angles, 0, 56) = 135.132791;	 MAT_ELEM(angles, 1, 56) = 75.219088;
    	MAT_ELEM(angles, 0, 57) = 152.867209;	 MAT_ELEM(angles, 1, 57) = 75.219088;
    	MAT_ELEM(angles, 0, 58) = 180.000000;	 MAT_ELEM(angles, 1, 58) = -26.565058;
    	MAT_ELEM(angles, 0, 59) = 167.188979;	 MAT_ELEM(angles, 1, 59) = -42.234673;
    	MAT_ELEM(angles, 0, 60) = 180.000000;	 MAT_ELEM(angles, 1, 60) = -58.282544;
    	MAT_ELEM(angles, 0, 61) = 161.533003;	 MAT_ELEM(angles, 1, 61) = -59.620797;
    	MAT_ELEM(angles, 0, 62) = 171.132791;	 MAT_ELEM(angles, 1, 62) = -75.219088;
    	MAT_ELEM(angles, 0, 63) = 108.000000;	 MAT_ELEM(angles, 1, 63) = -26.565058;
    	MAT_ELEM(angles, 0, 64) = 120.811021;	 MAT_ELEM(angles, 1, 64) = -42.234673;
    	MAT_ELEM(angles, 0, 65) = 95.188979;	 MAT_ELEM(angles, 1, 65) = -42.234673;
    	MAT_ELEM(angles, 0, 66) = 126.466996;	 MAT_ELEM(angles, 1, 66) = -59.620797;
    	MAT_ELEM(angles, 0, 67) = 108.000000;	 MAT_ELEM(angles, 1, 67) = -58.282544;
    	MAT_ELEM(angles, 0, 68) = 89.533003;	 MAT_ELEM(angles, 1, 68) = -59.620797;
    	MAT_ELEM(angles, 0, 69) = 108.000000;	 MAT_ELEM(angles, 1, 69) = 90.000000;
    	MAT_ELEM(angles, 0, 70) = 116.867209;	 MAT_ELEM(angles, 1, 70) = -75.219088;
    	MAT_ELEM(angles, 0, 71) = 99.132791;	 MAT_ELEM(angles, 1, 71) = -75.219088;
    	MAT_ELEM(angles, 0, 72) = 36.000000;	 MAT_ELEM(angles, 1, 72) = -26.565058;
    	MAT_ELEM(angles, 0, 73) = 48.811021;	 MAT_ELEM(angles, 1, 73) = -42.234673;
    	MAT_ELEM(angles, 0, 74) = 23.188979;	 MAT_ELEM(angles, 1, 74) = -42.234673;
    	MAT_ELEM(angles, 0, 75) = 54.466996;	 MAT_ELEM(angles, 1, 75) = -59.620797;
    	MAT_ELEM(angles, 0, 76) = 36.000000;	 MAT_ELEM(angles, 1, 76) = -58.282544;
    	MAT_ELEM(angles, 0, 77) = 17.533003;	 MAT_ELEM(angles, 1, 77) = -59.620797;
    	MAT_ELEM(angles, 0, 78) = 36.000000;	 MAT_ELEM(angles, 1, 78) = 90.000000;
    	MAT_ELEM(angles, 0, 79) = 44.867209;	 MAT_ELEM(angles, 1, 79) = -75.219088;
    	MAT_ELEM(angles, 0, 80) = 27.132791;	 MAT_ELEM(angles, 1, 80) = -75.219088;
    	}

    }


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
			double rot, double tilt, double ang_con)
    {
		int ZdimFT1=(int)ZSIZE(threeD_FSC);
		int YdimFT1=(int)YSIZE(threeD_FSC);
		int XdimFT1=(int)XSIZE(threeD_FSC);

		double maxFreq_2 =0.;
		maxFreq_2 = maxFreq;
		double x_dir, y_dir, z_dir, uz, uy, ux, cosAngle;
		cosAngle = cos(ang_con);
		x_dir = sin(tilt*PI/180)*cos(rot*PI/180);
		y_dir = sin(tilt*PI/180)*sin(rot*PI/180);
		z_dir = cos(tilt*PI/180);
		long n = 0;
		for (int k=0; k<ZdimFT1; k++)
		{
			double uz = VEC_ELEM(freq_fourier_z,k);
			uz *= z_dir;
			for (int i=0; i<YdimFT1; i++)
			{
				double uy = VEC_ELEM(freq_fourier_y,i);
				uy *= y_dir;
				for (int j=0; j<XdimFT1; j++)
				{
					double ux = VEC_ELEM(freq_fourier_x,j);
					ux *= x_dir;
					double iun = DIRECT_MULTIDIM_ELEM(freqMap, n);
					double f = 1/iun;
					iun *= (ux + uy + uz);
					double cosine =fabs(iun);

					if (cosine>=cosAngle)
						{
							if (f>maxFreq_2)
							{
								++n;
								continue;
							}
							int idx = (int) round(f * m1sizeX);
							cosine = exp( -((cosine -1)*(cosine -1))/pow(ang_con,6) );
							DIRECT_MULTIDIM_ELEM(threeD_FSC, n) += cosine*dAi(fsc, idx);
							DIRECT_MULTIDIM_ELEM(counterMap, n) += cosine;//1.0;

						}
					++n;
				}
			}
		}
    }


    void sortArr(double arr[], int n, std::vector<std::pair<double, int> > &vp)
    {
        // Inserting element in pair vector
        // to keep track of previous indexes
        for (int i = 0; i < n; ++i) {
            vp.push_back(std::make_pair(arr[i], i));
        }

        // Sorting pair vector
        sort(vp.begin(), vp.end());

        // Displaying sorted element
        // with previous indexes
        // corresponding to each element
        std::cout << "Element\t"
             << "index" << std::endl;
        for (int i = 0; i < vp.size(); i++) {
        	std::cout << vp[i].first << "\t"
                 << vp[i].second << std::endl;
        }
    }


    void anistropyParameter(const MultidimArray<double> FSC,
    		const MultidimArray<double> &freq,
			MultidimArray<double> &aniParam)
    {

    	double thrs = 0.143;
		for (size_t k = 0; k<aniParam.nzyxdim; k++)
			if (DIRECT_MULTIDIM_ELEM(FSC, k) >= thrs)
				DIRECT_MULTIDIM_ELEM(aniParam, k) += 1.0;
    }


    void prepareData(FileName &fnhalf1, FileName &fnhalf2,
    		MultidimArray<double> &half1, MultidimArray<double> &half2, bool test)
    {
    	MultidimArray<double> &phalf1 = half1, &phalf2 = half2;

    	Image<double> mask;
    	MultidimArray<double> &pmask = mask();

    	if (test)
		{
			Monogenic mono;
			std::cout << "Preparing test data ..." << std::endl;
			size_t xdim = 301, ydim = 301, zdim = 301;
			double wavelength = 5.0, mean = 0.0, std = 0.5;
			int maskrad = 125;
			half1 = mono.createDataTest(xdim, ydim, zdim, wavelength, mean, 0.0);
			half2 = half1;

			mono.addNoise(phalf1, 0, std);
			mono.addNoise(phalf2, 0, std);
			FileName fn;
			Image<double> saveImg;
			fn = formatString("inputVol1_large.vol");
			saveImg() = half1;
			saveImg.write(fn);
			fn = formatString("inputVol2_large.vol");
			saveImg() = half2;
			saveImg.write(fn);
		}
		else
		{
			std::cout << "Reading data..." << std::endl;
			Image<double> imgHalf1, imgHalf2;
			imgHalf1.read(fnhalf1);
			imgHalf2.read(fnhalf2);

			half1 = imgHalf1();
			half2 = imgHalf2();

			if (fnmask!="")
			{
				mask.read(fnmask);
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pmask)
				{
					double valmask = (double) DIRECT_MULTIDIM_ELEM(pmask, n);
					DIRECT_MULTIDIM_ELEM(phalf1, n) = DIRECT_MULTIDIM_ELEM(phalf1, n) * valmask;
					DIRECT_MULTIDIM_ELEM(phalf2, n) = DIRECT_MULTIDIM_ELEM(phalf2, n) * valmask;
				}
			}
			mask.clear();
			pmask.clear();
		}

		phalf1.setXmippOrigin();
		phalf2.setXmippOrigin();

    	std::cout << "Starting..." << std::endl;
    }

    void saveFourierAmplitudes(FourierTransformer &transformer1, const MultidimArray< std::complex< double > > &FT1)
    {
    	MultidimArray<double> FT1_ampl;
		FT1_ampl.initZeros(FT1);
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(FT1)
			DIRECT_MULTIDIM_ELEM(FT1_ampl, n) = abs(DIRECT_MULTIDIM_ELEM(FT1, n));

		std::complex<double> J(0,1);
		transformer1.fFourier.initZeros(FT1_ampl);
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(transformer1.fFourier)
		{
			DIRECT_MULTIDIM_ELEM(transformer1.fFourier, n) = -J*J*DIRECT_MULTIDIM_ELEM(FT1_ampl, n);
		}

		transformer1.getCompleteFourier(FT1_ampl);
		CenterFFT(FT1_ampl, true);

		Image<double> saveImgFT;
		saveImgFT() = FT1_ampl;
		saveImgFT.write("FT.vol");
    }

    void run()
    {
    	MultidimArray<double> half1, half2;
    	MultidimArray<double> &phalf1 = half1, &phalf2 = half2;

    	prepareData(fnhalf1, fnhalf2, half1, half2, test);

    	unsigned t0, t1;
    	t0=clock();


		//Defining Fourier transform
    	MultidimArray< std::complex< double > > FT1, FT2;

        FourierTransformer transformer2(FFTW_BACKWARD), transformer1(FFTW_BACKWARD);
        transformer1.FourierTransform(phalf1, FT1, false);
        transformer2.FourierTransform(phalf2, FT2, false);
        //This function modifies the Fourier transform
        //saveFourierAmplitudes(transformer1, FT1);

        int m1sizeX = XSIZE(phalf1), m1sizeY = YSIZE(phalf1), m1sizeZ = ZSIZE(phalf1);

        //Defining frequencies
        Matrix1D<double> freq_fourier_x, freq_fourier_y, freq_fourier_z;
        MultidimArray<double> freqMap;

        freqMap = defineFrequencies(FT1, phalf1,
        		freq_fourier_x,freq_fourier_y, freq_fourier_z);

        t1 = clock();

        double time = (double(t1-t0)/CLOCKS_PER_SEC);
        std::cout << "%Execution Time: " << time << std::endl;

        //TODO: check when they can be cleared
//        phalf2.clear(); // Free memory
//        phalf1.clear(); // Free memory

    	MultidimArray<double> fsc, freq, counterMap, threeD_FSC, aniParam;
    	counterMap.resizeNoCopy(FT1);
    	threeD_FSC.resizeNoCopy(counterMap);
    	counterMap.initConstant(1e-38);

    	std::cout << "ang_con = " << ang_con << std::endl;

    	ang_con = 20*PI/180;
    	bool alot = true;

    	generateDirections(angles, true);

    	std::cout << "angles.mdim = " << angles.mdim << std::endl;

    	aniParam.initZeros(m1sizeX/2+1);
    	for (size_t k = 0; k<angles.mdimx; k++)
			{
			double rot = MAT_ELEM(angles, 0, k);
			double tilt = MAT_ELEM(angles, 1, k);
			std::cout << "%rot " << rot << "  " << tilt << std::endl;
			fscDir(FT1, FT2, sampling, freq_fourier_x, freq_fourier_y, freq_fourier_z, freqMap, freq, fsc, 0.5,
					m1sizeX, m1sizeY, m1sizeZ, rot, tilt, ang_con);
			std::cout << "%------------------------------" <<  std::endl;

			std::cout << "B_" << k << "=[" << std::endl;
			FOR_ALL_ELEMENTS_IN_ARRAY1D(fsc)
				std::cout << dAi(freq, i) << "     " << dAi(fsc,i) << ";" << std::endl;
			std::cout << "];" << std::endl;

			anistropyParameter(fsc, freq, aniParam);


			interpolationCoarse(fsc, angles,
					freq_fourier_x, freq_fourier_y, freq_fourier_z,
		    		threeD_FSC, counterMap,
					freqMap, freq,
					0.5, m1sizeX, m1sizeY, m1sizeZ,
					rot, tilt, ang_con);
    	}
    	std::cout << "%------------------------------" <<  std::endl;
    	aniParam /= (double) angles.mdimx;

    	std::cout << "A_" << std::endl;
    	FOR_ALL_ELEMENTS_IN_ARRAY1D(aniParam)
    		std::cout << dAi(freq, i) << " " << dAi(aniParam, i) << ";" << std::endl;

    	std::cout << "];" << std::endl;


        unsigned t2 = clock();

        time = (double(t2-t0)/CLOCKS_PER_SEC);
        std::cout << "%Execution Time: " << time << std::endl;

        std::complex<double> J(0,1);
    	transformer1.fFourier.initZeros(threeD_FSC);
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(transformer1.fFourier)
		{
			DIRECT_MULTIDIM_ELEM(threeD_FSC, n) /= DIRECT_MULTIDIM_ELEM(counterMap, n);
			DIRECT_MULTIDIM_ELEM(transformer1.fFourier, n) = -J*J*DIRECT_MULTIDIM_ELEM(threeD_FSC, n);
		}

		MultidimArray<double> fscFull;
		transformer1.getCompleteFourier(fscFull);
		std::cout << "%dim = " << fscFull.getDim() <<std::endl;
//		CenterFFT(fscFull, true);
		Image<double> saveImg2;
		saveImg2() = fscFull;
		saveImg2.write(fn_out);

    }

};
