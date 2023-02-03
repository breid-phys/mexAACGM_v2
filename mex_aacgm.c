#include "mex.h"
#include "matrix.h"
#include "aacgmlib_v2.h"
#include <math.h>


#define sind(x) (sin(fmod((x),360) * M_PI / 180))
#define cosd(x) (cos(fmod((x),360) * M_PI / 180))
#define asind(x) (asin(x) * 180 / M_PI)
#define acosd(x) (acos(x) * 180 / M_PI)
#define atan2d(x,y) (atan2(x,y) * 180 / M_PI)
#define fix(x)  (double)(long)(x)

double subsolar(int yr, int mo, int dy, int hr, int mt, int sc) {
        double Y, M, D, A, B, JD, T, t0, GMST, L0, ep, C, v, R, Om, Alon, e0, RA;
        double Slon, Slat, err;
        double mlat, mlon, r;

        Y = (double)yr;
        M = (double)mo;
        D = (double)dy + (double)hr/24.0 + (double)mt/(24.0*60.0) + (double)sc/(24.0*60.0*60.0);

        if (Y<=2) {
            Y = Y - 1.0;
            M = M + 12.0;
        }
        
        A = fix(Y/100.0);
        B = 2.0 - A + fix(A/4.0);
        JD = fix(365.25*(Y+4716.0))+fix(30.6001*(M+1.0))+D+B-1524.5;

        T = (JD-2451545.0)/36525;

        t0 = 280.46061837 + 360.98564736629*(JD-2451545.0) + 0.000387933*pow(T,2) + (-1.0/38710000.0)*pow(T,3);

        GMST = fmod(t0,360.0);

        L0 = 280.46646 + 36000.76983*T + 0.0003032*pow(T,2);
        M = 357.52911+35999.05029*T - 0.0001537*pow(T,2);
        ep = 0.016708634 - 0.000042037*T - 0.0000001267*pow(T,2);

        C = (1.914602 - 0.004817*T - 0.000014*pow(T,2))*sind(M) + (0.019993 - 0.000101*T)*sind(2*M) + (0.000289)*sind(3*M);

        Slon = L0 + C;
        v = M + C;

        R = (1.000001018*(1-pow(ep,2)))/(1.0+ep*cosd(v));

        Om = 125.04 - 1934.136*T;
        Alon = Slon - 0.00569 - 0.00478*sind(Om);
        e0 = (23.0+(26.0/60.0)+(21.448/3600.0)) + (-46.8150/3600.0)*T + (-0.00059/3600.0)*pow(T,2) + (0.001813/3600.0)*pow(T,3);

        e0 = e0 + 0.00256*cosd(Om);

        RA = atan2d(cosd(e0)*sind(Alon),cosd(Alon));
        Slat = asind(sind(e0)*sind(Alon));

        Slon = fmod(RA-GMST,360.0);

        err = AACGM_v2_Convert(Slat, Slon, 700.0, &mlat, &mlon, &r, G2A);

        return mlon;
}
/*
 * xtimesy.c
 * Multiplies a real input matrix y by a given real input scalar x.
 *
 * The calling syntax is:
 *
 *		[result] = xtimesy(x, y)
 *
 * This is a MEX-file for MATLAB.
 * Copyright 1984-2018 The MathWorks, Inc.
 */

/* the gateway function */
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {

    double lat, lon, hgt;
    double *lla;

    int yr, mo, dy, hr, mt, sc;
    int prev_yr, prev_mo, prev_dy, prev_hr, prev_mt, prev_sc;
    double *datevec;

    double *lla_out;
    double *mlt_out;

    double *geo2mag;
    int code;
    int N_MAX_WARN; // maximum number of warnings to print
    int N_WARN;

    double mlt;
    double mlon_ref;
    double rtp[3];
    double mlat,mlon,r;
    int k, err, npts;
    double dd,jd,eqt,dec,ut,at;

    
    char *input_buf;
    size_t buflen;

    size_t mrows, ncols;

    N_MAX_WARN = 50;
    N_WARN = 0;

    /*  check for proper number of arguments */
    if (nrhs != 5)
        mexErrMsgIdAndTxt("AACGM_v2:InputError:invalidNumInputs", "Five inputs required.");
    if (nlhs > 2)
        mexErrMsgIdAndTxt("AACGM_v2:InputError:invalidNumOutputs", "One or two outputs required.");

    if ( mxIsChar(prhs[0]) != 1)
      mexErrMsgIdAndTxt( "AACGM_v2:InputError:inputNotString",
              "AACGM Prefix must be a string.");

     if ( mxIsChar(prhs[1]) != 1)
      mexErrMsgIdAndTxt( "AACGM_v2:InputError:inputNotString",
              "IGRF Coefficient path must be a string.");    

    /* get the length of the AACGM Prefix */
    buflen = (mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;

    /* copy the string data from prhs[0] into a C string input_ buf.    */
    input_buf = mxArrayToString(prhs[0]);

    setenv("AACGM_v2_DAT_PREFIX", input_buf, 1);

    /* get the length of the IGRF Coefficients */
    buflen = (mxGetM(prhs[1]) * mxGetN(prhs[1])) + 1;

    /* copy the string data from prhs[0] into a C string input_ buf.    */
    input_buf = mxArrayToString(prhs[1]);

    setenv("IGRF_COEFFS", input_buf, 1);

    mxFree(input_buf);

    if( !mxIsDouble(prhs[2]) || 
         mxIsComplex(prhs[2])) {
        mexErrMsgIdAndTxt("AACGM_v2:InputError:notDouble","Input lla must be type double.");
    }

    if(mxGetM(prhs[2])!=3) {
        mexErrMsgIdAndTxt("AACGM_v2:InputError:lla","Input lla must have 3 columns.");
    }
    
    if(mxGetN(prhs[2])!=mxGetN(prhs[3])) {
        mexErrMsgIdAndTxt("AACGM_v2:InputError:SizeMismatch","Input lla and date must have the same number of rows.");
    }
    
    if( !mxIsDouble(prhs[3]) || 
         mxIsComplex(prhs[3])) {
        mexErrMsgIdAndTxt("AACGM_v2:InputError:notInt","Input datevec must be type int32.");
    }

    if(mxGetM(prhs[3])!=6) {
        mexErrMsgIdAndTxt("AACGM_v2:InputError:date","Input datevec must have 6 columns.");
    }

#if MX_HAS_INTERLEAVED_COMPLEX
    lla = mxGetDoubles(prhs[2]);
#else
    lla = mxGetPr(prhs[2]);
#endif

#if MX_HAS_INTERLEAVED_COMPLEX
    datevec = mxGetDoubles(prhs[3]);
#else
    datevec = mxGetPr(prhs[3]);
#endif
    /*  get the dimensions of the matrix input y */
    ncols = mxGetN(prhs[2]);
    mrows = mxGetM(prhs[2]);

#if MX_HAS_INTERLEAVED_COMPLEX
    geo2mag = mxGetDoubles(prhs[4]);
#else
    geo2mag = mxGetPr(prhs[4]);
#endif

    code = (int)*(geo2mag);

/*  set the output pointer to the output matrix */
    plhs[0] = mxCreateDoubleMatrix((mwSize)3, (mwSize)ncols, mxREAL);

#if MX_HAS_INTERLEAVED_COMPLEX
    lla_out = mxGetDoubles(plhs[0]);
#else
    lla_out = mxGetPr(plhs[0]);
#endif

    if (nlhs > 1){
        plhs[1] = mxCreateDoubleMatrix((mwSize)1, (mwSize)ncols, mxREAL);
#if MX_HAS_INTERLEAVED_COMPLEX
        mlt_out = mxGetDoubles(plhs[1]);
#else
        mlt_out = mxGetPr(plhs[1]);
#endif
    }

    /* initialize prev times*/
    prev_yr = -1;
    prev_mo = -1;
    prev_dy = -1;
    prev_hr = -1;
    prev_mt = -1;
    prev_sc = -1;
    
    for (int i=0; i<ncols; i++) {
        
        yr = (int)*(datevec + i*6);
        mo = (int)*(datevec + i*6 + 1);
        dy = (int)*(datevec + i*6 + 2);
        hr = (int)*(datevec + i*6 + 3);
        mt = (int)*(datevec + i*6 + 4);
        sc = (int)*(datevec + i*6 + 5);

        if (sc!=prev_sc || mt!=prev_mt || hr!=prev_hr || dy!=prev_dy || mo!=prev_mo || yr!=prev_yr) {
            prev_yr = yr;
            prev_mo = mo;
            prev_dy = dy;
            prev_hr = hr;
            prev_mt = mt;
            prev_sc = sc; 
            AACGM_v2_SetDateTime(yr, mo, dy, hr, mt, sc);
            if (nlhs > 1){
                mlon_ref = subsolar(yr, mo, dy, hr, mt, sc);
            }
        }

        lat = *(lla + i*3);
        lon = *(lla + i*3 + 1);

        if (code & A2G) { // aacgm to geodetic
            hgt = 1.0 + (*(lla + i*3 + 2))/RE;
        } else { // geodetic to aacgm
            hgt = *(lla + i*3 + 2);
        }

        err = AACGM_v2_Convert(lat, lon, hgt, &mlat, &mlon, &r, code);

        if ((err < -1) && (N_WARN < N_MAX_WARN)) {
            N_WARN++;
            if (err == -2){
                mexWarnMsgIdAndTxt("AACGM_v2:ConversionError:HeightOutOfRange",
                "WARNING: coordinate transformations are not intended "
                "for altitudes < 0 km: %lf\n", hgt);
            } else if (err == -4) {
                mexWarnMsgIdAndTxt("AACGM_v2:ConversionError:HeightOutOfRange",
                "ERROR: coefficients are not valid for altitudes "
                    "above %d km: %lf.\n       You must either use field-line tracing "
                    "(TRACE or ALLOWTRACE) or\n"
                    "       indicate that you know this is a very bad idea "
                    "(BADIDEA)\n\n",                    
                    MAXALT, hgt);
            } else if (err == -8){
                mexWarnMsgIdAndTxt("AACGM_v2:ConversionError:LatOutOfRange",
                "ERROR: latitude must be in the range -90 to +90 degrees: "
                    "%lf\n", lat);
            } 
            else {
            mexWarnMsgIdAndTxt("AACGM_v2:ConversionError","aacgm_v2 reporting error code:%i", err);
            }
        }
      
        *(lla_out + i*3) = mlat;
        *(lla_out + i*3 + 1) = mlon;

        if (code & A2G) { // aacgm to geodetic
            *(lla_out + i*3 + 2) = r;
        } else { // geodetic to aacgm
            *(lla_out + i*3 + 2) = (r-1.)*RE;
        }
        
        if (nlhs > 1){
            if (code & A2G) { // aacgm to geodetic
                mlt = 12. + (lon - mlon_ref)/15.;  /* MLT based on subsolar point */
            } else { // geodetic to aacgm
                mlt = 12. + (mlon - mlon_ref)/15.;  /* MLT based on subsolar point */
            }

            if (mlt > 24.) mlt -= 24.;
            if (mlt <  0.) mlt += 24.;
            *(mlt_out + i) = mlt;
        }
        

    }

    

}

