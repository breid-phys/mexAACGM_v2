A MATLAB MEX wrapper for the Altitude Adjusted Corrected Geogmagnetic Coordinates

A MEX file allows users to execute C code from within MATLAB.

aacgm.m should automatically compile a MEX file when first run.


See aacgm.m for details.

http://superdarn.thayer.dartmouth.edu/aacgm.html

#######################################################################
AACGM-v2 Software
v2.6 20191228

This package include the following files:

AACGM C software:

README.txt            ; this file
release_notes.txt     ; details of changes to v2.6
aacgmlib_v2.c         ; AACGM-v2 functions
aacgmlib_v2.h         ; AACGM-v2 header file
genmag.c              ; general purpose functions
genmag.h              ; general purpose header file
igrflib.c             ; internal IGRF functions
igrflib.h             ; internal IGRF header file
rtime.c               ; internal date/time functions
rtime.h               ; internal date/time header file
astalg.c              ; Astronomical algorithms functions
astalg.h              ; Astronomical algorithms header file
mlt_v2.c              ; MLT-v2 functions
mlt_v2.h              ; MLT-v2 header file
igrf13coeffs.txt      ; IGRF13 coefficients (1900-2020)
magmodel_1590-2020.txt; magnetic field coefficients (1590-2020)
test_aacgm.c          ; testing and example program
LICENSE-AstAlg.txt    ; license file for Astro algrorithms

