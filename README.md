#SURF - Speeded Up Robust Features

##Requirements

The library has been compiled using g++, version 4.0.2, for usage on
a machine Pentium 4 or better. To use the library in your program,
you need to use the same compiler version.

If you require the library to be compiled using another compiler, or
another platform (such as Athlon XP), please contact us.

##Usage

Execute surf.ln without any argument in order to get more
information concerning the usage and possible parameters.

Use "make match.ln" to compile the matching demo application.

##Data Format

The output format of SURF is as follows:

 (1 + length of descriptor)
 number of points
 x y a b c l des
 x y a b c l des
 ...

 x, y = position of interest point
 a, b, c = [a b; b c] entries of second moment matrix.
   SURF only has circular regions, hence b = 0; a = c -> radius = 1 / a^2
 l = sign of laplacian (-1 or 1)
 des = descriptor vector itself

##Data Input Format

If only the SURF descriptor should be computed, the -p1 command can
be used. As an argument, it takes a file of the following format:

 (dummy byte)
 number of points
 x y a b c
 x y a b c

Where, as above, [a b; b c] forms the second moment matrix. Note that
SURF uses circular regions. Hence, a = c and b = 0.

##Licensing conditions

This software is being made available for research purposes only.  It
is necessary to obtain a license (see LICENSE file) for commercial
applications.
