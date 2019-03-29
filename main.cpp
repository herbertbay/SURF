/*
 * Speeded-Up Robust Features (SURF)
 * http://people.ee.ethz.ch/~surf
 *
 * Authors: Herbert Bay, Andreas Ess, Geert Willems
 * Windows port by Stefan Saur
 *
 * Copyright (2006): ETH Zurich, Switzerland
 * Katholieke Universiteit Leuven, Belgium
 * All rights reserved.
 *
 * For details, see the paper:
 * Herbert Bay,  Tinne Tuytelaars,  Luc Van Gool,
 *  "SURF: Speeded Up Robust Features"
 * Proceedings of the ninth European Conference on Computer Vision, May 2006
 *
 * Permission to use, copy, modify, and distribute this software and
 * its documentation for educational, research, and non-commercial
 * purposes, without fee and without a signed licensing agreement, is
 * hereby granted, provided that the above copyright notice and this
 * paragraph appear in all copies modifications, and distributions.
 *
 * Any commercial use or any redistribution of this software
 * requires a license from one of the above mentioned establishments.
 *
 * For further details, contact Herbert Bay (herbert.bay@gmail.com).
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string.h>

#ifdef WIN32
#include "surfWINDLL.h"
#endif

#include "imload.h"
#include "surflib.h"

#include "os_mapping.h"

using namespace std;
using namespace surf;

// Length of the descriptor vector
int VLength;

// Forward declaration of the functions to load/save the SURF points
void loadIpoints(string fn, vector< Ipoint >& keys, bool bVerbose = false);
void saveIpoints(string fn, const vector< Ipoint >& keys, bool bVerbose = false, bool bLaplacian = true);

int main (int argc, char **argv)
{
  // Initial sampling step (default 2)
  int samplingStep = 2;
  // Number of analysed octaves (default 4)
  int octaves = 4;
  // Blob response treshold
  double thres = 4.0;
  // Set this flag "true" to double the image size
  bool doubleImageSize = false;
  // Initial lobe size, default 3 and 5 (with double image size)
  int initLobe = 3;
  // Upright SURF or rotation invaraiant
  bool upright = false;
  // If the extended flag is turned on, SURF 128 is used
  bool extended = false;
  // Spatial size of the descriptor window (default 4)
  int indexSize = 4;
  // Variables for the timing measure
  osmapping::os_TIME tim1, tim2; //STS
  // verbose output
  bool bVerbose = true;
  // skip sign of laplacian
  bool bLaplacian = true;

  bool bLoadRegions  = false;
  string sRegionFile = "";

  // Print command line help
  if (argc==1) {
    cerr << "./surf -i img.pgm -o img.surf [options]\n"
         << "  blob response threshold          -thres 1000\n"
         << "  double image size:               -d\n"
         << "  custom lobe size:                -ms 3\n"
         << "  initial sampling step:           -ss 2\n"
         << "  number of octaves:               -oc 4\n"
         << "  U-SURF (not rotation invariant): -u\n"
         << "  extended descriptor (SURF-128):  -e\n"
         << "  descriptor size:                 -in 4\n"
         << "  input regions:                   -p1 <file>\n"
         << "  verbose output:                  -v\n"
         << "  don't write laplacian:           -nl\n"
         << "  quiet mode:                      -q\n";
    return(0);
  }

  // Read the arguments
  ImLoad ImageLoader;
  int arg = 0;
  string fn = "out.surf";
  Image *im=NULL;
  while (++arg < argc) { 
    if (! strcmp(argv[arg], "-i"))
      im = ImageLoader.readImage(argv[++arg]); 
    if (! strcmp(argv[arg], "-o"))
      fn = argv[++arg];
    if (! strcmp(argv[arg], "-thres"))
      thres = (atof(argv[++arg]))/10000;
    if (! strcmp(argv[arg], "-d"))
      doubleImageSize = true;
    if (! strcmp(argv[arg], "-ms"))
      initLobe = atoi(argv[++arg]);
    if (! strcmp(argv[arg], "-oc"))
      octaves = atoi(argv[++arg]);
    if (! strcmp(argv[arg], "-ss"))
      samplingStep = atoi(argv[++arg]);
    if (! strcmp(argv[arg], "-u"))
      upright = true;
    if (! strcmp(argv[arg], "-e"))
      extended = true;
    if (! strcmp(argv[arg], "-in"))
      indexSize = atoi(argv[++arg]);
    if (! strcmp(argv[arg], "-p1")) {
      bLoadRegions = true;
      sRegionFile  = argv[++arg];
    }
    if (! strcmp(argv[arg], "-v"))
      bVerbose = true;
    if (! strcmp(argv[arg], "-nl"))
      bLaplacian = false;
    if (! strcmp(argv[arg], "-q"))
      bVerbose = false;
  }

  // Start measuring the time
  osmapping::os_GetTime(&tim1);

  // Create the integral image
  Image iimage(im, doubleImageSize);

  // Start finding the SURF points
  if( bVerbose )
    cout << "Finding SURFs...\n";

  // These are the interest points
  vector< Ipoint > ipts;
  ipts.reserve(1000);

  // Extract interest points with Fast-Hessian
  FastHessian fh(&iimage, /* pointer to integral image */
                 ipts,
                 thres, /* blob response threshold */
                 doubleImageSize, /* double image size flag */
                 initLobe * 3 /* 3 times lobe size equals the mask size */, 
                 samplingStep, /* subsample the blob response map */
                 octaves /* number of octaves to be analysed */);

  if( bLoadRegions ) {
    // Load the interest points from disk
    loadIpoints( sRegionFile, ipts, bVerbose );
  } else {
    // Extract them and get their pointer
    fh.getInterestPoints();
  }

  // Initialise the SURF descriptor
  Surf des(&iimage, /* pointer to integral image */  
           doubleImageSize, /* double image size flag */ 
           upright, /* rotation invariance or upright */
           extended, /* use the extended descriptor */
           indexSize /* square size of the descriptor window (default 4x4)*/);

  // Get the length of the descriptor vector resulting from the parameters
  VLength = des.getVectLength();

  // Compute the orientation and the descriptor for every interest point
  for (unsigned n=0; n<ipts.size(); n++){
    //for (Ipoint *k = ipts; k != NULL; k = k->next){
    // set the current interest point
    des.setIpoint(&ipts[n]);
    // assign reproducible orientation
    des.assignOrientation();
    // make the SURF descriptor
    des.makeDescriptor();
  }

  // stop measuring the time, we're all done
  osmapping::os_GetTime(&tim2);

  // save the interest points in the output file
  saveIpoints(fn, ipts, bVerbose, bLaplacian);

  // print some nice information on the command prompt
  if( bVerbose )
    cout << "Detection time: " << osmapping::os_TimeDiff(&tim2, &tim1) << " ms" << endl;

  delete im;

  return 0;
}

// Save the interest points to a regular ASCII file
void saveIpoints(string sFileName, const vector< Ipoint >& ipts, bool bVerbose, bool bLaplacian) 
{
  ofstream ipfile(sFileName.c_str());
  if( !ipfile ) {
    cerr << "ERROR in loadIpoints(): "
         << "Couldn't open file '" << sFileName.c_str() << "'!" << endl; //STS
    return;
  }

  double sc;
  unsigned count = ipts.size();

  // Write the file header
  if (bLaplacian)
    ipfile << VLength + 1 << endl << count << endl;
  else
    ipfile << VLength << endl << count << endl;

  // In order to just save the interest points without descriptor, comment 
  // the above and uncomment the following command.
  // ipfile << 1.0 << endl << count << endl;
  // Save interest point with descriptor in the format of Krystian Mikolajczyk
  // for reasons of comparison with other descriptors. As our interest points 
  // are circular in any case, we use the second component of the ellipse to 
  // provide some information about the strength of the interest point. This is 
  // important for 3D reconstruction as only the strongest interest points are 
  // considered. Replace the strength with 0.0 in order to perform Krystian's 
  // comparisons.
  for (unsigned n=0; n<ipts.size(); n++){
    // circular regions with diameter 5 x scale
    sc = 2.5 * ipts[n].scale; sc*=sc;
    ipfile  << ipts[n].x /* x-location of the interest point */
            << " " << ipts[n].y /* y-location of the interest point */
            << " " << 1.0/sc /* 1/r^2 */
            << " " << 0.0     //(*ipts)[n]->strength /* 0.0 */
            << " " << 1.0/sc; /* 1/r^2 */

    if (bLaplacian)
      ipfile << " " << ipts[n].laplace;

    // Here comes the descriptor
    for (int i = 0; i < VLength; i++) {
      ipfile << " " << ipts[n].ivec[i];
    }
    ipfile << endl;
  }

  // Write message to terminal.
  if( bVerbose )
    cout << count << " interest points found" << endl;
}


// Load the interest points from a regular ASCII file
void loadIpoints(string sFileName, vector< Ipoint >& ipts, bool bVerbose) 
{
  ifstream ipfile(sFileName.c_str());
  if( !ipfile ) {
    cerr << "ERROR in loadIpoints(): "
         << "Couldn't open file '" << sFileName.c_str() << "'!" << endl; //STS
    return;
  }

  // Load the file header
  float    dummy;
  unsigned count;
  ipfile >> dummy >> count;

  // create a new interest point vector
  ipts.clear();
  ipts.resize(count);

  // Load the interest points in Mikolajczyk's format
  for (unsigned n=0; n<count; n++){
    // circular regions with diameter 5 x scale
    float x, y, a, b, c;
    ipfile >> x >> y >> a >> b >> c;

    float det = sqrt((a-c)*(a-c) + 4.0*b*b);
    float e1 = 0.5*(a+c + det);
    float e2 = 0.5*(a+c - det);
    float l1 = (1.0/sqrt(e1));
    float l2 = (1.0/sqrt(e2));
    float sc = sqrt( l1*l2 );

    ipts[n].x     = x;
    ipts[n].y     = y;
    ipts[n].scale = sc/2.5;
  }

  // close the interest point file again
  ipfile.close();

  // Write message to terminal.
  if( bVerbose )
    cout << "read in " << count << " interest points." << endl;
}
