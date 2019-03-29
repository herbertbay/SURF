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

/**
 * SURF library functions
 **/

#ifndef __SURFLIB_H
#define __SURFLIB_H

#include "ipoint.h"
#include "fasthessian.h"
#include "surf.h"
#include "image.h"

namespace surf {

/**
 * Identify interest points and calculate their descriptor
 *
 * @param im pointer to double image
 * @param ipts (return) vector of interest points
 * @param thres blob response threshold
 * @param doubleImageSize double image size
 * @param initLobe custom lobe size
 * @param samplingStep initial sampling step
 * @param octaves number of octaves
 * @param upright true to switch off rotation invariance
 * @param extended true for SURF-128 instead of SURF-64
 * @param indexSize descriptor size
 **/
inline void surfDetDes(Image *im, std::vector< Ipoint >& ipts,
				double thres = 4.0, bool doubleImageSize = false,
				int initLobe = 3, int samplingStep = 2, int octaves = 4,
				bool upright = false, bool extended = false, int indexSize = 4) {
  // Create the integral image
  Image iimage(im, doubleImageSize);

  // Extract interest points with Fast-Hessian
  FastHessian fh(&iimage, /* pointer to integral image */
                 ipts, /* interest point vector to be filled */
                 thres, /* blob response threshold */
                 doubleImageSize, /* double image size flag */
                 initLobe * 3 /* 3 times lobe size equals the mask size */, 
                 samplingStep, /* subsample the blob response map */
                 octaves /* number of octaves to be analysed */);

  // Extract them and get their pointer
  fh.getInterestPoints();

  // Initialise the SURF descriptor
  Surf des(&iimage, /* pointer to integral image */  
           doubleImageSize, /* double image size flag */ 
           upright, /* rotation invariance or upright */
           extended, /* use the extended descriptor */
           indexSize /* square size of the descriptor window (default 4x4)*/);

  // Compute the orientation and the descriptor for every interest point
  for (unsigned n=0; n<ipts.size(); n++){
    // set the current interest point
    des.setIpoint(&ipts[n]);
    // assign reproducible orientation
    des.assignOrientation();
    // make the SURF descriptor
    des.makeDescriptor();
  }
}

/**
 * Calculate descriptor for given interest points
 *
 * @param im pointer to double image
 * @param ipts (return) vector of interest points
 * @param doubleImageSize double image size
 * @param upright true to switch off rotation invariance
 * @param extended true for SURF-128 instead of SURF-64
 * @param indexSize descriptor size
 **/
inline void surfDes(Image *im, std::vector< Ipoint >& ipts,
			 bool doubleImageSize = false,
			 bool upright = false, bool extended = false, int indexSize = 4) {
  // Create the integral image
  Image iimage(im, doubleImageSize);

  // Initialise the SURF descriptor
  Surf des(&iimage, /* pointer to integral image */  
           doubleImageSize, /* double image size flag */ 
           upright, /* rotation invariance or upright */
           extended, /* use the extended descriptor */
           indexSize /* square size of the descriptor window (default 4x4)*/);

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
}

}

#endif
