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
 * For further details, contact Andreas Ess (aess@vision.ee.ethz.ch).
 */

#ifndef IPOINT_H
#define IPOINT_H

#include <stdlib.h>

namespace surf {

class Ipoint {
  public:
  // Constructor
  Ipoint(){
    ivec = NULL;
    ori = 0.0;
  };

  // Destructor
  ~Ipoint(){
    if (ivec)
      delete [] ivec;
  };

  // Allocate space
  void allocIvec(const int si){
    ivec = new double[si];
  };

    // x, y value of the interest point
    double x, y;
    // detected scale
    double scale;
    // strength of the interest point
    double strength;
    // orientation
    double ori;
    // sign of Laplacian
    int laplace;
    // descriptor
    double *ivec;
};

}

#endif // IPOINT_H
