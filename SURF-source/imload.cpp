/*
 * Speeded-Up Robust Features (SURF)
 * https://github.com/herbertbay/SURF
 *
 * Authors: Herbert Bay, Andreas Ess, Geert Willems
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
 * its documentation for educational, research, and not-for-profit
 * purposes, without fee and without a signed licensing agreement, is
 * hereby granted, provided that the above copyright notice and this
 * paragraph appear in all copies modifications, and distributions.
 *
 * Any commercial use or any redistribution of this software
 * requires a license from one of the above mentioned establishments.
 *
 * For further details, contact Herbert Bay (bay@vision.ee.ethz.ch).
 */

#include <iostream>
#include <fstream>

#include "imload.h"
#include "image.h"

namespace surf {

#define MAX(x,y)  (((x) > (y)) ? (x) : (y))
#define MIN(x,y)  (((x) < (y)) ? (x) : (y))
#define ABS(x)    (((x) > 0) ? (x) : (-(x)))

using namespace std;

Image *ImLoad::readImage(const char *fn){
  ifstream imfile(fn);
  if (!imfile.is_open()){
    cerr << "Sorry, could not open: " << fn << endl;
    exit(0);
  }

  // Reading file header
  char P;
  char num;
  imfile >> P >> num;
  imfile.ignore(256,'\n'); // ignore potential comments

  // Ignore some potential comments
  char c;
  imfile >> c;
  while (c == '#') {
    imfile.ignore(256,'\n');
    imfile >> c;
  }
  imfile.unget();

  // Read image dimensions and extremum value
  int width, height, extr;
  imfile >> width >> height;
  // ignore potential comments
  imfile.ignore(256,'\n');
  imfile >> extr;

  // Check whether the file is OK
  if (P != 'P' || num != '5' ||
      width <= 0 || height <= 0 ||
      extr > 255) {
    cerr << "Input image has to be PGM format" << endl;
    exit(0);
  }

  // Get the image intensities and normalise to 1.
  imfile.get();
  Image *im = new Image(width, height);
  for (int y = 0; y < height; y++)
    for (int x = 0; x < width; x++)
      im->setPix(x, y, ((double) imfile.get()) / extr);
  return im;
}

}
