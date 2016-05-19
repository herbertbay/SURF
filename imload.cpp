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

#include <iostream>
#include <fstream>

#include "imload.h"
#include "image.h"

namespace surf {

#define MAX(x,y)  (((x) > (y)) ? (x) : (y))
#define MIN(x,y)  (((x) < (y)) ? (x) : (y))
#define ABS(x)    (((x) > 0) ? (x) : (-(x)))

using namespace std;

void ignoreComments(ifstream& imfile) {
  char c;
  do {
    imfile >> c;
  } while (c == ' ');
  imfile.putback(c);

  imfile >> c;
  while (c == '#') {
    imfile.ignore(256, '\n');
    imfile >> c;
  }
  imfile.putback(c);
}

Image *ImLoad::readImage(const char *fn){
  ifstream imfile(fn, ios::binary);
  if (!imfile.is_open()) {
    cerr << "Sorry, could not open: " << fn << endl;
    exit(0);
  }

  // Reading file header
  char P;
  char num;
  imfile >> P >> num; 
  ignoreComments(imfile);

  // Read image dimensions and extremum value
  int width, height, extr;
  imfile >> width;
  ignoreComments(imfile);
  imfile >> height;
  ignoreComments(imfile);
  imfile >> extr;

  // Check whether the file is OK
  if (P != 'P' || num != '5' ||
      width <= 0 || height <= 0 ||
      extr > 255) {
    cerr << "Input image has to be PGM format" << endl;
    exit(0);
  }

  // Get the image intensities and normalise to 0 - 1.
  imfile.get();
  Image *im = new Image(width, height);
  for (int y = 0; y < height; y++)
    for (int x = 0; x < width; x++)
      im->setPix(x, y, ((double) imfile.get()) / extr);
  return im;
}

void ImLoad::saveImage(const char *fn, Image *im) {
  ofstream imfile(fn, ios::binary);
  if (!imfile.is_open()) {
    cerr << "Sorry, could not open: " << fn << endl;
    exit(0);
  }
  imfile << "P5" << endl;
  imfile << im->getWidth() << " " << im->getHeight() << " 255" << endl;
  for (int y = 0; y < im->getHeight(); y++)
    for (int x = 0; x < im->getWidth(); x++)
      imfile.put((unsigned char)(im->getPix(x, y) * 255));
}

}
