/*
 * Speeded-Up Robust Features (SURF)
 * http://people.ee.ethz.ch/~surf
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

#include "image.h"
#include <cassert>

namespace surf {

#define get_sum(x1, y1, x2, y2) (_pixels[y1+1][x1+1] + _pixels[y2][x2] - _pixels[y2][x1+1]   - _pixels[y1+1][x2])

// Constructor
Image::Image(const int w, const int h) {
  allocPixels(w, h);
}

// Destructor
Image::~Image() {
  if (!_ref) {
    delete[] _pixels;
    delete[] _buf;
  }
}

// Constructor from existing array
Image::Image(double **im, int w, int h) {
  _pixels = im;
  _width = w;
  _orihi = _height = h;
  _ref = true;
}

// Constructor for integral image
Image::Image(Image *im, bool doubleImSize) {
  if (doubleImSize){
    allocPixels(2 * im->_width - 1, 2 * im->_height - 1);

    const int hi = im->_height;
    const int wi = im->_width;

    // Indices
    int i, j, y2, x2;
    double p1, p2, p3, p4;

    // set first row and column to zero
    for (i=0; i < _height; i++)
      _pixels[i][0] = 0.0;
    for (i=1; i < _width; i++)
      _pixels[0][i] = 0.0;

    _pixels[1][1] = im->_pixels[0][0];
    _pixels[2][1] = _pixels[1][1] + 0.5*(im->_pixels[0][0]+im->_pixels[1][0]);
    _pixels[1][2] = _pixels[1][1] + 0.5*(im->_pixels[0][0]+im->_pixels[0][1]);
    _pixels[2][2] = _pixels[2][1] + _pixels[1][2] +
      0.25 * (im->_pixels[0][0] + im->_pixels[1][0] + im->_pixels[0][1] +
              im->_pixels[1][1]) - _pixels[1][1];

    // Calculate first row
    for (i = 1; i < hi-1; i++){
      y2 = 2*i;
      p1 = im->_pixels[i][0];
      p2 = 0.5*(im->_pixels[i][0]+im->_pixels[i][1]);
      p3 = 0.5*(im->_pixels[i+1][0]+im->_pixels[i][0]);
      p4 = 0.25 * (im->_pixels[i+1][0] + im->_pixels[i][0] +
                   im->_pixels[i][1] + im->_pixels[i+1][1]);
      _pixels[y2+1][1] = _pixels[y2][1] + p1;
      _pixels[y2+1][2] = _pixels[y2][2] + _pixels[y2+1][1]+ p2 - _pixels[y2][1];
      _pixels[y2+2][1] = _pixels[y2+1][1] + p3;
      _pixels[y2+2][2] = _pixels[y2+1][2] + _pixels[y2+2][1] + p4 -
        _pixels[y2+1][1];
    }

    // Calculate the rest of the image
    for (j = 1; j < wi-1; j++) {
      x2 = 2*j;
      p1 = im->_pixels[0][j];
      p2 = 0.5 * (im->_pixels[0][j] + im->_pixels[0][j+1]);
      p3 = 0.5 * (im->_pixels[0][j] + im->_pixels[1][j]);
      p4 = 0.25 * (im->_pixels[1][j] + im->_pixels[0][j] +
                   im->_pixels[0][j+1] + im->_pixels[1][j+1]);
      _pixels[1][x2+1] = _pixels[1][x2] + p1;
      _pixels[1][x2+2] = _pixels[1][x2+1] + p2;
      _pixels[2][x2+1] = _pixels[1][x2+1] + _pixels[2][x2] + p3 -
        _pixels[1][x2];
      _pixels[2][x2+2] = _pixels[1][x2+2] + _pixels[2][x2+1] + p4 -
        _pixels[1][x2+1];
      for (i = 1; i < hi-1; i++) {
        y2 = 2*i;
        p1 = im->_pixels[i][j];
        p2 = 0.5*(im->_pixels[i][j]+im->_pixels[i][j+1]);
        p3 = 0.5*(im->_pixels[i+1][j]+im->_pixels[i][j]);
        p4 = 0.25 * (im->_pixels[i+1][j] + im->_pixels[i][j] +
                     im->_pixels[i][j+1] + im->_pixels[i+1][j+1]);
        _pixels[y2+1][x2+1] = _pixels[y2][x2+1] + _pixels[y2+1][x2] + p1 -
          _pixels[y2][x2];
        _pixels[y2+1][x2+2] = _pixels[y2][x2+2] + _pixels[y2+1][x2+1] + p2 -
          _pixels[y2][x2+1];
        _pixels[y2+2][x2+1] = _pixels[y2+1][x2+1] + _pixels[y2+2][x2] + p3 -
          _pixels[y2+1][x2];
        _pixels[y2+2][x2+2] = _pixels[y2+1][x2+2] + _pixels[y2+2][x2+1] + p4 -
          _pixels[y2+1][x2+1];
      }
    }
  } else {
    allocPixels(im->_width + 1, im->_height + 1);

    // set first row and column to zero
/*    for (int i=0; i < _height; i++)
      _pixels[i][0] = 0.0;
    for (int i=1; i < _width; i++)
      _pixels[0][i] = 0.0;

    double s;
    // Indices
    int i, j;
    // Initialise first image value
    _pixels[1][1] = im->_pixels[0][0];
    // Calculate first row
    for (i = 1; i < _height-1; i++)
      _pixels[i+1][1] = _pixels[i][1] + im->_pixels[i][0];

    // Calculate the rest of the image
    for (j = 1; j < _width-1; j++) {
      _pixels[1][j+1] = _pixels[1][j] + im->_pixels[0][j];
      s = im->_pixels[0][j];
      for (i = 1; i < _height-1; i++) {
        s += im->_pixels[i][j];
        _pixels[i+1][j+1] = s + _pixels[i+1][j];
      }
    }*/

    // set first row and column to zero
    for (int i = 0; i < _height; i++)
      _pixels[i][0] = 0.0;
    for (int i = 1; i < _width; i++)
      _pixels[0][i] = 0.0;

    double s;
    double *b = im->_buf;
    for (int j = 1; j < _height; j++) {
      s = 0;
      for (int i = 1; i < _width; i++) {
        s += *b++;
        _pixels[j][i] = s + _pixels[j - 1][i];
      }
    }
  }
}

// Pass a single frame to the (pre-initialized) structure
void Image::setFrame(unsigned char *im)
{
  int wi = 640;
  // set first row and column to zero
  for( int i=0; i<_height; i++ )
    _pixels[i][0] = 0.0;
  for( int i=1; i<_width; i++ )
    _pixels[0][i] = 0.0;

  double s;
  // Indices
  int i, j;
  // Initialize first image value
  _pixels[1][1] = ((double)(im[0]))/255;
  // Calculate first row
  for( i=1; i<_height-1; i++ )
    _pixels[i+1][1] = _pixels[i][1] + (double (im[i*wi]))/255;

  // Calculate the rest of the image
  for( j=1; j<_width-1; j++ ) {
    _pixels[1][j+1] = _pixels[1][j] + (double (im[j]))/255;
    s = (double(im[j]))/255;
    for( i=1; i<_height-1; i++ ) {
      s += (double (im[i*wi+j]))/255;
      _pixels[i+1][j+1] = s + _pixels[i+1][j];
    }
  }
}

// Pass a single frame to the (pre-initialized) structure
void Image::setFrame(Image *im)
{
  int wi = im->getWidth();
  int hi = im->getHeight();
  assert( wi==_width && hi==_height );

  // set first row and column to zero
  for( int i=0; i<_height; i++ )
    _pixels[i][0] = 0.0;
  for( int i=1; i<_width; i++ )
    _pixels[0][i] = 0.0;

  double s;
  // Indices
  int i, j;
  // Initialize first image value
  _pixels[1][1] = im->_pixels[1][1];
  // Calculate first row
  for( i=1; i<_height-1; i++ )
    _pixels[i+1][1] = _pixels[i][1] + im->_pixels[i+1][1];

  // Calculate the rest of the image
  for( j=1; j<_width-1; j++ ) {
    _pixels[1][j+1] = _pixels[1][j] + im->_pixels[1][j+1];
    s = im->_pixels[1][j+1];
    for( i=1; i<_height-1; i++ ) {
      s += im->_pixels[i+1][j+1];
      _pixels[i+1][j+1] = s + _pixels[i+1][j];
    }
  }
}


// Divide the image size by two
Image *Image::HalfImage() {
  int nwidth = _width / 2;
  int nheight = _height / 2;
  int yi = 0;

  for (int y = 0; y < nheight; y++, yi += 2) {
    int xi = 0;
    for (int x = 0; x < nwidth; x++, xi += 2)
      _pixels[y][x] = _pixels[yi][xi];
  }
  _width = nwidth;
  _height = nheight;
  return this;
}

// Get the Hessian response
double Image::getHessian(int *x){
  double Lxx, Lyy, Lxy;
  /* Get second order derivatives */
  Lxx = (get_sum(x[5] + x[2], x[1] + x[3], x[6] - x[2], x[1] - x[3])-
         3 * get_sum(x[0] + x[2], x[1] + x[3], x[0] - x[2], x[1] - x[3]));
  Lyy = (get_sum(x[0] + x[3], x[7] + x[2], x[0] - x[3], x[8] - x[2])-
         3 * get_sum(x[0] + x[3], x[1] + x[2], x[0] - x[3], x[1] - x[2]));
  Lxy =  0.6*(get_sum(x[0] + x[4], x[1], x[0], x[1] - x[4]) +
              get_sum(x[0], x[1] + x[4], x[0] - x[4], x[1]) -
              get_sum(x[0] + x[4], x[1] + x[4], x[0], x[1]) -
              get_sum(x[0], x[1], x[0] - x[4], x[1] - x[4]) );
  return Lxx*Lyy-Lxy*Lxy;
}

// Get the Hessian trace response
int Image::getTrace(int *x){
  double Lxx, Lyy;
  /* Get second order derivatives */
  Lxx = (get_sum(x[5] + x[2], x[1] + x[3], x[6] - x[2], x[1] - x[3])-
         3 * get_sum(x[0] + x[2], x[1] + x[3], x[0] - x[2], x[1] - x[3]));
  Lyy = (get_sum(x[0] + x[3], x[7] + x[2], x[0] - x[3], x[8] - x[2])-
         3 * get_sum(x[0] + x[3], x[1] + x[2], x[0] - x[3], x[1] - x[2]));
  return (Lxx + Lyy > 0 ? 1 : -1);
}

// Return the pointer to the pixel array
double **Image::getPixels() const {
  return _pixels;
}

// get width
int Image::getWidth() {
  return _width;
}

// get height
int Image::getHeight() {
  return _height;
}

// set width
void Image::setWidth(int wi) {
  _width = wi;
}

// set height
void Image::setHeight(int hi) {
  _height = hi;
}

// protected functions
//

// 2D array of image _pixels
void Image::allocPixels(int w, int h) {
  _width = w;
  _orihi = _height = h;
  _buf = new double[_width * _height];
  _pixels = new double*[_height];
  for (int i = 0; i < _height; i++)
    _pixels[i] = &_buf[i * _width];
  _ref = false;

/*  double **m = new double*[_height];
  for (int i = 0; i < _height; i++)
    m[i] = new double[_width];
  return m;*/
}

}
