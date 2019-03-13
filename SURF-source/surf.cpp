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

#include <math.h>
#include <assert.h>
#include <vector>
#include <algorithm>

#include "ipoint.h"
#include "surf.h"
#include "image.h"

namespace surf {

#define OriHistTh 0.8
#define window M_PI/3
#define IndexSigma 1.0

#define get_sum(I, x1, y1, x2, y2) (I[y1+1][x1+1] + I[y2][x2] - I[y2][x1+1] - I[y1+1][x2])
#define get_wavelet1(IPatch, x, y, size) (get_sum(IPatch, x + size, y, x - size, y - size) - get_sum(IPatch, x + size, y + size, x - size, y))
#define get_wavelet2(IPatch, x, y, size) (get_sum(IPatch, x + size, y + size, x, y - size) - get_sum(IPatch, x, y + size, x - size, y - size))

#define MAX(x,y)  (((x) > (y)) ? (x) : (y))

using namespace std;

// Constructor
Surf::Surf(){
  _iimage = NULL;
  _Pixels = NULL;
  _doubleImage = false;
  _extended = false;
  _IndexSize = 4;
  _MagFactor = 3;
  _OriSize = 4;
  // calculate length of the descriptor vector
  _VecLength = _IndexSize * _IndexSize * _OriSize;
  // rotation invariance
  _upright = false;
  // allocate _index
  _index = new double**[_IndexSize];
  // Initialize _index array
  for (int i = 0; i < _IndexSize; i++) {
    _index[i] = new double*[_IndexSize];
    for (int j = 0; j < _IndexSize; j++) {
      _index[i][j] = new double[_OriSize];
      for (int k = 0; k < _OriSize; k++)
        _index[i][j][k] = 0.0;
    }
  }
  // create _lookup tables
  createLookups();
  _sine = 0.0;
  _cose = 1.0;
}

// Constructor with parameters
Surf::Surf(Image *im, bool dbl, bool usurf,
           bool ext, int insi){
  // set image
  _iimage = im;
  _Pixels = _iimage->getPixels();
  _doubleImage = dbl;

  _IndexSize = insi;
  _MagFactor = 12/insi;
  _extended = ext;
  _OriSize = 4 + _extended*4;
  // calculate length of the descriptor vector
  _VecLength = _IndexSize * _IndexSize * _OriSize;
  _upright = usurf;
  _width = _iimage->getWidth();
  _height = _iimage->getHeight();

  // create _lookup tables
  createLookups();

  // allocate _index
  _index = new double**[_IndexSize];
  // Initialize _index array
  for (int i = 0; i < _IndexSize; i++) {
    _index[i] = new double*[_IndexSize];
    for (int j = 0; j < _IndexSize; j++) {
      _index[i][j] = new double[_OriSize];
      for (int k = 0; k < _OriSize; k++)
        _index[i][j][k] = 0.0;
    }
  }

  // initial sine and cosine
  _sine = 0.0;
  _cose = 1.0;
}

// Destructor
Surf::~Surf() {
  for (int i = 0; i < _IndexSize; i++) {
    for (int j = 0; j < _IndexSize; j++) {
      delete[] _index[i][j];
    }
    delete[] _index[i];
  }
  delete[] _index;
}

// Get length of the descriptor vector
int Surf::getVectLength(){
  return _VecLength;
}

// set Ipoint for which a descriptor has to be computed
void Surf::setIpoint(Ipoint* ipt){
  _current = ipt;
}

// Assign orienationt
void Surf::assignOrientation() {
  if (_upright)
    return;

  double scale = (1.0+_doubleImage)* _current->scale;
  int x = (int)((1.0+_doubleImage) * _current->x + 0.5);
  int y = (int)((1.0+_doubleImage) * _current->y + 0.5);
  int pixSi = (int)(2*scale + 1.6);
  const int pixSi_2 = (int)(scale + 0.8);
  double weight;
  const int radius=9;
  double dx=0, dy=0, magnitude, angle, distsq;
  const double radiussq = 81.5;
  int y1, x1;
  int yy, xx;

  vector< pair< double, double > > values;
  for (yy = y - pixSi_2*radius, y1= -radius; y1 <= radius; y1++,
       yy+=pixSi_2){
    for (xx = x - pixSi_2*radius, x1 = -radius; x1 <= radius; x1++,
         xx+=pixSi_2) {
      // Do not use last row or column, which are not valid
      if (yy + pixSi + 2 < _height &&
          xx + pixSi + 2 < _width &&
          yy - pixSi > -1 &&
          xx - pixSi > -1) {
        distsq = (y1 * y1 + x1 * x1);
        if (distsq < radiussq) {
          weight = _lookup1[(int)distsq];
          dx = get_wavelet2(_iimage->getPixels(), xx, yy, pixSi);
          dy = get_wavelet1(_iimage->getPixels(), xx, yy, pixSi);

          magnitude = sqrt(dx * dx + dy * dy);
          if (magnitude > 0.0){
            angle = atan2(dy, dx);
            values.push_back( make_pair( angle, weight*magnitude ) );
          }
        }
      }
    }
  }

  double best_angle = 0;

  if (values.size()) {
    sort( values.begin(), values.end() );
    int N = values.size();

    float d2Pi = 2.0*M_PI;
    for( int i = 0; i < N; i++ ) {
      values.push_back( values[i] );
      values.back().first += d2Pi;
    }

    double part_sum = values[0].second;
    double best_sum = 0;
    double part_angle_sum = values[0].first * values[0].second;

    for( int i = 0, j = 0; i < N && j<2*N; ) {
      if( values[j].first - values[i].first < window ) {
        if( part_sum > best_sum ) {
          best_angle  = part_angle_sum / part_sum;
          best_sum = part_sum;
        }
        j++;
        part_sum += values[j].second;
        part_angle_sum += values[j].second * values[j].first;
      }
      else {
        part_sum -= values[i].second;
        part_angle_sum -= values[i].second * values[i].first;
        i++;
      }
    }
  }
  _current->ori = best_angle;
}

// Compute the robust features
void Surf::makeDescriptor() {
  _current->allocIvec(_VecLength);
  // Initialize _index array
  for (int i = 0; i < _IndexSize; i++) {
    for (int j = 0; j < _IndexSize; j++) {
      for (int k = 0; k < _OriSize; k++)
        _index[i][j][k] = 0.0;
    }
  }

  if (_upright){
    // Produce _upright sample vector
    createUprightVector(1.65*(1+_doubleImage)*_current->scale,
                        (1+_doubleImage)*_current->y,
                        (1+_doubleImage)*_current->x);
  } else {
    // calculate _sine and co_sine once
    _sine = sin(_current->ori);
    _cose = cos(_current->ori);

    // Produce _upright sample vector
    createVector(1.65*(1+_doubleImage)*_current->scale,
                 (1+_doubleImage)*_current->y,
                 (1+_doubleImage)*_current->x);
  }

  int v = 0;
  for (int i = 0; i < _IndexSize; i++){
    for (int j = 0; j < _IndexSize; j++){
      for (int k = 0; k < _OriSize; k++)
        _current->ivec[v++] = _index[i][j][k];
    }
  }
  normalise();
}

// -------------------------------------------------------------------------------
// protected:
// -------------------------------------------------------------------------------

// Create Descriptor vector
void Surf::createVector(double scale, double y, double x) {
  int i, j, iradius, iy, ix;
  double spacing, radius, rpos, cpos, rx, cx;
  int step = MAX((int)(scale/2 + 0.5),1);

  iy = (int) (y + 0.5);
  ix = (int) (x + 0.5);

  double fracy = y-iy;
  double fracx = x-ix;
  double fracr =   _cose * fracy + _sine * fracx;
  double fracc = - _sine * fracy + _cose * fracx;

  // The spacing of _index samples in terms of pixels at this scale
  spacing = scale * _MagFactor;

  // Radius of _index sample region must extend to diagonal corner of
  // _index patch plus half sample for interpolation.
  radius = 1.4 * spacing * (_IndexSize + 1) / 2.0;
  iradius = (int) (radius/step + 0.5);

  // Examine all points from the gradient image that could lie within the
  // _index square.
  for (i = -iradius; i <= iradius; i++)
    for (j = -iradius; j <= iradius; j++) {
      // Rotate sample offset to make it relative to key orientation.
      // Uses (x,y) coords.  Also, make subpixel correction as later image
      // offset must be an integer.  Divide by spacing to put in _index units.
      rpos = (step*(_cose * i + _sine * j) - fracr) / spacing;
      cpos = (step*(- _sine * i + _cose * j) - fracc) / spacing;

      // Compute location of sample in terms of real-valued _index array
      // coordinates.  Subtract 0.5 so that rx of 1.0 means to put full
      // weight on _index[1] (e.g., when rpos is 0 and _IndexSize is 3.
      rx = rpos + _IndexSize / 2.0 - 0.5;
      cx = cpos + _IndexSize / 2.0 - 0.5;

      // Test whether this sample falls within boundary of _index patch
      if (rx > -1.0 && rx < (double) _IndexSize  &&
          cx > -1.0 && cx < (double) _IndexSize)
        AddSample(iy + i*step, ix + j*step, rpos, cpos,
                  rx, cx, (int)(scale));
    }
}

// Create Upright Descriptor vector
void Surf::createUprightVector(double scale, double y, double x) {
  int i, j, iradius, iy, ix;
  double spacing, radius, rpos, cpos, rx, cx, dx, dy;
  int step = MAX((int)(scale/2 + 0.5),1);

  iy = (int) (y + 0.5);
  ix = (int) (x + 0.5);
  dx = x - ix;
  dy = y - iy;

  spacing = scale * _MagFactor;

  radius = spacing * (_IndexSize + 1) / 2.0;
  iradius = (int) (radius/step + 0.5);

  for (i = -iradius; i <= iradius; i++)
    for (j = -iradius; j <= iradius; j++) {
      rpos = (step * i - dy) / spacing;
      cpos = (step * j - dx) / spacing;

      rx = rpos + _IndexSize / 2.0 - 0.5;
      cx = cpos + _IndexSize / 2.0 - 0.5;

      if (rx > -1.0 && rx < (double) _IndexSize  &&
          cx > -1.0 && cx < (double) _IndexSize)
        AddUprightSample(iy + i*step, ix + j*step, rpos, cpos,
                         rx, cx, (int)(scale));
    }
}

// Add Sample in the descriptor vector
void Surf::AddSample(int r, int c, double rpos,
                     double cpos, double rx, double cx, int step) {
  double weight;
  double dx, dy;

  // Clip at image boundaries.
  if (r < 1+step  ||  r >= _height - 1-step  ||
      c < 1+step  ||  c >= _width - 1-step)
     return;

  weight = _lookup2[int(rpos * rpos + cpos * cpos)];
  double dxx, dyy;

  dxx = weight*get_wavelet2(_Pixels, c, r, step);
  dyy = weight*get_wavelet1(_Pixels, c, r, step);
  dx = _cose*dxx + _sine*dyy;
  dy = _sine*dxx - _cose*dyy;

  //64
  if(!_extended)
    PlaceInIndex(dx, (dx<0?0:1), dy, (dy<0?2:3), rx, cx);
  else {
    //128
    PlaceInIndex(dx, (dy<0?0:1), fabs(dx), (dy<0?2:3), rx, cx);
    PlaceInIndex(dy, (dx<0?4:5), fabs(dy), (dx<0?6:7), rx, cx);
  }
}

// Add Sample in the _upright descriptor vector
void Surf::AddUprightSample(int r, int c, double rpos,
                            double cpos, double rx, double cx, int step) {
  double weight;
  double dx, dy;

  // Clip at image boundaries.
  if (r < 1+step || r >= _height - 1-step  ||
      c < 1+step || c >= _width - 1-step)
     return;

  weight = _lookup2[int(rpos * rpos + cpos * cpos)];
  dx = weight*get_wavelet2(_Pixels, c, r, step);
  dy = weight*get_wavelet1(_Pixels, c, r, step);

  //64
  if(!_extended)
    PlaceInIndex(dx, (dx<0?0:1), dy, (dy<0?2:3), rx, cx);
  else {
    //128
    PlaceInIndex(dx, (dy<0?0:1), fabs(dx), (dy<0?2:3), rx, cx);
    PlaceInIndex(dy, (dx<0?4:5), fabs(dy), (dx<0?6:7), rx, cx);
  }
}

// Place in _index
void Surf::PlaceInIndex(double mag1, int ori1,
                        double mag2, int ori2, double rx, double cx) {
  int ri, ci, r_index, c_index;
  double rfrac, cfrac, cfrac1, rweight1, cweight1, rweight2, cweight2;

//  fprintf(stderr,"%f %f\n", rx, cx);
  ri = (int)((rx >= 0.0) ? rx : rx - 1.0);
  ci = (int)((cx >= 0.0) ? cx : cx - 1.0);
  rfrac = rx - ri;
  cfrac = cx - ci;
  cfrac1 = 1 - cfrac;
  //assert(ri >= -1  &&  ri < _IndexSize  && rfrac >= 0.0 && rfrac <= 1.0);

  r_index = ri;
  if (r_index >=0) {
    rweight1 = mag1 * (1.0 - rfrac);
    rweight2 = mag2 * (1.0 - rfrac);

    c_index = ci;
    if (c_index >=0) {
      cweight1 = rweight1 * (cfrac1);
      cweight2 = rweight2 * (cfrac1);
      _index[r_index][c_index][ori1] += cweight1;
      _index[r_index][c_index][ori2] += cweight2;
    }
    c_index++;
    if (c_index < _IndexSize) {
      cweight1 = rweight1 * cfrac;
      cweight2 = rweight2 * cfrac;
      _index[r_index][c_index][ori1] += cweight1;
      _index[r_index][c_index][ori2] += cweight2;
    }
  }
  r_index++;
  if (r_index < _IndexSize) {
    rweight1 = mag1 * rfrac;
    rweight2 = mag2 * rfrac;

    c_index = ci;
    if (c_index >=0) {
      cweight1 = rweight1 * (cfrac1);
      cweight2 = rweight2 * (cfrac1);
      _index[r_index][c_index][ori1] += cweight1;
      _index[r_index][c_index][ori2] += cweight2;
    }
    c_index++;
    if (c_index < _IndexSize) {
      cweight1 = rweight1 * cfrac;
      cweight2 = rweight2 * cfrac;
      _index[r_index][c_index][ori1] += cweight1;
      _index[r_index][c_index][ori2] += cweight2;
    }
  }
}


// Normalise descriptor vector for illumination invariance for
// Lambertian surfaces
void Surf::normalise() {
  double val, sqlen = 0.0, fac;
  for (int i = 0; i < _VecLength; i++){
    val = _current->ivec[i];
    sqlen += val * val;
  }
  fac = 1.0/sqrt(sqlen);
  for (int i = 0; i < _VecLength; i++)
    _current->ivec[i] *= fac;
}

// Create _lookup tables
void Surf::createLookups(){
  for (int n=0;n<83;n++)
    _lookup1[n]=exp(-((double)(n+0.5))/12.5);

  for (int n=0;n<40;n++)
    _lookup2[n]=exp(-((double)(n+0.5))/8.0);
}

}
