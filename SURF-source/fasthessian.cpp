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

#include <iostream>
#include <math.h>
#include <assert.h>

#include "fasthessian.h"
#include "image.h"
#include "ipoint.h"

namespace surf {

using namespace std;

void SolveLinearSystem(double *solution, double sq[3][3], int size);
double DotProd(double *v1, double *v2, int len);
//--------------------- Least-squares solutions ---------------------------
// Give a least-squares solution to the system of linear equations given in
// the jacobian and errvec arrays.  Return result in solution.
// This uses the method of solving the normal equations.
/*void SolveLeastSquares(double *solution, int _height, int _width,
                       double **jacobian, double *errvec, double **sqarray)
{
    int r, c, i;
    double sum;

    assert(_height >= _width);
    // Multiply Jacobian transpose by Jacobian, and put result in sqarray
    for (r = 0; r < _width; r++)
  for (c = 0; c < _width; c++) {
      sum = 0.0;
      for (i = 0; i < _height; i++)
    sum += jacobian[i][r] * jacobian[i][c];
      sqarray[r][c] = sum;
  }
    // Multiply transpose of Jacobian by errvec, and put result in solution
    for (c = 0; c < _width; c++) {
  sum = 0.0;
  for (i = 0; i < _height; i++)
      sum += jacobian[i][c] * errvec[i];
  solution[c] = sum;
    }
    // Now, solve square system of equations
    SolveLinearSystem(solution, sqarray, _width);
}*/


// Solve the square system of linear equations, Ax=b, where A is given
// in matrix "sq" and b in the vector "solution".  Result is given in
// solution.  Uses Gaussian elimination with pivoting
void SolveLinearSystem(double *solution, double sq[3][3], int size)
{
    int row, col, c, pivot = 0, i;
    double maxc, coef, temp, mult, val;

    // Triangularize the matrix
    for (col = 0; col < size - 1; col++) {
  // Pivot row with largest coefficient to top
  maxc = -1.0;
  for (row = col; row < size; row++) {
      coef = sq[row][col];
      coef = (coef < 0.0 ? - coef : coef);
      if (coef > maxc) {
    maxc = coef;
    pivot = row;
      }
  }
  if (pivot != col) {
      // Exchange "pivot" with "col" row (this is no less efficient
            // than having to perform all array accesses indirectly)
      for (i = 0; i < size; i++) {
    temp = sq[pivot][i];
    sq[pivot][i] = sq[col][i];
    sq[col][i] = temp;
      }
      temp = solution[pivot];
      solution[pivot] = solution[col];
      solution[col] = temp;
  }
  // Do reduction for this column
  for (row = col + 1; row < size; row++) {
      mult = sq[row][col] / sq[col][col];
      for (c = col; c < size; c++) //Could start with c=col+1
    sq[row][c] -= mult * sq[col][c];
      solution[row] -= mult * solution[col];
  }
    }

    // Do back substitution.  Pivoting does not affect solution order
    for (row = size - 1; row >= 0; row--) {
  val = solution[row];
  for (col = size - 1; col > row; col--)
      val -= solution[col] * sq[row][col];
  solution[row] = val / sq[row][row];
    }
}


// Return dot product of two vectors with given length
double DotProd(double *v1, double *v2, int len)
{
    int i;
    double sum = 0.0;

    for (i = 0; i < len; i++)
      sum += v1[i] * v2[i];
    return sum;
}

// Constructor
/*FastHessian::FastHessian() {
  // assign variables
  _Iimage = NULL;
  threshold = 0.2;
  _doubled = false;
  _initLobe = 3;
  _maxScales = _initLobe + 2;
//      _maxScales = _initLobe + 1;
  _maxOctaves = 4;
  //ipoints = NULL;
  _sampling = 2 * (1 + _doubled);
  _width = 0;
  _height = 0;
  // allocate variables for interpolation
  allocateVars();
}*/

// Destructor
FastHessian::~FastHessian() {
  if (_scaleLevel) {
    for( int n=0; n<_maxScales; n++ )
      delete _scaleLevel[n];
    delete[] _scaleLevel;
  }
}


// Constructor with parameters
FastHessian::FastHessian(Image *im, vector< Ipoint >& ip, double thres, bool doub,
                         short int initMasksize,
                         short int _samplingStep,
                         short int octaves) : _ipts(ip) {
  // assign variables
  _Iimage = im;
  _threshold = thres;
  _doubled = doub;
  _initLobe = initMasksize / 3;
  _maxScales = _initLobe + 2;
  _maxOctaves = octaves;
  _sampling = _samplingStep * (1 + _doubled);
  _width = _Iimage->getWidth();
  _height = _Iimage->getHeight();

  // allocate octave
  allocateOctave();
}


// Pass the integral image
void FastHessian::setIimage( Image *iim )
{
  _Iimage = iim;
  _width  = iim->getWidth();
  _height = iim->getHeight();
  for( int i=0; i<_maxScales; i++ ) {
    _scaleLevel[i]->setWidth((_width-1)/_sampling);
    _scaleLevel[i]->setHeight((_height-1)/_sampling);
  }
}


// Detect the interest Points
void FastHessian::getInterestPoints()
{
  /* Intensity values of the integral image and the determinant */
  int maskSize = _initLobe-2;
  /* Indices */
  int octave = 1;
  int s=0, k, l;

  // border for non-maximum suppression
  int *borders = new int[_maxScales];

  for (int o=0; o < _maxOctaves; o++){
    // distance to border
    int border1;

    if (o > 0) {
      // divide the last and the third last image in order to save the
      // blob response computation for those two levels
      Image *temp0 = _scaleLevel[0];
      Image *temp1 = _scaleLevel[1];
      _scaleLevel[0] = _scaleLevel[_maxScales-3]->HalfImage();
      _scaleLevel[1] = _scaleLevel[_maxScales-1]->HalfImage();
      _scaleLevel[_maxScales-3] = temp0;
      _scaleLevel[_maxScales-1] = temp1;

      // Calculate border and store in vector
      border1 = ((3*(maskSize + 4*octave) )/2)/(_sampling*octave)+1;
      borders[0] = border1; borders[1] = border1;
      s=2;
    }
    else {
      // Save borders for non-maximum suppression
      border1 = ((3*(maskSize + 6*octave) )/2)/(_sampling*octave)+1;
    }

    // Calculate blob response map for all scales in the octave
    for (; s < _maxScales; s++) {
      borders[s] = border1;

      maskSize += 2*octave;
      if (s>2)
        border1 = ((3*(maskSize))/2)/(_sampling*octave)+1;

      _scaleLevel[s]->setWidth(_scaleLevel[0]->getWidth());
      _scaleLevel[s]->setHeight(_scaleLevel[0]->getHeight());

      _vas[2] = (int)(0.5*maskSize); _vas[3]=2*_vas[2]; _vas[4]=_vas[3]+_vas[2];
      double norm = 9.0/(maskSize*maskSize);
      norm *= norm;

      /* Calculate border size */
      const int border2 = _sampling*octave*border1;
      const int detrows = _scaleLevel[s]->getHeight()-border1;
      const int detcols = _scaleLevel[s]->getWidth()-border1;
      const int delt = _sampling*octave;

      _vas[7] = border2 + maskSize; _vas[8] = border2 - maskSize;
      for ( k=border1, _vas[1] = border2; k < detrows; _vas[1]+=delt, k++){
        _vas[5] = border2 + maskSize; _vas[6] = border2 - maskSize;
        for ( l=border1, _vas[0] = border2; l < detcols; _vas[0]+=delt, l++){
          _scaleLevel[s]->setPix(l, k, norm*_Iimage->getHessian(_vas));
          _vas[5]+=delt; _vas[6]+=delt;
        }
        _vas[7]+=delt; _vas[8]+=delt;
      }
    }
    findMaximum(borders, o, octave);
    octave *= 2;
  }
  delete [] borders;
}

// Create a new ipoint and return list of ipoints with new one added
void FastHessian::makeIpoint(double x, double y, double scale,
                                double strength) {
  double divisor = 1.0;
  if (_doubled)
    divisor = 2.0;

  Ipoint ip;
  ip.x = x*_sampling/divisor;
  ip.y = y*_sampling/divisor;
  ip.scale = 1.2*scale/divisor;
  ip.strength = strength;
  ip.ori = 0.0;
  _vas[0] = int(x*_sampling + 0.5);
  _vas[1] = int(y*_sampling + 0.5);
  _vas[2] = int(3*scale+0.5)/2;
  _vas[3] = 2 * _vas[2];
  _vas[4] = _vas[3] + _vas[2];
  _vas[5] = _vas[0] + int(3*scale+0.5);
  _vas[6] = _vas[0] - int(3*scale+0.5);
  _vas[7] = _vas[1] + int(3*scale+0.5);
  _vas[8] = _vas[1] - int(3*scale+0.5);
  ip.laplace = _Iimage->getTrace(_vas);

  _ipts.push_back(ip);
}

void FastHessian::allocateOctave() {
  _scaleLevel = new Image*[_maxScales];
  for( int n = 0; n < _maxScales; n++)
    _scaleLevel[n] = new Image((_width-1)/_sampling, (_height-1)/_sampling);
}

void FastHessian::findMaximum( int *borders, int o, int octave ) {
  Image *map = NULL;

  const int rows = _scaleLevel[0]->getHeight();
  const int cols = _scaleLevel[0]->getWidth();
  const double thres = 0.8*_threshold;
  double best;
        int r, c, s, ss;
        int dr, dc, ds;
  //map = new Image(cols, rows);

  int cas;
  for( int k = 1; k < _maxScales-1; k+=2 )
    for( int i = borders[k+1]+1; i < rows-(borders[k+1]+1); i+=2 )
      for( int j = borders[k+1]+1; j < cols-(borders[k+1]+1); j+=2 ) {
        best = _scaleLevel[k]->getPix(j, i);
        cas = 0;
        if( _scaleLevel[k]->getPix(j+1, i) > best ) {
          best = _scaleLevel[k]->getPix(j+1, i); cas = 1;}
        if( _scaleLevel[k]->getPix(j, i+1) > best ) {
          best = _scaleLevel[k]->getPix(j, i+1); cas = 2;}
        if( _scaleLevel[k]->getPix(j+1, i+1) > best ) {
          best = _scaleLevel[k]->getPix(j+1, i+1); cas = 3;}

        if( _scaleLevel[k+1]->getPix(j, i) > best ) {
          best = _scaleLevel[k+1]->getPix(j, i); cas = 4;}
        if( _scaleLevel[k+1]->getPix(j+1, i) > best ) {
          best = _scaleLevel[k+1]->getPix(j+1, i); cas = 5;}
        if( _scaleLevel[k+1]->getPix(j, i+1) > best ) {
          best = _scaleLevel[k+1]->getPix(j, i+1); cas = 6;}
        if( _scaleLevel[k+1]->getPix(j+1, i+1) > best ) {
          best = _scaleLevel[k+1]->getPix(j+1, i+1); cas = 7;}
        if (best < thres) continue;
        if (k+1 == _maxScales-1 && cas>3) continue;


        c = j;
        r = i;
        s = k;
        dc = -1;
        dr = -1;
        ds = -1;
        if(cas != 0) {
          if(cas == 1){ c = j+1; dc = 1;} else
          if(cas == 2){ r = i+1; dr = 1;} else
          if(cas == 3){ c = j+1; r=i+1; dc = 1; dr = 1;} else {
            s++;
            ds = 1;
            if(cas == 5){ c = j+1; dc = 1;} else
            if(cas == 6){ r = i+1; dr = 1;} else
            if(cas == 7){ c = j+1; r=i+1; dc = 1; dr = 1;}
          }
        }

        ss = s+ds;
        if( best < _scaleLevel[ss]->getPix(c-1, r-dr)) continue;
        if( best < _scaleLevel[ss]->getPix(c, r-dr)) continue;
        if( best < _scaleLevel[ss]->getPix(c+1, r-dr)) continue;
        if( best < _scaleLevel[ss]->getPix(c-1, r)) continue;
        if( best < _scaleLevel[ss]->getPix(c, r)) continue;
        if( best < _scaleLevel[ss]->getPix(c+1, r)) continue;
        if( best < _scaleLevel[ss]->getPix(c-1, r+dr)) continue;
        if( best < _scaleLevel[ss]->getPix(c, r+dr)) continue;
        if( best < _scaleLevel[ss]->getPix(c+1, r+dr)) continue;
        if( best < _scaleLevel[s]->getPix(c-1, r+dr)) continue;
        if( best < _scaleLevel[s]->getPix(c, r+dr)) continue;
        if( best < _scaleLevel[s]->getPix(c+1, r+dr)) continue;
        if( best < _scaleLevel[s]->getPix(c+dc, r)) continue;
        if( best < _scaleLevel[s]->getPix(c+dc, r-dr)) continue;
        ss = s-ds;
        if( best < _scaleLevel[ss]->getPix(c-1, r+dr)) continue;
        if( best < _scaleLevel[ss]->getPix(c, r+dr)) continue;
        if( best < _scaleLevel[ss]->getPix(c+1, r+dr)) continue;
        if( best < _scaleLevel[ss]->getPix(c+dc, r)) continue;
        if( best < _scaleLevel[ss]->getPix(c+dc, r-dr)) continue;

        interpFeature(s, r, c, map, o, octave, 5, borders);
      }
  //delete map;
}

void FastHessian::interpFeature( int s, int row, int col, Image *map,
                                 int o, int octave, int movesRemain,
                                 int *borders) {
  int newr = row, newc = col;
  double newScale;

  // Interpolate the detected maximum in order to
  // get a more accurate location
  const double strength = fitQuadrat(s, row, col);

  if (_offset[1] > 0.6 && row < _scaleLevel[0]->getHeight() - borders[s])
    newr ++;
  if (_offset[1] < -0.6 && row > borders[s])
    newr --;
  if (_offset[2] > 0.6 && col < _scaleLevel[0]->getWidth() - borders[s])
    newc ++;
  if (_offset[2] < -0.6 && col > borders[s])
    newc --;

  if (movesRemain > 0  &&  (newr != row || newc != col))
    return interpFeature(s, newr, newc, map, o, octave,
                         movesRemain - 1, borders);
  // Do not create a keypoint if interpolation still remains far
  // outside expected limits, or if magnitude of peak value is below
  // threshold (i.e., contrast is too low)
  if(isnan(_offset[0]) || isnan(_offset[1]) || isnan(_offset[2]))
    return;
  if (fabs(_offset[0]) > 1.5  || fabs(_offset[1]) > 1.5  ||
      fabs(_offset[2]) > 1.5 || strength <  _threshold)
    return;

  // Check that no ipoint has been created at this location (to avoid
  // duplicates). Otherwise, mark this map location
  /*if (map->getPix(col, row) > 0.0)
    return;
    map->setPix(col, row, 1.0);
  */

  newScale = (_initLobe+(octave-1)*(_maxScales)+(s+_offset[0])*2*octave)/3.0;

  makeIpoint( octave * (col + _offset[2]), octave * (row + _offset[1]),
              newScale, strength);
}

double FastHessian::fitQuadrat(int s, int r, int c)
{
    // variables to fit quadratic
    double _g[3], _H[3][3];

    // Fill in the values of the gradient from pixel differences
    _g[0] = (_scaleLevel[s+1]->getPix(c, r) -
            _scaleLevel[s-1]->getPix(c, r) ) / 2.0;
    _g[1] = (_scaleLevel[s]->getPix(c, r+1) -
            _scaleLevel[s]->getPix(c, r-1) ) / 2.0;
    _g[2] = (_scaleLevel[s]->getPix(c+1, r) -
            _scaleLevel[s]->getPix(c-1, r) ) / 2.0;

    // Fill in the values of the Hessian from pixel differences
    _H[0][0] = _scaleLevel[s-1]->getPix(c, r) -
        2.0 * _scaleLevel[s]->getPix(c, r) +
              _scaleLevel[s+1]->getPix(c, r);
    _H[1][1] = _scaleLevel[s]->getPix(c, r-1) -
        2.0 * _scaleLevel[s]->getPix(c, r) +
              _scaleLevel[s]->getPix(c, r+1);
    _H[2][2] = _scaleLevel[s]->getPix(c-1, r) -
        2.0 * _scaleLevel[s]->getPix(c, r) +
              _scaleLevel[s]->getPix(c+1, r);
    _H[0][1] = _H[1][0] = ((_scaleLevel[s+1]->getPix(c, r+1) -
                          _scaleLevel[s+1]->getPix(c, r-1)) -
       (_scaleLevel[s-1]->getPix(c, r+1) -
                          _scaleLevel[s-1]->getPix(c, r-1))) / 4.0;
    _H[0][2] = _H[2][0] = ((_scaleLevel[s+1]->getPix(c+1, r) -
                          _scaleLevel[s+1]->getPix(c-1, r)) -
       (_scaleLevel[s-1]->getPix(c+1, r) -
                          _scaleLevel[s-1]->getPix(c-1, r))) / 4.0;
    _H[1][2] = _H[2][1] = ((_scaleLevel[s]->getPix(c+1, r+1) -
                          _scaleLevel[s]->getPix(c-1, r+1)) -
       (_scaleLevel[s]->getPix(c+1, r-1) -
                          _scaleLevel[s]->getPix(c-1, r-1))) / 4.0;

    // Solve the 3x3 linear sytem, Hx = -g.  Result gives peak _offset.
    // Note that SolveLinearSystem destroys contents of H
    _offset[0] = - _g[0];
    _offset[1] = - _g[1];
    _offset[2] = - _g[2];

    SolveLinearSystem(_offset, _H, 3);

    // Also return value of the determinant at peak location using initial
    // value plus 0.5 times linear interpolation with gradient to peak
    // position (this is correct for a quadratic approximation).
    return (_scaleLevel[s]->getPix(c, r) + 0.5 * DotProd(_offset, _g, 3));
}

}
