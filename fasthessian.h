/*
 * Speeded-Up Robust Features (SURF)
 * https://github.com/herbertbay/SURF
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

#ifndef __FASTHESSIAN_H
#define __FASTHESSIAN_H

#include <vector>

namespace surf {

class Ipoint;
class Image;

class FastHessian {
  public:
    //! Destructor
    ~FastHessian();

    //! Constructor with parameters
    FastHessian(Image *im, std::vector< Ipoint >& ip, double thres = 0.2, bool doub = false, 
                short int initMasksize = 9, short int samplingStep = 2,
                short int octaves = 4);

    //! Pass the integral image
    void setIimage( Image *iim );

    //! Detect the interest Points, write into ipts
    void getInterestPoints();

    //! Create a new ipoint at location (x, y),  at a certain scale 
    //! and corner response strength
    void makeIpoint(double x, double y, double scale, double strength=0);

  protected:
    //! Allocate scale layers for one octave
    void allocateOctave();

    //! Fast non-maximum-suppression
    void findMaximum(int *borders, int o, int octave);
    void interpFeature(int s, int row, int col, Image *map,
                       int o, int octave, int movesRemain,
                       int *borders);
    int fitQuadrat(int s, int r, int c, double &res);

  private:
    //! Integral image
    Image *_Iimage;

    //! Octaves
    Image **_scaleLevel;

    //! Vector of variables
    int _vas[9];

    //! Threshold for interest point detection
    double _threshold;

    //! Indicates whether the image size was doubled or not
    //! default is false
    bool _doubled;

    //! Reference to vector of interest points passed from outside
    std::vector< Ipoint >& _ipts;

    //! Initial lobe size for the second derivative in one direction
    //! default is 3
    short int _initLobe;

    //! Number scales
    short int _maxScales;

    //! Number octaves
    short int _maxOctaves;

    //! The sampling step
    short int _sampling;

    //! Integral image dimensions
    int _width;
    int _height;

    //! Result of fitting quadratic
    double _offset[3];
};

}

#endif // FASTHESSIAN_H
