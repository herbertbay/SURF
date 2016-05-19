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

#ifndef __SURF_H
#define __SURF_H

namespace surf {

class Image;

class Surf {
  public:
    //! Constructor
    Surf();

    //! Constructor with parameters
    Surf(Image *im, bool dbl=false, bool usurf=false, 
         bool ext=false, int insi=4);

    //! Destructor
    ~Surf();

    //! Get length of the descriptor vector
    int getVectLength();

    //! set Ipoint for which a descriptor has to be computed
    void setIpoint(Ipoint *ipt);

    //! Assign reproducible orienation
    void assignOrientation();

    //! Compute the robust features
    void makeDescriptor();

  protected:
    //! Create the vector 
    void createVector(double scale, 
                      double row, double col);

    //! Create the vector 
    void createUprightVector(double scale,
                             double row, double col);

    //! Add sample to the vector
    void AddSample(int r, int c, double rpos, 
                   double cpos, double rx, double cx, int step);

    //! Add upright sample to the vector
    void AddUprightSample(int r, int c, double rpos,
                          double cpos, double rx, double cx, int step);

    //! Place sample to index in vector
    void PlaceInIndex(double mag1, int ori1,
                      double mag2, int ori2, double rx, double cx);

    //! Normalise descriptor vector for illumination invariance for
    //! Lambertian surfaces
    void normalise();

    //! Create Lookup tables 
    void createLookups();

  private:
    Image *_iimage;
    Ipoint *_current;
    double ***_index;
    bool _doubleImage;
    bool _upright;
    bool _extended;
    int _VecLength;
    int _IndexSize;
    double _MagFactor;
    int _OriSize;
    int _width, _height;

    double _sine, _cose;
    double **_Pixels;

    double _lookup1[83], _lookup2[40];
};

}

#endif // SURF_H
