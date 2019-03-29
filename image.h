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

#ifndef __IMAGE_H
#define __IMAGE_H

namespace surf {

class Image {

  public:
    //! Constructor
    Image(const int w, const int h);

    //! Destructor
    ~Image();

    //! Constructor from existing double array
    Image(double **pixels, int w, int h);

    //! constructor for integral image
    Image(Image *im, bool doubleImSize=false);

    //! Pass a single frame to the (pre-initialized) structure
    void setFrame(unsigned char *im);
    void setFrame(Image *im);

    //! Divide the image size by two
    Image *HalfImage();

    //! Get Hessian response at a certain point
    double getHessian(int *x);

    //! Get Trace of the Hessian
    int getTrace(int *x);

    //! Get the pointer to the image pixels
    double **getPixels() const;

    //! Get the pixel intensity at location (\a x, \a y)
    double getPix(const int x, const int y) const {
      return _pixels[y][x];
    }

    //! Overload of getPix returning the reference
    double &getPix(const int x, const int y) {
      return _pixels[y][x];
    }

    //! Set the Pixel at location (\a x, \a y) to the value "\a val"
    void setPix(const int x, const int y, const double val) {
      _pixels[y][x] = val;
    }

    //! get width
    int getWidth();

    //! get height
    int getHeight();

    //! set width
    void setWidth(int wi);

    //! set height
    void setHeight(int hi);

  protected:
    //! Allocate 2D array of image pixels
    void allocPixels(int w, int h);

  private:
    //! Actual image buffer
    double *_buf;

    //! 2D array of image pixels
    double **_pixels;

    //! Image height and width
    int _height, _width;

    //! Original image height
    int _orihi;

    //! Flag if this image is just a reference
    bool _ref;
};

}

#endif //IMAGE_H
