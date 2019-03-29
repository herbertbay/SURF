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

#ifndef __os_mapping_h
#define __os_mapping_h

#ifdef WIN32
#include <windows.h>
#else
#include <sys/types.h>
#include <sys/time.h>
#endif

namespace osmapping {

#ifdef  __cplusplus
extern "C" {
#endif


	#ifdef WIN32
	typedef DWORD os_TIME;
	#else
	typedef struct timeval os_TIME;
	#endif

	void os_GetTime(os_TIME *time);
	int os_TimeDiff(os_TIME *time1, os_TIME *time2);


#ifdef  __cplusplus
}
#endif
}

#endif

