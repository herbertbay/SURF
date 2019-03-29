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

#ifdef WIN32
#include <windows.h>
#include <process.h>
#else
#include <pthread.h>
#endif

#include "os_mapping.h"

namespace osmapping {

	void os_GetTime(os_TIME *time) { 
	#ifdef WIN32	
		*time = GetTickCount();
	#else
		struct timezone tz;
		gettimeofday(time, &tz);
	#endif
	}

	int os_TimeDiff(os_TIME *time1, os_TIME *time2) { 
	#ifdef WIN32	
		return *time1 - *time2;
	#else
		return (int)((double)time1->tv_sec*1000 + ((double)time1->tv_usec)*1e-3 -
			(double)time2->tv_sec*1000 - ((double)time2->tv_usec)*1e-3);
	#endif
	}
}

