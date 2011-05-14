/***************************************************************************
 *   Copyright (C) 2008 by Eduard S. Fomin                                 *
 *   fomin@bionet.nsc.ru                                                   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef _PRGCONFIG___0077A726_E6E3_58a1_C16D_CE436AB90500__H
#define _PRGCONFIG___0077A726_E6E3_58a1_C16D_CE436AB90500__H

/// I want <math.h> constants definitions
#define _USE_MATH_DEFINES

/// использование точных часов (с разрешением до 1 нс)
#define USE_POSIX_2001

#ifdef __GNUC__

	#define USE_LINUX_SPECIFIC_CODE
	#define _DEBUG

	#ifdef __OPTIMIZE__
		#undef _DEBUG
	#endif

	#define INLINE inline __attribute__((__always_inline__))
	#include <stdint.h>
#else
	#define INLINE inline
#endif

#if defined (__WIN32__) || defined (WIN32)

	#define USE_WINDOWS_SPECIFIC_CODE
	#pragma warning (disable : 4996)
		// warning C4996: 'sprintf' was declared deprecated

#endif

// delete bad standart macro
#undef min
#undef max

#define COUT_MUTEX_LOCK
#define COUT_MUTEX_UNLOCK
	// don't lock output by default

#define SYSTEM_MUTEX_LOCK
#define SYSTEM_MUTEX_UNLOCK
	// don't lock output by default

#define CYCLE_EXPANSION_SIZE  4

#endif
