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
#ifndef _DEBUG___0077A726_E6E3_58a1_C16D_CE436AB90500__H
#define _DEBUG___0077A726_E6E3_58a1_C16D_CE436AB90500__H

#include <iostream>

#define _MSG(s) std::cout << "[DEBUG] " << s << std::endl;

#define _MSG1(s, s1)           std::cout << "[DEBUG] " << s << " " \
	<< #s1 << "=" << prgkern::make_string(s1) << std::endl;

#define _MSG2(s, s1, s2)       std::cout << "[DEBUG] " << s << " " \
	<< #s1 << "=" << prgkern::make_string(s1) << " " \
	<< #s2 << "=" << prgkern::make_string(s2) << std::endl;

#define _MSG3(s, s1, s2, s3)   std::cout << "[DEBUG] " << s << " " \
	<< #s1 << "=" << prgkern::make_string(s1) << " " \
	<< #s2 << "=" << prgkern::make_string(s2) << " " \
	<< #s3 << "=" << prgkern::make_string(s3) << std::endl;

#define _MSG4(s, s1, s2, s3, s4)  std::cout << "[DEBUG] " << s << " " \
	<< #s1 << "=" << prgkern::make_string(s1) << " " \
	<< #s2 << "=" << prgkern::make_string(s2) << " " \
	<< #s3 << "=" << prgkern::make_string(s3) << " " \
	<< #s4 << "=" << prgkern::make_string(s4) << std::endl;

#define _MSG5(s, s1, s2, s3, s4, s5)  std::cout << "[DEBUG] " << s << " " \
	<< #s1 << "=" << prgkern::make_string(s1) << " " \
	<< #s2 << "=" << prgkern::make_string(s2) << " " \
	<< #s3 << "=" << prgkern::make_string(s3) << " " \
	<< #s4 << "=" << prgkern::make_string(s4) << " " \
	<< #s5 << "=" << prgkern::make_string(s5) << std::endl;

#define _VAL(s)               std::cout << "[DEBUG] " \
	<< #s << "=" << prgkern::make_string(s) << std::endl;

#define _VAL1(s)               std::cout << "[DEBUG] " \
	<< #s << "=" << prgkern::make_string(s) << std::endl;

#define _VAL2(s1, s2)          std::cout << "[DEBUG] " \
	<< #s1 << "=" << prgkern::make_string(s1) << " " \
	<< #s2 << "=" << prgkern::make_string(s2) << std::endl;

#define _VAL3(s1, s2, s3)      std::cout << "[DEBUG] " \
	<< #s1 << "=" << prgkern::make_string(s1) << " " \
	<< #s2 << "=" << prgkern::make_string(s2) << " " \
	<< #s3 << "=" << prgkern::make_string(s3) << std::endl;

#define _VAL4(s1, s2, s3, s4)  std::cout << "[DEBUG] " \
	<< #s1 << "=" << prgkern::make_string(s1) << " " \
	<< #s2 << "=" << prgkern::make_string(s2) << " " \
	<< #s3 << "=" << prgkern::make_string(s3) << " " \
	<< #s4 << "=" << prgkern::make_string(s4) << std::endl;

#define _VAL5(s1, s2, s3, s4, s5)  std::cout << "[DEBUG] " \
	<< #s1 << "=" << prgkern::make_string(s1) << " " \
	<< #s2 << "=" << prgkern::make_string(s2) << " " \
	<< #s3 << "=" << prgkern::make_string(s3) << " " \
	<< #s4 << "=" << prgkern::make_string(s4) << " " \
	<< #s5 << "=" << prgkern::make_string(s5) << std::endl;

#endif
