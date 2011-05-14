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
#include "prgkern/_prgconfig.h"

#include "prgkern/_assert.h"
#include "prgkern/_average.h"
#include "prgkern/_blas1.h"
#include "prgkern/_box.h"
#include "prgkern/_debug.h"
#include "prgkern/_dense.h"
#include "prgkern/_for.h"
#include "prgkern/_fstring.h"
#include "prgkern/_graph.h"
#include "prgkern/_index.h"
#include "prgkern/_iterator.h"
#include "prgkern/_math.h"
#include "prgkern/_m3x3dense.h"
#include "prgkern/_mdense.h"
#include "prgkern/_mesh.h"
#include "prgkern/_minimize.h"
#include "prgkern/_minimize_lbfgs.h"
#include "prgkern/_minimize_round.h"
#include "prgkern/_os.h"
#include "prgkern/_pproc.h"
#include "prgkern/_random.h"
#include "prgkern/_regex.h"
#include "prgkern/_rotator.h"
#include "prgkern/_stencil.h"
#include "prgkern/_string.h"
#include "prgkern/_time.h"
#include "prgkern/_transforms.h"
#include "prgkern/_type.h"
#include "prgkern/_v3dense.h"
#include "prgkern/_vdense.h"
#include "prgkern/_sse.h"

