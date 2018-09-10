/***************************************************************************
 *   Copyright (C) 2005 by Tian-Li Yu,,,                                   *
 *   tianliyu@fishlaptop.ytgroup                                           *
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

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <float.h>
#include <limits.h>
#include "global.h"
#include "myrand.h"

using namespace std;



bool SHOW_POPULATION = false;
bool SHOW_SELECTION_INDEX = false;
bool SHOW_REPLACEMENT = false;
double R1 = 1.0; // random number for sweeping
double R2 = 1.0; // random number for sweeping
int NUM_CUTPOINTS = 3;
double Bratio = 0.0;

bool isRTR = true;



int numOfReplacement = 0;
MyRand myRand;


void outputErrMsg (const char *errMsg)
{
    printf ("%s\n", errMsg);
    exit (EXIT_FAILURE);
}
