/*****************************************************************************/
/**Author: Gaojie Sun                                                       **/
/**Date: 21/07/2017                                                         **/
/**Supervisor: Dongryeol Ryu                                                **/
/**Originazation: The University of Melbourne                               **/
/*****************************************************************************/

/*****************************************************************************/
/**Documentation Description:                                               **/
/**   This file is used to define some struct to abstract the position data.**/
/*****************************************************************************/
#ifndef _POSITION_DATA_STRUCTURE_H
#define _POSITION_DATA_STRUCTURE_H
#include <stdint.h>

typedef struct Coordinate {
	double longtitude;
	double latitude;
} Coordinate;

typedef struct Position {
	int64_t i, j;
} Position;

typedef struct Polygon {
	Coordinate* coordinates;
	Position* position;
	int pointNumber;
} Polygon;

#endif