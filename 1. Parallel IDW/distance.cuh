
/*
* distance.h
*
*  Created on: 15/03/2018
*      Author: lanase
*/

#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <iterator>
#include <set>
#include <climits>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

/*
struct cell {
unsigned short int x, y;
float distance;
//  cell(int x, int y, float distance) :
//      x(x), y(y), distance(distance) {}
};
*/

class costDistance
{
  public:
  //structs
/*  struct cell2 {
  int x, y;
  cell2(int x, int y) :
      x(x), y(y) {}
  };
*/
/*
  struct cell {
    int x, y;
    float distance;
  }
  struct cell cell;
  */
/*
  struct cell {
  unsigned short int x, y;
  float distance;
//  cell(int x, int y, float distance) :
//      x(x), y(y), distance(distance) {}
  };
*/
  //methods
  __device__ bool isInsideGrid(int i, int j, float * friction);
  __device__ bool isInsideLimit(int i, int j, int source_x, int source_y, int bound_value);
  __device__ inline void acumulados(int * active_costs_x, int * active_costs_y, float * active_costs_distance, int origen_x, int origen_y, double max, double min, float* grid, bool diagonals_cost_more, float * cost_distance, bool * active_raster, int uniqueThread, int * rutas, int * dxdy, int bound_value);
  __device__ void inicio_cost_distance(float* grid, int srcX, int srcY, double projection, bool is_relative, double max, double min, bool flag, bool diagonals_cost_more, float * cost_distance, bool * active_raster, int * active_costs_x, int * active_costs_y, float * active_costs_distance, int uniqueThread, int * rutas, int * dxdy, int bound_value);
  __device__ void freeMem();
  __device__ bool validate_limit(double max, double min, float friction);
  __device__ void insert(int * active_costs_x, int * active_costs_y, float * active_costs_distance, int x, int y, float distance, int uniqueThread, int bound_value);
	__device__ void sort(int * active_costs_x, int * active_costs_y, float * active_costs_distance, int tam, int uniqueThread, int bound_value);
	__device__ void begin(int * active_costs_x, int * active_costs_y, float * active_costs_distance, int uniqueThread, int src_x, int src_y, int bound_value);
	__device__ void erase_begin(int * active_costs_x, int * active_costs_y, float * active_costs_distance, int uniqueThread, int bound_value);
  __device__ float calculateMin(float dis1, float dis2, float dis3);

  //variables
  int ROW, COL;
//  int** locations;
//  float** output_raster;
//  float** cost_raster;
//  bool** active_raster;
//  int rutas[8][8]{
//  {0,0,0,0,0,0,0,0},
//  {-1,0,0,1,0,1,-1,0},
//  {0,0,0,0,0,0,0,0},
//  {0,1,1,0,1,0,0,1},
//  {0,0,0,0,0,0,0,0},
//  {0,-1,1,0,1,0,0,-1},
//  {0,0,0,0,0,0,0,0},
//  {-1,0,0,-1,0,-1,-1,0}
//  };

  int absolute_size, index_float;
  int kx;
  int ky;
  int kx_2;
  int ky_2;
  float kdistance;
  int counter_order;
};
