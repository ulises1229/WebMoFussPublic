/*
 * costDistance.cpp
 *
 *  Created on: 04/10/2017
 *      Author: lanase
 */

#include "distance.cuh"
#include <math.h>

	// Utility method for comparing two cells

	__device__ bool costDistance::isInsideGrid(int i, int j, float * friction){
		int z = COL*i + j;
		return (i >= 0 && i < ROW && j >= 0 && j < COL && friction[z] >= 0);
	}

	__device__ bool costDistance::isInsideLimit(int i, int j, int source_x, int source_y, int bound_value){
		int max_X = source_x + bound_value/2;
		int min_X = source_x - bound_value/2;
		int max_Y = source_y + bound_value/2;
		int min_Y = source_y - bound_value/2;
		return (i >= min_X && i < max_X && j >= min_Y && j < max_Y);
	}

	__device__ void costDistance::inicio_cost_distance(float*friction, int srcX, int srcY, double projection, bool is_relative, double max, double min, bool flag, bool diagonals_cost_more, float * cost_distance, bool * active_raster, int * active_costs_x, int * active_costs_y, float * active_costs_distance, int uniqueThread, int * rutas, int * dxdy, int bound_value)
	{
		int z = 0;
		int r = 0;
		counter_order = 0;

		for(int x = 0; x < ROW; x++){
			for(int y = 0; y < COL; y++){
				if (is_relative)
				{
					friction[z] = friction[z] * projection;
				}
				z++;
			}
		}
		for(int x = 0; x < bound_value; x++){
			for(int y = 0; y < bound_value; y++){
				r = bound_value*bound_value*uniqueThread + x*bound_value + y;
				cost_distance[r] = INT_MAX;
				active_raster[r] = false;
//				r = 0;
			}
		}

		//the position of a cell that you want to display its neighbours
		int srcX_2 = (int)bound_value/2;
		int srcY_2 = (int)bound_value/2;
		int idraster = bound_value*bound_value*uniqueThread + srcX_2*bound_value + srcY_2;
		active_raster[idraster] = 1;

		//se obtienen los vecinos proximos al origen y sus distancias calculadas. ordenas de menor a mayor
		absolute_size = 0;
		acumulados(active_costs_x, active_costs_y, active_costs_distance, srcX, srcY, max, min,friction, diagonals_cost_more, cost_distance, active_raster, uniqueThread, rutas, dxdy, bound_value);//metodo para calcular demas vecinos.
		int za = 0;
		for(int r=0; r<bound_value; r++){
			for(int c=0; c<bound_value; c++){
				za = bound_value*bound_value*uniqueThread + r*bound_value + c;
				if(cost_distance[za] == INT_MAX)
					cost_distance[za] = 0;
			}
		}

	}

	__device__ void costDistance::acumulados(int * active_costs_x, int * active_costs_y, float * active_costs_distance, int origen_x, int origen_y, double max, double min, float*friction, bool diagonals_cost_more, float * cost_distance, bool * active_raster, int uniqueThread, int * rutas, int * dxdy, int bound_value)
	{
		int cont = 1;
		int z = 0;
		float sqtwo = (float)sqrtf(2);
		int origen_z = origen_x*COL + origen_y;
		float distance = INT_MAX;
		kx = origen_x;
		ky = origen_y;
		//los k del array redimensionado
		kx_2 = bound_value/2;
		ky_2 = bound_value/2;
		kdistance = 0;
		int x = 0, y = 0, x2 = 0, y2 = 0, z2 = 0;
		//PRIMERA CELDA
		for (int i = 0; i < 8; i++)
		{
			//Control sobre fricc
			x = kx + dxdy[i];
			y = ky + dxdy[8 + i];
			z = y + x*COL;

			//Control sobre costo distancia y costos activos
			x2 = kx_2 + dxdy[i];
			y2 = ky_2 + dxdy[8 + i];
			z2 = y2 + x2*bound_value;

			if ( (!isInsideGrid(x, y, friction)) && (!isInsideLimit(x, y, origen_x, origen_y, bound_value)) )
				continue;
			if(i % 2 == 0)
				distance = (friction[ky + kx*COL] + friction[z])/2;
			else
				distance = sqtwo*(friction[ky + kx*COL] + friction[z])/2;
			insert(active_costs_x, active_costs_y, active_costs_distance, x, y, distance, uniqueThread, bound_value);
			int idraster = ROW*COL*uniqueThread + z;
			int idraster2 = bound_value*bound_value*uniqueThread + z2;
			active_raster[idraster2] = true;
			if( distance < cost_distance[bound_value*bound_value*uniqueThread + z2] )
				cost_distance[idraster2] = distance;
		}
		distance = 0;
		cost_distance[bound_value*bound_value*uniqueThread + ky_2 + kx_2*bound_value] = 0;
		//RESTO DE LAS CELDAS
		while( absolute_size>0 )
		{
			begin(active_costs_x, active_costs_y, active_costs_distance, uniqueThread, origen_x, origen_y, bound_value);//inicio con la distancia menor
			erase_begin(active_costs_x, active_costs_y, active_costs_distance, uniqueThread, bound_value);
				cont++;
				int x, y, x_2, y_2, z_2;
				for (int i = 0; i < 8; i++){
					x = kx + dxdy[i];
					y = ky + dxdy[i+8];
					z = COL*x + y;

					//De los redimensionados
					x_2 = kx_2 + dxdy[i];
					y_2 = ky_2 + dxdy[i+8];
					z_2 = bound_value * x_2 + y_2;

					index_float = 0;

					if( (isInsideGrid(x,y, friction)) && (isInsideLimit(x, y, origen_x, origen_y, bound_value)) ) {
						if(i % 2 != 0){// si es movimiento diagonal
							int idraster = ROW*COL*uniqueThread + z;
							int idraster_2 = bound_value*bound_value*uniqueThread + z_2;
							if((x != origen_x || y != origen_y) && active_raster[idraster_2] == 0){
								float dist=0;
								int cont_x = kx; int cont_y = ky;
								int cont_x2 = kx; int cont_y2 = ky;
								int caminos = 1;
								float distancias1=0;
								float distancias2=0;
								float distancias3=0;
								int count = 0;
								for(int z = 0; z < 8; z=z+2){
									int mov_x = rutas[i*8 + z];
									int mov_y = rutas[i*8 + z+1];
									int new_x = cont_x+mov_x;
									int new_y = cont_y+mov_y;
									int new_z = COL*new_x + new_y;
									int bef_z = COL*cont_x2 + cont_y2;
									dist += (friction[new_z] + friction[bef_z])/2;
									if(caminos == 2 || caminos == 4) {
										if(count == 0)
										{
											distancias1 = dist+cost_distance[bound_value*bound_value*uniqueThread + ky_2 + kx_2*bound_value];
											count++;
										}
										else
										{
											distancias2 = dist+cost_distance[bound_value*bound_value*uniqueThread + ky_2 + kx_2*bound_value];
											count = 0;
										}
										cont_x = cont_x2 = kx;
										cont_y = cont_y2 = ky;
										caminos = 0;
										dist = 0;
									} else {
										cont_x2 = cont_x+mov_x;
										cont_y2 = cont_y+mov_y;
										cont_x = cont_x+mov_x;
										cont_y = cont_y+mov_y;
									}
									caminos++;
								}
								if (diagonals_cost_more)
									distancias3 = sqtwo * (friction[ky + kx*COL]+friction[z]) / 2 + cost_distance[bound_value*bound_value*uniqueThread + ky_2 + kx_2*bound_value];
								else
									distancias3 = ((friction[ky + kx*COL]+friction[z]) / 2) + cost_distance[bound_value*bound_value*uniqueThread + ky_2 + kx_2*bound_value];

									float minimo = calculateMin(distancias1, distancias2, distancias3);
									distancias1 = distancias2 = distancias3 = 0;

									if(minimo < cost_distance[bound_value*bound_value*uniqueThread + y_2 + x_2*bound_value] && minimo >= 0){
									bool validate_min = validate_limit(max, min, friction[z]);
										if(validate_min)
										{
												cost_distance[bound_value*bound_value*uniqueThread + y_2 + x_2*bound_value] = minimo;
												insert( active_costs_x, active_costs_y, active_costs_distance, x, y, minimo, uniqueThread, bound_value);
										}
								}
							}
						}
						else{
							float dist = 0;
							int idraster = 0;
							if((x != origen_x || y != origen_y )&& active_raster[idraster] == 0  && isInsideGrid(x,y, friction) && isInsideLimit(x, y, origen_x, origen_y, bound_value)){
								dist = ((friction[ky+ COL * kx] + friction[z])/2) + cost_distance[bound_value*bound_value*uniqueThread + ky_2 + kx_2*bound_value];
								if(dist < cost_distance[bound_value*bound_value*uniqueThread + y_2 + x_2*bound_value] && dist >= 0){
									bool validate_dist = validate_limit(max, min, friction[z]);
									if(validate_dist)
									{
										cost_distance[bound_value*bound_value*uniqueThread + y_2 + x_2*bound_value] = dist;
										insert( active_costs_x, active_costs_y, active_costs_distance, x, y, dist, uniqueThread, bound_value);
									}
								}
							}

						}
					}
				}

				int idraster = bound_value*bound_value*uniqueThread + kx_2*bound_value + ky_2;
				active_raster[idraster] = true;
			}
	}

	__device__ bool costDistance::validate_limit(double max, double min, float friction) {
		bool res = false;
		float limit = (float)(max - min)*0.75;
		if(friction<=limit)
			res = true;
		return res;
	}

	//Algoritmos de la clase stl para estructura set

	__device__ void costDistance::insert(int * active_costs_x, int * active_costs_y, float * active_costs_distance, int x, int y, float distance, int uniqueThread, int bound_value)
	{
		int init = uniqueThread*bound_value*bound_value;
		active_costs_x[init + absolute_size] = x;
		active_costs_y[init + absolute_size] = y;
		active_costs_distance[init + absolute_size] = distance;
		absolute_size++;
		counter_order++;
		if(counter_order%100==0)
			sort(active_costs_x, active_costs_y, active_costs_distance, absolute_size, uniqueThread, bound_value);
	}

	__device__ void costDistance::sort(int * active_costs_x, int * active_costs_y, float * active_costs_distance, int tam, int uniqueThread, int bound_value)
	{
		int init = uniqueThread*bound_value*bound_value;
		float temp = 0;
		int j = 0, temp_x = 0, temp_y = 0, i = 0;
		for (i=init; i<init+tam; i++)
			{
				temp = active_costs_distance[i];
				temp_x = active_costs_x[i];
				temp_y = active_costs_y[i];
				j = i - 1;
				while ( (active_costs_distance[j] > temp) && (j >= init) )
				{
					active_costs_distance[j+1] = active_costs_distance[j];
					active_costs_x[j+1] = active_costs_x[j];
					active_costs_y[j+1] = active_costs_y[j];
					j--;
				}
				active_costs_distance[j+1] = temp;
				active_costs_x[j+1] = temp_x;
				active_costs_y[j+1] = temp_y;
			}
	}

	__device__ void costDistance::begin(int * active_costs_x, int * active_costs_y, float * active_costs_distance, int uniqueThread, int src_x, int src_y, int bound_value)
	{
		int init = uniqueThread*bound_value*bound_value;
		kx = active_costs_x[init];
		ky = active_costs_y[init];
		kdistance = active_costs_distance[init];
		//se determina la posicion de la esq superior izquierda basado en centroide
		int ini_x = src_x - bound_value/2;
		int ini_y = src_y - bound_value/2;
		//actualizo valores de kx_2 y ky_2
		kx_2 = active_costs_x[init] - ini_x;
		ky_2 = active_costs_y[init] - ini_y;
	}

	__device__ void costDistance::erase_begin(int * active_costs_x, int * active_costs_y, float * active_costs_distance, int uniqueThread, int bound_value)
	{
		int init = uniqueThread*bound_value*bound_value;
		active_costs_distance[init] = 0;
		active_costs_x[init] = 0;
		active_costs_y[init] = 0;
		for(int i = init+1; i < init+absolute_size; i++)
		{
			active_costs_distance[i-1] = active_costs_distance[i];
			active_costs_x[i-1] = active_costs_x[i];
			active_costs_y[i-1] = active_costs_y[i];

			if(i == absolute_size-1)
			{
				active_costs_distance[i] = 0;
				active_costs_x[i] = 0;
				active_costs_y[i] = 0;
			}

		}
		absolute_size--;
	}

	  __device__ float costDistance::calculateMin(float dis1, float dis2, float dis3)
		{
			float min=0;
			if ( dis1 <= dis2 && dis1 <= dis3)
				min = dis1;
			else if (dis2 <= dis1 && dis2 <=dis3)
				min = dis2;
			else
				min = dis3;
			return min;
		}
