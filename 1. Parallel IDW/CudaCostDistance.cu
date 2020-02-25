/*
============================================================================
Name        : CudaCostDistance.cu
Author      : Rodrigo Gutierrez Vazquez | LANASE
Version     : 1.0
Copyright   : @2019
Description : CudaCostDistance methods
============================================================================
*/
//#pragma once

#include "CudaCostDistance.cuh"

#include "device_launch_parameters.h"
#include "cuda.h"
#include <thrust/device_vector.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>

#include "distance.cuh"
#include "DisplayImage.cu"
#include "gputimer.cuh"

//struct centroid {
//	int x;
//	int y;
//};

CudaCostDistance::CudaCostDistance(vector<centroid> localities_a, int rows, int cols, int real_rows, int real_cols)
{
	cout << "Valor de bound: " << rows << " x " << cols << endl;
	cudaError err = cudaGetLastError();
			if ( cudaSuccess != err )
			{
				cout << "Ya entro con error" << endl;
					fprintf( stderr, "cudaCheckError() failed at : %s\n"
									 , cudaGetErrorString( err ) );
					exit( -1 );
			}
	numLocalities = localities_a.size();

	locality = CopyVectorToArray(localities_a);
	cudaError err10 = cudaGetLastError();
	    if ( cudaSuccess != err10 )
	    {
				cout << "Error 10" << endl;
	        fprintf( stderr, "cudaCheckError() failed at : %s\n"
	                 , cudaGetErrorString( err10 ) );
	        exit( -1 );
	    }

	gridBlockSize = calculateGridBlockSize(numLocalities);
	gridBlockSize2 = calculateGridBlockSize(real_rows*real_cols);

	output_matrices = createOutputMatrix(rows, cols);
	cudaError err1 = cudaGetLastError();
	    if ( cudaSuccess != err1 )
	    {
				cout << "Error 1" << endl;
	        fprintf( stderr, "cudaCheckError() failed at : %s\n"
	                 , cudaGetErrorString( err1 ) );
	        exit( -1 );
	    }
	active_raster = createActiveMatrix(rows, cols);
	cudaError err2 = cudaGetLastError();
	    if ( cudaSuccess != err2 )
	    {
				cout << "Error 2" << endl;
	        fprintf( stderr, "cudaCheckError() failed at : %s\n"
	                 , cudaGetErrorString( err2 ) );
	        exit( -1 );
	    }

//	active_costs = createActiveCosts(rows, cols);
	active_costs_x = createActiveCostsXY(rows, cols);
	cudaError err3 = cudaGetLastError();
			if ( cudaSuccess != err3 )
			{
				cout << "Error 3" << endl;
					fprintf( stderr, "cudaCheckError() failed at : %s\n"
									 , cudaGetErrorString( err3 ) );
					exit( -1 );
			}
	active_costs_y = createActiveCostsXY(rows, cols);
	cudaError err4 = cudaGetLastError();
			if ( cudaSuccess != err4 )
			{
				cout << "Error 4" << endl;
					fprintf( stderr, "cudaCheckError() failed at : %s\n"
									 , cudaGetErrorString( err4 ) );
					exit( -1 );
			}
	active_costs_distance = createActiveCostsDistance(rows, cols);
	cudaError err5 = cudaGetLastError();
			if ( cudaSuccess != err5 )
			{
				cout << "Error 5" << endl;
					fprintf( stderr, "cudaCheckError() failed at : %s\n"
									 , cudaGetErrorString( err5 ) );
					exit( -1 );
			}
	distancias_float = createFloatCosts();
	cudaError err6 = cudaGetLastError();
			if ( cudaSuccess != err6 )
			{
				cout << "Error 6" << endl;
					fprintf( stderr, "cudaCheckError() failed at : %s\n"
									 , cudaGetErrorString( err6 ) );
					exit( -1 );
			}
	rutas = createRutas();
	cudaError err7 = cudaGetLastError();
			if ( cudaSuccess != err7 )
			{
				cout << "Error 7" << endl;
					fprintf( stderr, "cudaCheckError() failed at : %s\n"
									 , cudaGetErrorString( err7 ) );
					exit( -1 );
			}
	dxdy = CreateDxDy();
	cudaError err8 = cudaGetLastError();
			if ( cudaSuccess != err8 )
			{
				cout << "Error 8" << endl;
					fprintf( stderr, "cudaCheckError() failed at : %s\n"
									 , cudaGetErrorString( err8 ) );
					exit( -1 );
			}

	rows2 = rows;
	cols2 = cols;
}

void RunDistance(CudaCostDistance obj, Display_image di, float** m_friction,
					bool is_relative, double max, double min, bool diagonals_cost_more, int rows, int cols, int bound_value, string consumption, float exp_value, std::vector<std::vector<std::string> > dataList)
{
	//	Num. de bloques en (x,y,z) por grid
	dim3 dimGrid(obj.gridBlockSize[0]);
	//	NUm. de hilos en (x, y, z) por bloque
	dim3 dimBlock(obj.gridBlockSize[1]);
	//	Num. de bloques en (x,y,z) por grid
	dim3 dimGrid2(obj.gridBlockSize2[0]);
	//	NUm. de hilos en (x, y, z) por bloque
	dim3 dimBlock2(obj.gridBlockSize2[1]);
	//Kernel

	//Timer para medir tiempo de ejecucion del GPU
	GpuTimer timer;

	//Se declara el array en device para hacer el pre-conteo de friccion

	//Se copia el array de friccion a dos diferentes objetos: uno para el kernel y otro para su uso normal de ordenamiento
	float * friction = obj.CopyMatrixToArray(m_friction, rows, cols);

	float * consumption_vec = obj.CreateConsumptionArray(dataList);

	cudaError err11 = cudaGetLastError();
			if ( cudaSuccess != err11 )
			{
				cout << "Error 11" << endl;
					fprintf( stderr, "cudaCheckError() failed at : %s\n"
									 , cudaGetErrorString( err11 ) );
					exit( -1 );
			}

	//Se manda llamar la rutina para ordenar el arreglo y obtener el valor minimo de friccion del raster completo

//	int radius = ceil( bound_value / min );
//	int bound_va = radius*2;

	cout << "antes de entrar a kernel" << endl;
	timer.Start();
	CostDistanceKernel << < dimGrid, dimBlock >> >(obj, di, friction, is_relative, max, min, diagonals_cost_more, rows, cols, bound_value);
	cudaDeviceSynchronize();
	cout << "despues de entrar a kernel" << endl;
	cudaError err = cudaGetLastError();
	    if ( cudaSuccess != err )
	    {
	        fprintf( stderr, "cudaCheckError() failed at : %s\n"
	                 , cudaGetErrorString( err ) );
	        exit( -1 );
	    }

	cout << "Lanza el segundo kernel para la sumatoria " << endl;
	SumDistancesKernel << < dimGrid2, dimBlock2 >> >(obj, di, friction, rows, cols, bound_value, consumption_vec);
	cout << "despues de entrar a segundo kernel" << endl;
	cudaError err2 = cudaGetLastError();
	    if ( cudaSuccess != err2 )
	    {
	        fprintf( stderr, "cudaCheckError() failed at : %s\n"
	                 , cudaGetErrorString( err ) );
	        exit( -1 );
	    }
	timer.Stop();
 	float timeValue = timer.Elapsed();
	printf("Tiempo de ejecucion en paralelo: %f ms\n", timeValue);
	//Copiar vector de resultados de device a cudaMemcpyHostToDevice
	float ** raster_output = obj.CopyFinalOutputToHost(friction, rows, cols);

	di.matrix_to_tiff(raster_output, rows, cols, 2);
	//Generar la suma de matrices finales

}

__global__ void SumDistancesKernel(CudaCostDistance obj, Display_image di, float* m_friction, int rows, int cols, int bound_value, float * consumption)
{
	int bx = blockIdx.x;
	int by = blockIdx.y;
	int bz = blockIdx.z;
	int thx = threadIdx.x;
	int thy = threadIdx.y;
	int thz = threadIdx.z;
	int nThreads = blockDim.x * blockDim.y * blockDim.z;
	int idThread = (thy*blockDim.x + thx) + (blockDim.x * blockDim.y)*thz;
	int idBlock = (by*gridDim.x + bx) + (gridDim.x * gridDim.y)*bz;
	int uniqueThread = nThreads*idBlock + idThread;
	if (uniqueThread < rows*cols)
	{

		m_friction[uniqueThread] = 0;

		//Se divide el posicion uniqueThread en los correspondientes i j, posicion en matrix rectangular
		int i = uniqueThread / cols;
		int j = uniqueThread - i*cols;

		//Itera sobre las localidades
		for(int c=0; c < obj.numLocalities; c++)
		{

			//Obtiene origen de la localidad
			int origen_x = obj.locality[c * 2];
			int origen_y = obj.locality[c * 2 + 1];

			//Delimita el area de exploracion dentro de la localidad
			int max_X = origen_x + bound_value/2;
			int min_X = origen_x - bound_value/2;
			int max_Y = origen_y + bound_value/2;
			int min_Y = origen_y - bound_value/2;

			//	Verifica si coordenadas caen dentro del area de exploracion
			if (i >= min_X && i < max_X && j >= min_Y && j < max_Y)
			{
				//Calculo la distancia entre el origen y el actual del mapa de friccion
				int distancia_i = origen_x - i;
				int distancia_j = origen_y - j;

				//En el mapa de costos, sumo la distancia anterior, a partir del centro, para calcular la posicion correspondiente al de costos
				int i_cost = bound_value/2 + distancia_i;
				int j_cost = bound_value/2 + distancia_j;

				//Calculo la posicion absoluta dentro del vector de costos
				int z_cost = bound_value*bound_value*c + j_cost + i_cost*bound_value;

				//Realiza la sumatoria
				if( obj.output_matrices[z_cost] > 0 )
					m_friction[uniqueThread] += consumption[c]/pow(obj.output_matrices[z_cost], 1);
				//printf("Valor de consumo: %f\n", consumption[c]);
				//printf("[%f]\n", obj.output_matrices[z_cost]);
				//printf("Valor de parte baja de division: %f\n", pow(obj.output_matrices[z_cost], 1));
				//printf("Valor final de la division: %f\n", consumption[c]/pow(obj.output_matrices[z_cost], 1));
				//printf("Valor actual de la casilla (%d, %d) en sumatoria: %f\n", i, j, m_friction[uniqueThread]);
			}

		}

	}
}

__global__ void CostDistanceKernel(CudaCostDistance obj, Display_image di, float* m_friction, bool is_relative, double max, double min, bool diagonals_cost_more, int rows, int cols, int bound_value)
{
	int bx = blockIdx.x;
	int by = blockIdx.y;
	int bz = blockIdx.z;
	int thx = threadIdx.x;
	int thy = threadIdx.y;
	int thz = threadIdx.z;
	int nThreads = blockDim.x * blockDim.y * blockDim.z;
	int idThread = (thy*blockDim.x + thx) + (blockDim.x * blockDim.y)*thz;
	int idBlock = (by*gridDim.x + bx) + (gridDim.x * gridDim.y)*bz;
	int uniqueThread = nThreads*idBlock + idThread;
	if (uniqueThread < obj.numLocalities)
	{
		costDistance d;
		d.COL = cols;
		d.ROW = rows;
		int coord_x = obj.locality[uniqueThread * 2];
		int coord_y = obj.locality[uniqueThread * 2 + 1];
//		costDistance::cell * active_cos = (costDistance::cell *)((char*)obj.active_costs + uniqueThread * obj.pitch_3);
		d.inicio_cost_distance(m_friction, coord_x, coord_y, di.projection, is_relative, max, min, false, diagonals_cost_more, obj.output_matrices, obj.active_raster, obj.active_costs_x, obj.active_costs_y, obj.active_costs_distance, uniqueThread, obj.rutas, obj.dxdy, bound_value);
	}
}

int * CudaCostDistance::calculateGridBlockSize(int requiredThreads)
{
	int i = 0;
	int numThreads = 0;
	int * gridBlock = new int[2];
	float numBlock = 0.0;
	for (i = 1; numThreads <= requiredThreads; i++)
	{
		numThreads = i * 32;
	}
	numBlock = 512.0;
	gridBlock[0] = ceil(numThreads / numBlock);
	gridBlock[1] = (int)numBlock;
	return gridBlock;
}

int * CudaCostDistance::CopyVectorToArray(std::vector<centroid> localities)
{
	int * locali_array = (int *)malloc(numLocalities * sizeof(int) * 2);
	int * locali_device;
	for (int i = 0; i < numLocalities; i++)
	{
		locali_array[i*2] = localities.at(i).x;
		locali_array[i*2+1] = localities.at(i).y;
	}
	cudaMalloc((void **)&locali_device, numLocalities * sizeof(int) * 2);
	cudaMemcpy(locali_device, locali_array, numLocalities * sizeof(int) * 2, cudaMemcpyHostToDevice);
	return locali_device;
}

float * CudaCostDistance::CreateConsumptionArray(std::vector<std::vector<std::string> > dataList)
{
//		cout << "Creando valor de consumo. " << endl;
		size_localities = dataList.size() - 1;
		float * consumption = (float *)malloc(size_localities*sizeof(double));
		for(int i = 1, j=0; i < size_localities + 1 ; i++, j++)
		{
			consumption[j] = 0;
//			cout << "Valor de consumo en array: " << dataList.at(i).at(1) << endl;
			consumption[j] = (float) std::stof(dataList.at(i).at(1));
//			consumption[i] = atof(dataList.at(i).at(1).c_str());
//			cout << "Valor de consumo: " << consumption[j] << endl;
		}

		float * consumption_device;
		cudaMalloc((void **)&consumption_device, size_localities * sizeof(float));
		cudaMemcpy(consumption_device, consumption, size_localities * sizeof(float), cudaMemcpyHostToDevice);
		return consumption_device;

}

		float * CudaCostDistance::CopyMatrixToArray(float ** friction_m, int rows, int cols)
		{
			int z = 0;
			float * friction_host = new float[rows*cols];
			for(int x = 0; x < rows; x++)
			{
				for(int y = 0; y < cols; y++)
				{
					friction_host[z] = friction_m[x][y];
		//			DevFric[z] = friction_m[x][y];
					z++;
				}
			}
			float * friction_device;
			cudaMalloc((void **)&friction_device, rows*cols * sizeof(float));
			cudaMemcpy(friction_device, friction_host, rows*cols * sizeof(float), cudaMemcpyHostToDevice);
			return friction_device;
}

float CudaCostDistance::CalculateMinFric(float * DevFric, int rows, int cols)
{
	float min = 0;
	thrust::sort(DevFric, DevFric + rows*cols);
	for(int i = 0; i < rows*cols; i++)
	{
		if( DevFric[i]>0 )
		{
			min = DevFric[i];
			break;
		}
	}
	return min;
}

float * CudaCostDistance::createOutputMatrix(int rows, int cols)
{
	float * output_device;
	size_t width = cols * rows * sizeof(float);
	size_t height = numLocalities;
	cudaMalloc((void **)&output_device, width*height);
	return output_device;
}

bool * CudaCostDistance::createActiveMatrix(int rows, int cols)
{
	bool * output_device;
	size_t width = cols * rows * sizeof(bool);
	size_t height = numLocalities;
	cudaMalloc((void **)&output_device, width*height);
	return output_device;
}
/*
cell * CudaCostDistance::createActiveCosts(int rows, int cols)
{
	cell * output_device;
	size_t width = cols * sizeof(cell);
	size_t height = rows;
	cudaMalloc((void **)&output_device, width*height);
	return output_device;
}
*/
int * CudaCostDistance::createActiveCostsXY(int rows, int cols)
{
	int * output_device;
	size_t width = cols * sizeof(int) * rows;
	size_t height = numLocalities;
	cudaMalloc((void **)&output_device, width*height);
	return output_device;
}

float * CudaCostDistance::createActiveCostsDistance(int rows, int cols)
{
	float * output_device;
	size_t width = cols * sizeof(float) * rows;
	size_t height = numLocalities;
	cudaMalloc((void **)&output_device, width*height);
	return output_device;
}

float * CudaCostDistance::createFloatCosts()
{
	float * output_device;
	size_t width = 3 * sizeof(float);
	size_t height = numLocalities;
	cudaMalloc((void **)&output_device, width*height);
	return output_device;
}

int * CudaCostDistance::createRutas()
{
	int rutas[64] = {0,0,0,0,0,0,0,0,-1,0,0,1,0,1,-1,0,0,0,0,0,0,0,0,0,0,1,1,0,1,0,0,1,0,0,0,0,0,0,0,0,0,-1,1,0,1,0,0,-1,0,0,0,0,0,0,0,0,-1,0,0,-1,0,-1,-1,0};
	int * rutas_device;
	cudaMalloc((void **)&rutas_device, 8 * 8 * sizeof(int));
	cudaMemcpy(rutas_device, rutas, 8 * 8 * sizeof(int), cudaMemcpyHostToDevice);
	return rutas_device;
}

int * CudaCostDistance::CreateDxDy()
{
	int dxdy[16] = { -1, -1, 0, 1, 1,  1,  0, -1, 0,  1, 1, 1, 0, -1, -1, -1 };
	int * dxdy_device;
	cudaMalloc((void **)&dxdy_device, 8 * 2 * sizeof(int));
	cudaMemcpy(dxdy_device, dxdy, 8 * 2 * sizeof(int), cudaMemcpyHostToDevice);
	return dxdy_device;
}

void CudaCostDistance::DestroyLocalArray(int * localities)
{
	cudaFree(localities);
}

float ** CudaCostDistance::SumCostDistance(float * costo, int rows, int cols)
{
	int z = 0, consumo_lenia = 2, coeficiente = 1;
	float ** output_costo = new float*[rows];
	for(int i = 0; i < rows; i++)
	{
		output_costo[i] = new float[cols];
	}

	for( int j = 0; j < rows; j++)
	{
		for( int k = 0; k < cols; k++ )
		{
			output_costo[j][k] = 0;
		}
	}

	for( int j = 0; j < rows; j++)
	{
		for( int k = 0; k < cols; k++ )
		{
			z = j * cols + k;
			output_costo[j][k] = output_costo[j][k] + consumo_lenia/pow(costo[z], coeficiente);
		}
	}

	return output_costo;

}

float ** CudaCostDistance::CopyFinalOutputToHost(float * friction, int rows, int cols)
{
	float * final_friction = (float * ) malloc(rows*cols*sizeof(float));
	cudaMemcpy(final_friction, friction, rows*cols*sizeof(float), cudaMemcpyDeviceToHost);
	float ** friction_f = (float ** ) malloc(rows*sizeof(float * ));
	int z = 0;
	for(int i=0; i < rows; i++)
	{

		friction_f[i] = (float * ) malloc(cols*sizeof(float));

		for(int j = 0; j < cols; j++)
		{
			friction_f[i][j] = final_friction[z];
//			cout << "[" << friction_f[i][j] << "]" << endl;
			z++;
		}

	}

	return friction_f;
}

void CudaCostDistance::freeArrays()
{
	cudaFree(locality);
	cudaFree(active_raster);
	cudaFree(active_costs_x);
	cudaFree(active_costs_y);
	cudaFree(active_costs_distance);
	cudaFree(distancias_float);
	cudaFree(rutas);
	cudaFree(dxdy);
}
