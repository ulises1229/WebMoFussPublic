/*
============================================================================
Name        : CudaCostDistance.cuh
Author      : Rodrigo Gutierrez Vazquez | LANASE
Version     : 1.0
Copyright   : @2019
Description : Header with class and CudaCostDistance methods
============================================================================
*/

#ifndef CudaCostDistance_cuh
#define CudaCostDistance_cuh

#include <stdio.h>
#include <thrust/device_vector.h>
#include <iostream>

#include "distance.cuh"
#include "DisplayImage.cu"

struct centroid {
	int x;
	int y;
};
/*
struct cell {
unsigned short int x, y;
float distance;
//  cell(int x, int y, float distance) :
//      x(x), y(y), distance(distance) {}
};
*/
class CudaCostDistance
{
public:
	CudaCostDistance(vector<centroid> localities_a, int rows, int cols, int real_rows, int real_cols);
	__device__ void CopyRasterToOutput(costDistance d, int rows, int cols, int uniqueThread);
	float ** CopyFinalOutputToHost(float * friction, int rows, int cols);
	void DestroyLocalArray(int * localities);
	float * CopyMatrixToArray(float ** friction_m, int rows, int cols);
	float * CreateConsumptionArray(std::vector<std::vector<std::string> > dataList);
	float CalculateMinFric(float * DevFric, int rows, int cols);
	float ** SumCostDistance(float * costo, int rows, int cols);
	int * gridBlockSize;
	int * gridBlockSize2;
	int * locality;
	float * output_matrices;
	float * output_raster_host;
	vector<centroid> localities;
	int numLocalities, rows2, cols2;
	size_t pitch, pitch_2, pitch_3;
	bool * active_raster;
//	cell * active_costs;
	int * active_costs_x;
	int * active_costs_y;
	float * active_costs_distance;
	float * distancias_float;
	int * rutas;
	int * dxdy;
	int size_localities;
private:
	int * CopyVectorToArray(std::vector<centroid> localities);
	int * calculateGridBlockSize(int requiredThreads);
	float * createOutputMatrix(int rows, int cols);
	bool * createActiveMatrix(int rows, int cols);
//	cell * createActiveCosts(int rows, int cols);
	int * createActiveCostsXY(int rows, int cols);
	float * createActiveCostsDistance(int rows, int cols);

	float * createFloatCosts();
	int * createRutas();
	int * CreateDxDy();

	void freeArrays();

};

void RunDistance(CudaCostDistance obj, Display_image di, float** m_friction,
					bool is_relative, double max, double min, bool diagonals_cost_more, int rows, int cols, int bound_value,
					string consumption, float exp_value, std::vector<std::vector<std::string> > dataList);
__global__ void CostDistanceKernel(CudaCostDistance obj, Display_image di, float* m_friction, bool is_relative, double max, double min, bool diagonals_cost_more, int rows, int cols, int bound_value);
__global__ void SumDistancesKernel(CudaCostDistance obj, Display_image di, float* m_friction, int rows, int cols, int bound_value, float * consumption);
//__global__ void ExampleKernel(CudaCostDistance obj);

#endif // !CudaCostDistance.cuh
