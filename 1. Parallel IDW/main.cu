/*
 * main.cpp
 *
 *  Created on: 17/08/2017
 *      Author: lanase
 */

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
//#include <float.h>
#include <vector>
#include <iterator>
#include <map>
#include "DisplayImage.cu"
#include "distance.cuh"

#include "tclap/CmdLine.h"
#include "CudaCostDistance.cuh"
#include "gputimer.cuh"
#include <omp.h>
#include <boost/algorithm/string.hpp>
/*
struct centroid {
	int x;
	int y;
};

*/
char is_usable;

void compare_matrices(float ** m_friction, float ** m_fricged, int r, int c)
{
  int count = 0;
  ofstream ofs("matrix_diff.txt", std::ofstream::app);
  for(int i = 0; i < r; i++)
  {
    for(int j = 0; j < c; j++)
    {
      if((m_friction[i][j] > m_fricged[i][j]) && (m_fricged[i][j] >= 0 && m_friction[i][j] >= 0))
      {
        ofs << "("<< i <<","<< j << "), F=" << setprecision(10) << m_friction[i][j] << ", G = " << m_fricged[i][j] << ", " << m_friction[i][j] - m_fricged[i][j] << endl;
        count++;
      }
    }
  }
  ofs << count;
  ofs.close();
}
/*
void SequentialExec(vector<centroid> localities, costDistance d, Display_image di, float** m_friction,
					bool is_relative, double max, double min, bool diagonals_cost_more, int rows, int cols)
{
	//Func dist is called from the localities vector
	for (int k = 0; k < localities.size(); k++)
	{
		int coord_x = localities.at(k).x;
		int coord_y = localities.at(k).y;

		clock_t begin_cd = clock();

		d.inicio_cost_distance(m_friction, coord_x, coord_y, di.projection, is_relative, max, min, false, diagonals_cost_more);

		//d.inicio_cost_distance(d.output_raster, coord_x, coord_y, di.projection, is_relative, max, min, true);
		//for (int i = 0; i < 25; i++) {
		//d.inicio_cost_distance(d.output_raster, coord_x, coord_y, di.projection, is_relative, max, min, true);
		//}

		clock_t end_cd = clock();
		double elapsed_secs_cd = double(end_cd - begin_cd) / CLOCKS_PER_SEC;
		cout << "Cost Distance time = " << elapsed_secs_cd << " secs." << endl;

		clock_t begin3 = clock();
		di.matrix_to_tiff(d.output_raster, rows, cols, k + 1);
		clock_t end3 = clock();
		double elapsed_secs3 = double(end3 - begin3) / CLOCKS_PER_SEC;
		cout << "Creating the final route image took " << elapsed_secs3 << " secs." << endl;

		//compare_matrices(d.output_raster, m_fricged, rows, cols);
		d.freeMem();
	}
}
*/
int main(int argc, const char** argv) {
  string map_friction, map_locality, is_r, ex_m, d_c, file_consumption /*map_fric_gdal,*/;
  bool is_relative = false;
  bool diagonals_cost_more = false;
  float bound_value = 0;
  float exp_value = 0;
  //Command line messages
  TCLAP::CmdLine cmd("Command description message", ' ', "1");
  TCLAP::ValueArg<std::string> friction("f","friction","Absolute path to friction_map.tif",true,"/path/to/image.tif","string");
  TCLAP::ValueArg<std::string> locality("l","locality","Absolute path to locality_map.tif",true,"/path/to/image.tif","string");
  TCLAP::ValueArg<std::string> consumption("c","consumption","Absolute path to consumption_file.csv",true,"/path/to/file.csv","string");
  //TCLAP::ValueArg<std::string> gdal("g","gdal","Absolute path to gdal_map.tif",true,"/path/to/image.tif","string");
  TCLAP::ValueArg<std::string> diagonals_cost("d","diagonals_cost","If true, diagonals have extra cost.",true,"true or false","string");
  TCLAP::ValueArg<std::string> is_rel("r","relative","Determines if the friction on map is relative",true,"true or false","string");
  TCLAP::ValueArg<std::string> ex_mod("e", "execution", "Program execution mode (sequential, OpenMP or CUDA)", true, "pseq, pomp or pcuda", "string");
  TCLAP::ValueArg<std::string> bo_val("t", "time", "Maximum analysis limit (time hours)", true, "enter a number", "string");
  TCLAP::ValueArg<std::string> e_val("x", "exponent", "Exponent value", true, "enter a number", "string");

  // Add arguments to command line output
  cmd.add(friction);
  cmd.add(locality);
  cmd.add(consumption);
  cmd.add(diagonals_cost);
  cmd.add(is_rel);
  cmd.add(ex_mod);
  cmd.add(bo_val);
  cmd.add(e_val);
  // Parse the argv array.
  cmd.parse( argc, argv );
  // Get the value parsed by each arg.
  map_friction = friction.getValue();
  map_locality = locality.getValue();
  file_consumption = consumption.getValue();
  //map_fric_gdal = gdal.getValue();
  is_r = is_rel.getValue();
  ex_m = ex_mod.getValue();
  d_c = diagonals_cost.getValue();
  bound_value = atof(bo_val.getValue().c_str());
  exp_value = atof(e_val.getValue().c_str());
  //Check if friction is relative
  if (is_r == "true" || is_r == "t")
    is_relative = true;

  if (d_c == "true" || d_c == "t")
    diagonals_cost_more = true;
  //Printing the input values
  cout << "Friction Map: " << map_friction << endl;
  cout << "Is Relative?: " << is_relative << endl;

  //Is usable set value
  is_usable = 'n';
  // Declare image object
  Display_image di;
  // Import friction and points
  // Import biomass and friction
	clock_t begin = clock();
  float** m_friction = di.tiff_to_matrix_gdal(map_friction, true);
  float** m_locality = di.tiff_to_matrix_gdal(map_locality, false);
  //float** m_fricged = di.tiff_to_matrix_gdal(map_fric_gdal, false);

  //cout << "Done" << endl;
	cout << di.epsg << endl;
	di.reproject_coords(map_friction);
	//exit(0);
	clock_t end = clock();
  int rows = di.ROWS;
  int cols = di.COLS;
  cout << "rows: " << rows;
  cout << "cols: " << cols;
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout << "TIFF to matrix = " << elapsed_secs << " secs." << endl;
  int local_count = 0;

  //Load file of distance limit

  std::ifstream file(file_consumption);
  cout << "Archivo de consumo: " << file_consumption << endl;
  std::vector<std::vector<std::string> > dataList;
  std::string line = "";
  // Iterate through each line and split the content using delimeter
  while (getline(file, line))
  {
    std::vector<std::string> vec;
    boost::algorithm::split(vec, line, boost::is_any_of(","));
    dataList.push_back(vec);
    cout << "Columnas del arreglo: " << dataList.back().size() << endl;
    cout << "Elemento insertado: " << dataList.back().at(1) << ", " << dataList.back().at(1) ;
  }
  cout << "Tamaño del vector: " << dataList.size() << endl;
  // Close the File
  file.close();
  int size_localities = dataList.size() - 1;

  //Generate vector of localities
  vector<centroid> localities;
  for(int i=0; i<rows; i++)
  {
    for(int j=0; j<cols; j++)
    {
     if(m_locality[i][j] >= 1 && m_locality[i][j] <= size_localities)
      {
        localities.push_back({i, j});
        local_count++;
      }
    }
  }
  cout << "termina de añadir localidades: " << local_count << endl;

  //Function distance is declared
  costDistance d;
  d.COL = cols;
  d.ROW = rows;

  double min = di.adfMinMax[0];
  double max = di.adfMinMax[1];

  double proj = di.projection;
  int resized_bound = ceil( ( (bound_value*3600) / (proj*min) )*2 );

  cout << "objeto por crear" <<  endl;

  CudaCostDistance obj(localities, resized_bound, resized_bound, rows, cols);
  cout << "objeto creado" << endl;

  //	Execution mode selector

  if (ex_m == "pseq")
  {
	  //	Sequential func dist is called
	//  SequentialExec(localities, d, di, m_friction, is_relative, max, min, diagonals_cost_more, rows, cols);
	  cout << "Secuencial execution finished" << endl;
  }
  else if (ex_m == "pomp")
  {
	  //	Parallel execution with OpenMP is called
  }
  else if (ex_m == "pcuda")
  {
	  //	Parallel execution with CUDA is called
	  //	CudaCostDistance obj(localities, rows, cols);
	  RunDistance(obj, di, m_friction, is_relative, max, min, diagonals_cost_more,rows, cols, resized_bound, file_consumption, exp_value, dataList);
	  cout << "Cuda execution finished" << endl;
  }
  else
	  cout << "Undefined option selected" << endl;
}
