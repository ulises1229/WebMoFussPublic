//#include <Python.h>
//#include <gdal_utils.h>
//#include "/usr/lib64/python2.7/site-packages/numpy/core/include/numpy/arrayobject.h"
#pragma once

#include <stdlib.h>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/core/mat.hpp>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <array>
#include <memory>
#include <stdexcept>
#include <cstdio>
#include <algorithm>
#include <map>
#include <set>
#include "/usr/include/gdal/gdal.h"
#include "/usr/include/gdal/gdal_priv.h"
#include "/usr/include/gdal/gdalwarper.h"
#include "/usr/include/gdal/ogr_spatialref.h"
#include "/usr/include/gdal/ogr_geometry.h"
 
using namespace cv;
using namespace std;

struct Point2D {
	int x, y;
	Point2D () {};
	Point2D(int x, int y) : x(x), y(y) {}
	bool operator () (Point2D const &e) const {
		return (x == e.x && y == e.y);
	}
};

struct Grid {
	vector<Point2D> elements;
	int noElements = 0;
	int invalidCells = 0;
	int id;
	float sum = 0.0;
	float biomass = 0;
	float friction = 0;
	float value;
	float biomassAvg = 0;

	Grid() {}

	Grid(vector<Point2D> elements, int noElements, int invalidCells, float sum) :
		elements(elements),noElements(noElements), invalidCells(invalidCells), sum(sum) {}

};

struct cellVecinos {
	int x, y;
	float distance;
	cellVecinos(int x, int y, float distance) :
			x(x), y(y), distance(distance) {}
};

// Utility method for comparing two cells
static bool operator<(const cellVecinos& a, const cellVecinos& b){
	if (a.distance == b.distance){
		if (a.x != b.x)
			return (a.x < b.x);
		else
			return (a.y < b.y);
	}
	return (a.distance < b.distance);
}

class Display_image{
public:
	float avg_biomasa, pixels_necesarios = 0, xMax = FLT_MIN, xMin = FLT_MAX, yMax = FLT_MIN, yMin = FLT_MAX;
	int ROWS, COLS, intervals = 0, totValidGrids = 0, totGridsAvg = 0, valid_points = 0, bGotMin, bGotMax;
	bool flag = true;
	double projection, adfGeoTransform[6], adfMinMax[2];
	string proj, epsg;
	vector<Point2D> active_raster;
	float** costos;
	int** costosAux;
	map<float,Grid> gridsMap;
	vector<float> tokens;
	Display_image()
	{
		//nothing to do here
	}

	bool isInsideGrid(int i, int j){
			return (i >= 0 && i < ROWS && j >= 0 && j < COLS);
	}


	float** tiff_to_matrix_gdal(string file, bool flag) {
		GDALDataset *dataset;
		char **MD;
		char *info;
		GDALAllRegister();
		string ds = file;
		dataset = (GDALDataset *) GDALOpen(ds.c_str(), GA_ReadOnly);
		double ndv;


		if(dataset == NULL) {
			cout << "Null dataset" << endl;
			exit(0);
		}

		GDALRasterBand  *poBand;

		poBand = dataset->GetRasterBand(1);

		dataset->GetGeoTransform( adfGeoTransform );
		projection = adfGeoTransform[1];

		proj = dataset->GetProjectionRef();

		cout << projection << endl;

		string tail = proj.substr(proj.size() - 20);

		int d1, d2;

		d1 = tail.find(",");
		d2 = tail.find("\"", d1 + 2);

		epsg = tail.substr(d1 + 2, d2 - d1 - 2);


		int nXSize = poBand->GetXSize();
		int nYSize = poBand->GetYSize();

		ROWS = nYSize; COLS = nXSize;

		//costosAux = new int*[ROWS];
		costos = new float*[ROWS];
		for(int i = 0; i< ROWS; ++i) {
			//costosAux[i] = new int[COLS];
			costos[i] = new float[COLS];
		}

		//adfMinMax[0] = poBand->GetMinimum( &bGotMin );
		//adfMinMax[1] = poBand->GetMaximum( &bGotMax );
		if( flag )
    	GDALComputeRasterMinMax((GDALRasterBandH)poBand, TRUE, adfMinMax);

		cout << fixed << "Min: " << adfMinMax[0] << "\t" << "Max:" << adfMinMax[1] << endl;

		//GDALDataType dataType = poBand->GetRasterDataType();
		float *pBuf = new float[nYSize * nXSize];

		poBand->RasterIO(GF_Read, 0, 0, nXSize, nYSize, pBuf, nXSize, nYSize, GDT_Float32, 0, 0);

		float biomass = 0;
		float secMaxPxVal = 0;
		int cCols = 0, cRows = 0;
		for (int i = 0; i < ROWS; i++) {
			for (int j = 0; j < COLS; j++) {
				int location = (nXSize * (i)) + j;
				if (cCols == COLS - 1) {
					cCols = 0;
					cRows++;
				}
				else {
					costos[cRows][cCols] = *(pBuf+location);
					//if(costos[cRows][cCols])

					if(costos[cRows][cCols] > 0 && flag){
						valid_points++;
						biomass += costos[cRows][cCols];
						if ( (costos[cRows][cCols] > secMaxPxVal) && (costos[cRows][cCols] < adfMinMax[1]) ) {
							secMaxPxVal = costos[cRows][cCols];
							cout << "NEW 2nd Max Val: " << secMaxPxVal << endl;
						}
						//cout << costos[cRows][cCols] << endl;
					}
					cCols++;
				}
			}
		}
		if(flag){
			cout << "The Second Max Val Is: " << fixed << secMaxPxVal << endl;
			//cout << fixed << "Total B: " << biomass << endl;
			//cout << "valid: " << valid_points << endl;
			this->avg_biomasa = biomass / (valid_points);
			cout << "Avg Frict: " << avg_biomasa << endl;
		}
		//cout << avg_biomasa << endl;
		//exit(0);
		return costos;
	}

	void write_image(float** output_raster, int rows, int cols, string heuristic, int stop, string map, string algName) {
		float* arr = new float[rows * cols * 3];
		int n_pixels = rows * cols, channelDiv = round(stop /255);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				arr[i * cols + j] = output_raster[i][j];
				/*if(output_raster[i][j] > 0 && output_raster[i][j] < 0.1)
					arr[i * cols + j] = output_raster[i][j] * 3;
				else if(output_raster[i][j] >= 0.1)
					arr[i * cols + j] = output_raster[i][j];
				else
					arr[i * cols + j] = 0;*/
			}
		}
		Mat matrix = Mat(rows, cols, CV_32FC1, arr);
		Mat outMatrix;
		Mat im_color = Mat::zeros(rows, cols, CV_32FC4);

		vector<Mat> planes;
		vector<Mat>::iterator itMat;
		for(int i = 0; i < 4; i++)
			planes.push_back(matrix);
		merge(planes, im_color);

		for(int i = 0; i < n_pixels; i++) {
			Vec4f &v = im_color.at<Vec4f>(i/cols, i%cols);
			if(v[0] > 0) {
				v.val[0] = v[0] / (projection * projection) * 1050;
				v.val[1] = v[1] / (projection * projection) * 450;// / channelDiv;
				v.val[2] = v[2] / (projection * projection) * 150;// / channelDiv;
				v.val[3] = 255;
				//cout << v.val[2] << endl;
			}
		}

		std::ostringstream ostr;
		ostr << stop;
		string sStop = ostr.str();
		string fileName = "final_route_"+map+"_"+algName+"_"+sStop+"_"+heuristic+".png";
		imwrite(fileName, im_color);
	}




	void define_intervals(int stop, int &xIntervals, int &yIntervals) {

		//cout << "Stop: " << stop << "  -  Avg Biom: " << avg_biomasa << endl;
		pixels_necesarios = ceil(stop / avg_biomasa);

		//cout << "Pixels necesarios: " << pixels_necesarios << endl;

		intervals = ceil(sqrt(pixels_necesarios));

		//cout << "Intervals: " << intervals << endl;

		yIntervals = ceil(COLS / (double) intervals);
		xIntervals = ceil(ROWS / (double) intervals);

		//cout << "xIntervals: " << xIntervals << "  -  yIntervals: " << yIntervals << endl;
		//cout << "R: " << ROWS << "  C: " << COLS << endl;
		//exit(0);
	}

	map<float,Grid> define_grids(int rows, int cols, const int &xIntervals, const int &yIntervals, float** biomass, float** friction) {
			int xPosGrid, yPosGrid, id = 1, c = 0, cont = 0, contValid = 0;
			Grid** totalGrids = new Grid*[xIntervals];
			for (int i = 0; i< xIntervals; i++) {
				totalGrids[i] = new Grid[yIntervals];
			}

			Point2D tmp;
			for(int i = 0; i < rows; i++)
				for(int j = 0; j < cols; j++) {
					xPosGrid = floor(i / intervals);
					yPosGrid = floor(j / intervals);
					//FIXME: Change tmp
					tmp.x = i;
					tmp.y = j;
					totalGrids[xPosGrid][yPosGrid].noElements++;
					if (biomass[i][j] >= 0 && friction[i][j] > 0 && biomass[i][j] > friction[i][j]) {
						totalGrids[xPosGrid][yPosGrid].elements.push_back(tmp);
						totalGrids[xPosGrid][yPosGrid].sum += biomass[i][j] / friction[i][j];
						totalGrids[xPosGrid][yPosGrid].value = biomass[i][j];
						totalGrids[xPosGrid][yPosGrid].biomass += biomass[i][j];
						totalGrids[xPosGrid][yPosGrid].friction += friction[i][j];
						//cout << totalGrids[xPosGrid][yPosGrid].sum << endl;
						contValid++;

					}
					else {
						totalGrids[xPosGrid][yPosGrid].invalidCells = totalGrids[xPosGrid][yPosGrid].invalidCells + 1;
						totalGrids[xPosGrid][yPosGrid].value = -9999;
					}

					c++;
					//cout << intervals << endl;
					/*if(totalGrids[xPosGrid][yPosGrid].elements.size() + totalGrids[xPosGrid][yPosGrid].invalidCells == 974169) {
						cout << totalGrids[xPosGrid][yPosGrid].elements.size() + totalGrids[xPosGrid][yPosGrid].invalidCells << endl;
						cout << totalGrids[xPosGrid][yPosGrid].elements.size() << " - " << totalGrids[xPosGrid][yPosGrid].invalidCells << endl;
					}*/
					if (totalGrids[xPosGrid][yPosGrid].elements.size() + totalGrids[xPosGrid][yPosGrid].invalidCells == (intervals * intervals)
					    && totalGrids[xPosGrid][yPosGrid].biomass > 0 && totalGrids[xPosGrid][yPosGrid].friction > 0
						&& totalGrids[xPosGrid][yPosGrid].biomass > totalGrids[xPosGrid][yPosGrid].friction
						/*&& totalGrids[xPosGrid][yPosGrid].elements.size() > totalGrids[xPosGrid][yPosGrid].invalidCells*/) {
							//cout << totalGrids[xPosGrid][yPosGrid].elements.size() + totalGrids[xPosGrid][yPosGrid].invalidCells << endl;
							//cout << totalGrids[xPosGrid][yPosGrid].sum << endl;
							totalGrids[xPosGrid][yPosGrid].id = id;
							totalGrids[xPosGrid][yPosGrid].biomassAvg = totalGrids[xPosGrid][yPosGrid].biomass / totalGrids[xPosGrid][yPosGrid].elements.size();
							if(totalGrids[xPosGrid][yPosGrid].biomassAvg >= avg_biomasa)
								totGridsAvg++;
							id++;
							cont++;
							float gridSum = totalGrids[xPosGrid][yPosGrid].biomass / totalGrids[xPosGrid][yPosGrid].friction;
							gridsMap.insert(pair<float,Grid>(gridSum, totalGrids[xPosGrid][yPosGrid]));
					}
				}
			totValidGrids = cont;
			return gridsMap;
		}

	Point2D find_centroid(map<float,Grid> grids, float** biomass, float** friction) {
		map<float,Grid>::iterator it;
		float xMax = FLT_MIN, xMin = FLT_MAX, yMax = FLT_MIN, yMin = FLT_MAX;
		/*map<float,Grid>::iterator it2;
		for ( it = gridsMap.begin(); it != gridsMap.end(); ++it) {
			float xMax = FLT_MIN, xMin = FLT_MAX, yMax = FLT_MIN, yMin = FLT_MAX;
			cout << it->second.elements.size() + it->second.invalidCells << "\t Relation: " << it->first  << "\t Biomass: " << it->second.biomass << "\t Friction: " << it->second.friction << endl;
		}
		cout << "Finished. " << gridsMap.size() << endl;
		//exit(0);*/
		if (!grids.empty()) {
			it = (--grids.end());
		} else {
			flag = false;
		}

		Point2D centroid;

		if(flag){
			cout << "Relation: " << it->first << endl;
			for (int i = 0; i < it->second.elements.size(); i++) {
				if(it->second.elements.at(i).x > xMax)
					xMax = it->second.elements.at(i).x;

				if(it->second.elements.at(i).x < xMin)
					xMin = it->second.elements.at(i).x;

				if(it->second.elements.at(i).y > yMax)
					yMax = it->second.elements.at(i).y;

				if(it->second.elements.at(i).y < yMin)
					yMin = it->second.elements.at(i).y;
			}

			centroid.x = xMin + std::round((xMax - xMin) / 2) ;
			centroid.y = yMin + std::round((yMax - yMin) / 2);
			//cout << xMin << " - " << xMax << " - " << yMin << " - " << yMax << endl;
			this->xMax = xMax; this->xMin = xMin;
			this->yMax = yMax; this->yMin = yMin;
			if(biomass[centroid.x][centroid.y] < 0 || biomass[centroid.x][centroid.y] < friction[centroid.x][centroid.y]){
				cout << "Invalid Centroid: " << centroid.x << ", " << centroid.y << endl;
				bool found = true;
				set<cellVecinos> vecinos = vecinos2(centroid.x, centroid.y);
				set<cellVecinos> celdas;
				celdas.insert(cellVecinos(centroid.x, centroid.y, 0));
				set <cellVecinos> :: iterator itr;

					while(found) {
						for (itr = vecinos.begin(); itr != vecinos.end(); ++itr){
							//cout << (*itr).x << ", " << (*itr).y << endl;
							if(biomass[(*itr).x][(*itr).y] > 0 &&  biomass[(*itr).x][(*itr).y] > friction[(*itr).x][(*itr).y]){
								//cout << biomass[(*itr).x][(*itr).y] << " " << friction[(*itr).x][(*itr).y] << endl;
								centroid.x = (*itr).x;
								centroid.y = (*itr).y;
								found = false;
								break;
							}
						}
						vecinos = vecinos3(vecinos);
					}
			}
		}
		else{
			centroid.x = -1;
			centroid.y = -1;
		}
		return centroid;
	}



	set<cellVecinos>vecinos2(int origen_x, int origen_y){//busca vecinos iniciales nada mas
			set<cellVecinos>distancias;

			int dx[] = { -1, -1, 0, 1, 1, 1, 0,-1 };
			int dy[] = {  0,  1, 1, 1, 0, -1, -1,-1 };

			set<cellVecinos> st;
			Point2D tmp;

			// insert (0, 0) cell with 0 distance
			st.insert(cellVecinos(origen_x, origen_y, 0));

			cellVecinos k = *st.begin();
			st.erase(st.begin());
			//cout << "k.x = " << k.x << " k.y = " << k.y << " k.distance = " << k.distance <<endl;
			// looping through all neighbours
			for (int i = 0; i < 8; i++){
				int x = k.x + dx[i];
				int y = k.y + dy[i];
				//cout << "x = " << x << " y = " << y << endl;
				// if not inside boundry, ignore them
				if (!isInsideGrid(x, y)){
					//cout << "fuera del grid" << endl;
					continue;
				}

				distancias.insert(cellVecinos(x, y, 0));
				tmp.x = x; tmp.y=y;
				active_raster.push_back(tmp);
			}


			return distancias;
		}

		set<cellVecinos> vecinos3(set<cellVecinos> vecinos){//obtiene vecinos que aun no han sido considerados
			set<cellVecinos>distancias;
			//cout << "\t recibo estos:" << endl;
			set <cellVecinos> :: iterator itr2;
			for (itr2 = vecinos.begin(); itr2 != vecinos.end(); ++itr2){
				//cout << "\t " << (*itr2).x << ", " << (*itr2).y << endl;
			}
			int dx[] = { -1, -1, 0, 1, 1, 1, 0,-1 };
			int dy[] = {  0,  1, 1, 1, 0, -1, -1,-1 };

			Point2D tmp;
			vector<Point2D>::iterator it;
			set <cellVecinos> :: iterator itr;
			for (itr = vecinos.begin(); itr != vecinos.end(); ++itr){
				cellVecinos k = *itr;
				//cout << "\t -k" << k.x << ", "<< k.y << endl;
			// looping through all neighbours
				for (int i = 0; i < 8; i++){
					int x = k.x + dx[i];
					int y = k.y + dy[i];
					it = find_if(active_raster.begin(), active_raster.end(), Point2D(x, y));
					//cout << "\t x = " << x << " y = " << y << " active_raster= " << active_raster[x][y] << endl;
					// if not inside boundry, ignore them
					if (!isInsideGrid(x, y) || it != active_raster.end()){
						//cout << "\t fuera del grid" << endl;
						continue;
					}
					distancias.insert(cellVecinos(x, y, 0));
					tmp.x = x; tmp.y = y;
					active_raster.push_back(tmp);

				}
			}

			return distancias;
		}

		void matrix_to_tiff(float** output_raster, int rows, int cols, int count) {
			GDALDataset *poDstDS;
			//char **papszOptions = NULL;
			GDALDriver *poDriver;
			OGRSpatialReference oSRS;
			string fileName = "cost_distance_output_" + to_string(count) + ".tiff";
			string proyeccion = "EPSG:" + epsg;
			cout << proyeccion << endl;
			poDriver = GetGDALDriverManager()->GetDriverByName("Gtiff");
			poDstDS = poDriver->Create( fileName.c_str(), cols, rows, 1, GDT_Float32, NULL);
			poDstDS->SetGeoTransform(adfGeoTransform);
			oSRS.SetWellKnownGeogCS(proyeccion.c_str());

			GDALRasterBand *poBand;
			float *pBuf = new float[rows * cols], maxVal = 0;
			for(int i = 0; i < rows; i++) {
				for (int j = 0; j < cols; j++) {
					pBuf[i * cols + j] = output_raster[i][j];
					/*if(output_raster[i][j] <= 2000)
						pBuf[i * cols + j] = output_raster[i][j];
					else
						pBuf[i * cols + j] = -9999;*/
					if(output_raster[i][j] > maxVal)
						maxVal = output_raster[i][j];
				}
			}

			poBand = poDstDS->GetRasterBand(1);
			poDstDS->GetRasterBand(1)->SetNoDataValue(-9999);
			poBand->RasterIO( GF_Write, 0, 0, cols, rows,
			                  pBuf, cols, rows, GDT_Float32, 0, 0 );
			GDALClose( (GDALDatasetH) poDstDS );
			cout << fixed << "Max Val: " << maxVal << endl;
		}

		void check_npa(float** npa_matrix, float** &biomass_matrix) {
			for (int i = 0; i < ROWS; i++) {
				for (int j = 0; j < COLS; j++) {
					if (npa_matrix[i][j] > 0) {
						biomass_matrix[i][j] = biomass_matrix[i][j] / 80;
					}
				}
			}
		}

		void check_waterbodies(float** water_matrix, float** &biomass_matrix) {
			for (int i = 0; i < ROWS; i++) {
				for (int j = 0; j < COLS; j++) {
					if (water_matrix[i][j] > 0) {
						biomass_matrix[i][j] = biomass_matrix[i][j] / 80;
					}
				}
			}
		}

		void reproject_coords(string map_biomass) {
			string coords = "coords.txt", cmd = "gdaltransform " + map_biomass + " -t_srs EPSG:4326 -s_srs EPSG:" + epsg + " < " + coords;
			array<char, 128> buffer;
			string result;
			shared_ptr<FILE> pipe(popen(cmd.c_str(), "r"), pclose);
			if (!pipe) throw std::runtime_error("popen() failed!");
			while (!feof(pipe.get())) {
				if (fgets(buffer.data(), 128, pipe.get()) != nullptr)
					result += buffer.data();
			}
			cout << "Source = " << result << endl;
			//memset(&buffer[0], 0, sizeof(buffer));
			//pipe.reset();
			//exit(0);
		}

		/*void divide_biomass(float** &biomass_matrix) {
			for (int i = 0; i < ROWS; i++) {
				for (int j = 0; j < COLS; j++) {
					biomass_matrix[i][j] = biomass_matrix[i][j] / 40;
				}
			}
		}*/

/*	void matrix_to_tiff(float** output_raster, int rows, int cols) {
			setenv("PYTHONPATH",".",1);
			Py_Initialize();
			PyObject *pName, *pModule, *pDict, *pFunc;
			PyObject* pArgs = PyTuple_New(rows*cols + 2);
			PyTuple_SetItem(pArgs, 0, Py_BuildValue("i", rows));
			PyTuple_SetItem(pArgs, 1, Py_BuildValue("i", cols));
			int c = 2;
			for (int i = 0; i < rows; i++)
				for (int j = 0; j < cols; j++, c++)
					PyTuple_SetItem(pArgs, c, Py_BuildValue("f", output_raster[i][j]));
			pName = PyString_FromString((char*)"write_array");
			pModule = PyImport_Import(pName);
			pDict = PyModule_GetDict(pModule);
			pFunc = PyDict_GetItemString(pDict, (char*)"writeArray");
		   if (PyCallable_Check(pFunc)){
			   PyErr_Print();
			   PyObject_CallObject(pFunc, pArgs);
			   //cout << "Done" << endl;
			   //PyObject_CallFunctionObjArgs(pFunc, pRows, pCols, pArgs);
			   //PyErr_Print();
		   } else {
			   printf("Err\n");
			   PyErr_Print();
		   }
		   //cout << "Done" << endl;
		}  */
};
