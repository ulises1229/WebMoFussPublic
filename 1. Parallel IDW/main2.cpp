/*
 * main.cpp
 *
 *  Created on: 17/08/2017
 *      Author: lanase
 */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
//#include <float.h>
#include <vector>
#include <iterator>
#include <map>
#include "DisplayImage.cpp"
#include "distance.h"
//#include "gdal.cpp"
//#include "exploracion.h" //A*
//#include "tree.h" //Nodos
//#include "dijkstra.cpp" //Dijkstra
//#include "bellford.cpp" //Bellman-Ford
#include <tclap/CmdLine.h>
#include <omp.h>

struct points_export{
	int x, y, xMin, xMax, yMin, yMax;
};
char is_usable;

int main(int argc, const char** argv){
	string map_biomass, map_friction, hName, region, map_npa, map_waterbodies;
	char algorithm, heuristic;
	int demanda, optValidation, grids_to_validate, export_points;


  		// Define the command line object, and insert a message
  		// that describes the program. The "Command description message"
  		// is printed last in the help text. The second argument is the
  		// delimiter (usually space) and the last one is the version number.
  		// The CmdLine object parses the argv array based on the Arg objects
  		// that it contains.
  		TCLAP::CmdLine cmd("Command description message", ' ', "1");

  		// Define a value argument and add it to the command line.
  		// A value arg defines a flag and a type of value that it expects,
  		// such as "-n Bishop".
  		//TCLAP::ValueArg<std::string> nameArg("n","name","Name to print",true,"homer","string");
  		//TCLAP::ValueArg<std::string> biomass("b","biomass","Absolute path to biomass_map.tif",true,"/path/to/image.tif","string");
  		TCLAP::ValueArg<std::string> friction("f","friction","Absolute path to friction_map.tif",true,"/path/to/image.tif","string");
      TCLAP::ValueArg<std::string> friction("i","init","Absolute path to init_map.tif",true,"/path/to/image.tif","string");
/**
      TCLAP::ValueArg<std::string> algor("a","algorithm","Searching algorithm",true,"A","char");
  		TCLAP::ValueArg<std::string> stop("s","demand","Biomass demand",true,"5000","float");
  		TCLAP::ValueArg<std::string> watts("w","watts","Watts demand",true,"20","float");
  		TCLAP::ValueArg<std::string> humedad("i","humidity","Humidity content (%)",false,"40","float");
  		TCLAP::ValueArg<std::string> merma("l","loss","Percentage of expected loss in plant (%)",false,"30","float");
  		TCLAP::ValueArg<std::string> produccion("p","production","Annual percentage of production / operation of plant (%)",false,"20","float");
  		TCLAP::ValueArg<std::string> heur("u","heuristic","Searching heuristic",false,"d","char");
  		TCLAP::ValueArg<std::string> reg("r","region","Name of the region/country",false,"Haiti","string");
  		TCLAP::ValueArg<std::string> validation("v","validation","Validation option",true,"1","int");
  		TCLAP::ValueArg<std::string> grids_validate("g","grids_to_validate","Number of grids to validate",false,"50","int");
  		TCLAP::ValueArg<std::string> export_p("o","export_points","Export a number of points",false,"4","int");
  		TCLAP::ValueArg<std::string> npa("n","npa","Map of Natural Protected Areas (NPA)",false,"/path/to/npa_map.tif","int");
  		TCLAP::ValueArg<std::string> waterbody("m","waterbody","Map of Water bodies",false,"/path/to/water_bodies_map.tif","int");
  		TCLAP::ValueArg<std::string> usable("e","usable","Biomass' info represents usable biomass",true,"y/n","char");
  		// Add the argument nameArg to the CmdLine object. The CmdLine object
  		// uses this Arg to parse the command line.
*/
	try {
		cmd.xorAdd(stop, watts);
		cmd.add(biomass);
		cmd.add(friction);
		cmd.add(algor);
		cmd.add(heur);
		cmd.add(reg);
		cmd.add(validation);
		cmd.add(grids_validate);
		cmd.add(humedad);
		cmd.add(merma);
		cmd.add(produccion);
		cmd.add(export_p);
		cmd.add(npa);
		cmd.add(waterbody);
		cmd.add(usable);

		// Parse the argv array.
		cmd.parse( argc, argv );

		// Get the value parsed by each arg.
		map_biomass = biomass.getValue();
		map_friction = friction.getValue();
		string algoritmo = algor.getValue();
		algorithm = algoritmo[0];
		string heuristic2 = heur.getValue();
		heuristic = heuristic2[0];
		region = reg.getValue();
		string validacion = validation.getValue();
		optValidation = atoi(validacion.c_str());
		string grids_a_validar = grids_validate.getValue();
		grids_to_validate = atoi(grids_a_validar.c_str());
		string exp = export_p.getValue();
		export_points = atoi(exp.c_str());
		string demand;
		string usa = usable.getValue();
		is_usable = usa[0];

		if(optValidation > 4 || optValidation < 1) {
			cerr << "Please verify the validation option. (-v):\n 1 -- Best-candidate search.\n 2 -- Candidates where relation >= average validation.\n 3 -- Custom validation.\n 4 -- All-points validation." << endl;
			exit(0);
		}

		if(optValidation == 3 && !grids_validate.isSet()) {
			cerr << "Please indicate the number of grids to validate. (-g)." << endl;
			exit(0);
		}

		if(optValidation == 3 && grids_to_validate == 0) {
			cerr << "Please indicate a positive number of grids to validate. (-g)." << endl;
			exit(0);
		}

		if(optValidation == 2 || optValidation == 3){
			if(export_p.isSet() && export_points <= 0){
				cerr << "Please select a positive number for points to export. (-o)" << endl;
				exit(0);
			}else if(!export_p.isSet()){
				export_points = 1;
			}else if(export_p.isSet() && export_points > grids_to_validate){
				cerr << "Please select a number of points to export less than or equal to the number of grids to validate. (-o)" << endl;
				exit(0);
			}
		} else {
			export_points = 1;
		}


		if(!heur.isSet())
			heuristic = 'x';

		if(!grids_validate.isSet())
			grids_to_validate = 1;

		if(npa.isSet())
			map_npa = npa.getValue();

		if(waterbody.isSet())
			map_waterbodies = waterbody.getValue();

		if ( stop.isSet() ){
			demand = stop.getValue();
			demanda = strtof(demand.c_str(),0);
			if (demanda <= 0) {
				cerr << "Biomass demand must be greater than 0 (-s)." << endl;
				exit(0);
			}
		}
		else if ( watts.isSet() ){
			if(merma.isSet() && humedad.isSet() && produccion.isSet()){
				demand = watts.getValue();
				float w = strtof(demand.c_str(),0);
				if (w <= 0) {
					cerr << "Watts demand must be greater than 0 (-w)." << endl;
					exit(0);
				}
				string mer = merma.getValue();
				string hum = humedad.getValue();
				string prod = produccion.getValue();
				float humedad = strtof(hum.c_str(),0);
				float merm = strtof(mer.c_str(),0);
				float produ = strtof(prod.c_str(),0);
				float v1 = 1 - (humedad/100), v2 = 1 - (merm / 100), v3 = (produ/100);
				float hpa = 8760 * v3;
				float sp = w * hpa;
				float eb = sp / v2;
				float biomasa = eb / (5 * v1);
				demanda = biomasa;
			}else{
				cerr << "Please verify that the following flags are set: -i [--humidity] -l [--loss] -p [--production]" << endl;
				exit(0);
			}
		}

	} catch (TCLAP::ArgException &e)  // catch any exceptions
		{ std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }




	if(heuristic == 'e')
		hName = "Euclidean";
	else if(heuristic == 'm')
		hName = "Manhattan";
	else if(heuristic == 'd')
		hName = "Diagonal";
	else if(heuristic == 'a')
		hName = "More_than_median";
	else if(heuristic == 'b')
		hName = "Best_2_nodes";
	else
		hName = "No_heuristic";

	Display_image di;

	// Import biomass and friction
	clock_t begin = clock();
	float** biomass = di.tiff_to_matrix_gdal(map_biomass, true);
	float** friction = di.tiff_to_matrix_gdal(map_friction, false);
	cout << di.epsg << endl;
	//di.reproject_coords(map_biomass);
	//exit(0);
	clock_t end = clock();
	int rows = di.ROWS;
	int cols = di.COLS;

	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout << "TIFF to matrix = " << elapsed_secs << " secs." << endl;
	Point2D centroid;
	ofstream info, bestInfo, coords;
	ostringstream ss, bestSs;
	ss << "centroids_" << demanda << "_" << hName << ".csv";
	bestSs << "Exported_points.csv";
	//string fileName = "puntos_" + stop + ".csv"
	info.open(ss.str().c_str());
	info << "X, Y, Size, Biomass, Cost" << endl;
	bestInfo.open(bestSs.str().c_str());

	if(map_npa != "") {
		float** npa_matrix = di.tiff_to_matrix_gdal(map_npa, false);
		di.check_npa(npa_matrix, biomass);
	} else if(map_waterbodies != "") {
		float** water_matrix = di.tiff_to_matrix_gdal(map_waterbodies, false);
		di.check_waterbodies(water_matrix, biomass);
	}

	int xIntervals = 0, yIntervals = 0;
	di.define_intervals(demanda, xIntervals, yIntervals);

	if(optValidation == 4) {
		costDistance d;
		d.COL = cols;
		d.ROW = rows;

		float bestCost = 0, bestxMin, bestxMax, bestyMin, bestyMax;
		int bestX, bestY, cont = 1;

		int i=0, j=0;

            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    centroid.x = i; centroid.y = j;
                    if (biomass[i][j] > 0) {
                        cout << biomass[i][j] << endl;
                        coords.open("coords.txt");
                        cout << centroid.x << ", " << centroid.y << endl;
                        cout << "No. " << cont << " / " << di.valid_points << endl;
                        coords << centroid.y << " " << centroid.x;
                        coords.close();
                        di.reproject_coords(map_biomass);
                        clock_t begin_cd = clock();

                        d.inicio_cost_distance(friction, centroid.x, centroid.y, biomass, di.intervals, i - 80, i + 80, j - 80, j + 80, di.projection);


                        clock_t end_cd = clock();
                        double elapsed_secs_cd = double(end_cd - begin_cd) / CLOCKS_PER_SEC;

                        cout << "Cost Distance time = " << elapsed_secs_cd << " secs." << endl;
                        switch(algorithm) {
                            case 'B': //Binary search
                            case 'b': {
                                string algName = "Breadth_First_Search";
                                Tree rn;
                                rn.COL = cols;
                                rn.ROW = rows;
                                clock_t begin2 = clock();
                                rn.inicio_rutas(biomass, d.output_raster, centroid.x, centroid.y, demanda, info, heuristic);
                                clock_t end2 = clock();
                                double elapsed_secs2 = double(end2 - begin2) / CLOCKS_PER_SEC;
                                cout << "Nodes time = " << elapsed_secs2 << " secs." << "\n\n";

                                if(rn.cost > bestCost) {
                                    bestCost = rn.cost;
                                    bestX = rn.x;
                                    bestY = rn.y;
                                    bestxMin = i - 50;
                                    bestxMax = i + 50;
                                    bestyMin = j - 50;
                                    bestyMax = j + 50;
                                }

                                d.freeMem();
                                rn.freeMem();
                                break;
                            }

                            case 'A': //A* search
                            case 'a': {
                                Explore e;
                                string algName = "AStar";
                                e.COL = cols;
                                e.ROW = rows;
                                e.inicio(biomass);
                                clock_t begin2 = clock();
                                e.explore(d.output_raster, centroid.x, centroid.y, demanda, info, heuristic);
                                clock_t end2 = clock();
                                double elapsed_secs2 = double(end2 - begin2) / CLOCKS_PER_SEC;
                                cout << "A* Search time = " << elapsed_secs2 << " secs." << endl;

                                if(e.cost > bestCost) {
                                    bestCost = e.cost;
                                    bestX = e.X;
                                    bestY = e.Y;
                                    bestxMin = i - 50;
                                    bestxMax = i + 50;
                                    bestyMin = j - 50;
                                    bestyMax = j + 50;
                                }

                                e.freeMem();
                                d.freeMem();
                                break;

                            }
                        }
                        cont++;
                    }
                }
		}
		cout << "*** Best point ***" << endl;
		coords.open("coords.txt");
		coords << bestY << " " << bestX;
		coords.close();
		di.reproject_coords(map_biomass);
		//cout << "Source = " << bestX << ", " << bestY << endl;
		bestInfo << bestX << ", " << bestY << ",";
		if(algorithm == 'a' || algorithm == 'A') {
			d.inicio_cost_distance(friction, bestX, bestY, biomass, di.intervals, bestxMin, bestxMax, bestyMin, bestyMax, di.projection);
			Explore e;
			string algName = "AStar";
			e.COL = cols;
			e.ROW = rows;
			e.inicio(biomass);
			clock_t begin2 = clock();
			e.explore(d.output_raster, bestX, bestY, demanda, bestInfo, heuristic);
			clock_t end2 = clock();
			double elapsed_secs2 = double(end2 - begin2) / CLOCKS_PER_SEC;
			cout << "A* Search time = " << elapsed_secs2 << " secs." << endl;
			clock_t begin3 = clock();
			di.write_image(e.matrix_path, rows, cols, hName, demanda, region, algName);
			di.matrix_to_tiff(e.matrix_path, rows, cols, hName, demanda, region, algName);
			clock_t end3 = clock();
			double elapsed_secs3 = double(end3 - begin3) / CLOCKS_PER_SEC;
			cout << "Creating the final route image took " << elapsed_secs3 << " secs." << endl;
			e.freeMem();
		}
		else if(algorithm == 'b' || algorithm == 'B') {
			d.inicio_cost_distance(friction, bestX, bestY, biomass, di.intervals, bestxMin, bestxMax, bestyMin, bestyMax, di.projection);
			string algName = "Breadth_First_Search";
			Tree rn;
			rn.COL = cols;
			rn.ROW = rows;
			clock_t begin2 = clock();
			rn.inicio_rutas(biomass, d.output_raster, bestX, bestY, demanda, bestInfo, heuristic);
			clock_t end2 = clock();
			double elapsed_secs2 = double(end2 - begin2) / CLOCKS_PER_SEC;
			cout << "Nodes time = " << elapsed_secs2 << " secs." << "\n\n";

			clock_t begin3 = clock();
			cout << "Image creation..." << endl;
			di.write_image(rn.matrix_path, rows, cols, hName, demanda, region, algName);
			di.matrix_to_tiff(rn.matrix_path, rows, cols, hName, demanda, region, algName);
			clock_t end3 = clock();
			double elapsed_secs3 = double(end3 - begin3) / CLOCKS_PER_SEC;
			cout << "Creating the final route image took " << elapsed_secs3 << " secs." << endl;
			rn.freeMem();
		}
	}
	else {
		map<float,Grid> grids;
		grids = di.define_grids(rows, cols, xIntervals, yIntervals, biomass, friction);

		if(optValidation == 2)
			grids_to_validate = di.totGridsAvg;
		else if(optValidation == 1)
			grids_to_validate = 1;

		//bestInfo << "X, Y, Size, Biomass, Cost" << endl;
		string validac;
		if(optValidation == 1)
			validac = "Best-candidate search";
		else if(optValidation == 2)
			validac = "Candidates where relation >= average validation";
		else if(optValidation == 3)
			validac = "Custom validation";
		if(algorithm == 'a')
			bestInfo << "A*,";
		else if(algorithm == 't')
			bestInfo << "Depth Search";
		bestInfo << demanda << "," << hName << "," << region << "," << validac << "," << grids_to_validate << "," << export_points << endl;
		costDistance d;
		d.COL = cols;
		d.ROW = rows;

		int cont = 1;
		float bestCost = 0, bestxMin, bestxMax, bestyMin, bestyMax;
		int bestX, bestY;

		map <float, points_export> mapa_puntos;

		while(cont <= grids_to_validate) {
			coords.open("coords.txt");
			cout << "\n\n" << endl;
			centroid = di.find_centroid(grids, biomass, friction);
			if(!grids.empty())
				grids.erase(--grids.end());

			cout << "No. " << cont << " / " << grids_to_validate << endl;
			//centroid.x = 49; centroid.y = 93;
			//cout << "Source = " << centroid.x << ", " << centroid.y << endl;
			info << centroid.x << ", " << centroid.y << ",";
			coords << centroid.y << " " << centroid.x;
			coords.close();
			di.reproject_coords(map_biomass);

			clock_t begin_cd = clock();
			d.inicio_cost_distance(friction, centroid.x, centroid.y, biomass, di.intervals, di.xMin, di.xMax, di.yMin, di.yMax, di.projection);
			//cout << centroid.x << ", " << centroid.y << " - " << di.xMin << " - " << di.xMax << " - " << di.yMin << " - " << di.yMax << endl;
			//d.inicio_cost_distance(friction, centroid.x, centroid.y, biomass, di.intervals, 0, 498, 0, 256, di.projection);
			clock_t end_cd = clock();
			double elapsed_secs_cd = double(end_cd - begin_cd) / CLOCKS_PER_SEC;

			cout << "Cost Distance time = " << elapsed_secs_cd << " secs." << endl;

			//di.matrix_to_tiff(d.output_raster, rows, cols, "", 0, "", "cost_distance");

			//exit(0);

			switch(algorithm) {
				case 'B': //Binary search
				case 'b': {
					string algName = "Breadth_First_Search";
					Tree rn;
					rn.COL = cols;
					rn.ROW = rows;
					clock_t begin2 = clock();
					if(optValidation == 1) {
						bestInfo << centroid.x << ", " << centroid.y << ",";
						rn.inicio_rutas(biomass, d.output_raster, centroid.x, centroid.y, demanda, bestInfo, heuristic);
					}
					else
						rn.inicio_rutas(biomass, d.output_raster, centroid.x, centroid.y, demanda, info, heuristic);
					clock_t end2 = clock();
					double elapsed_secs2 = double(end2 - begin2) / CLOCKS_PER_SEC;
					cout << "Nodes time = " << elapsed_secs2 << " secs." << "\n\n";

					if (optValidation == 1) {
						clock_t begin3 = clock();
						cout << "Image creation..." << endl;
						di.write_image(rn.matrix_path, rows, cols, hName, demanda, region, algName);
						di.matrix_to_tiff(rn.matrix_path, rows, cols, hName, demanda, region, algName);
						clock_t end3 = clock();
						double elapsed_secs3 = double(end3 - begin3) / CLOCKS_PER_SEC;
						cout << "Creating the final route image took " << elapsed_secs3 << " secs." << endl;
						exit(0);
					}

					if(export_points == 1){
						if(rn.cost > bestCost) {
							bestCost = rn.cost;
							bestX = rn.x;
							bestY = rn.y;
							bestxMin = di.xMin;
							bestxMax = di.xMax;
							bestyMin = di.yMin;
							bestyMax = di.yMax;
						}
					}else{
						points_export pp; pp.x = centroid.x; pp.y = centroid.y; pp.xMin = di.xMin; pp.xMax = di.xMax; pp.yMin = di.yMin; pp.yMax = di.yMax;
						mapa_puntos.insert(make_pair(rn.cost, pp));
					}

					d.freeMem();
					rn.freeMem();
					break;
				}

				case 'A': //A* search
				case 'a': {
					Explore e;
					string algName = "AStar";
					e.COL = cols;
					e.ROW = rows;
					e.inicio(biomass);
					clock_t begin2 = clock();
					if(optValidation == 1) {
						bestInfo << centroid.x << ", " << centroid.y << ",";
						e.explore(d.output_raster, centroid.x, centroid.y, demanda, bestInfo, heuristic);
					}
					else
						e.explore(d.output_raster, centroid.x, centroid.y, demanda, info, heuristic);
					clock_t end2 = clock();
					double elapsed_secs2 = double(end2 - begin2) / CLOCKS_PER_SEC;
					cout << "A* Search time = " << elapsed_secs2 << " secs." << endl;

					if (optValidation == 1) {
						clock_t begin3 = clock();
						di.write_image(e.matrix_path, rows, cols, hName, demanda, region, algName);
						di.matrix_to_tiff(e.matrix_path, rows, cols, hName, demanda, region, algName);
						clock_t end3 = clock();
						double elapsed_secs3 = double(end3 - begin3) / CLOCKS_PER_SEC;
						cout << "Creating the final route image took " << elapsed_secs3 << " secs." << endl;
						exit(0);
					}

					if(export_points == 1){
						if(e.cost > bestCost) {
							bestCost = e.cost;
							bestX = e.X;
							bestY = e.Y;
							bestxMin = di.xMin;
							bestxMax = di.xMax;
							bestyMin = di.yMin;
							bestyMax = di.yMax;
						}
					}else{
						points_export pp; pp.x = centroid.x; pp.y = centroid.y; pp.xMin = di.xMin; pp.xMax = di.xMax; pp.yMin = di.yMin; pp.yMax = di.yMax;
						mapa_puntos.insert(make_pair(e.cost, pp));
					}

					e.freeMem();
					d.freeMem();
					break;
				}

				case 'D': //Dijkstra search
				case 'd': {
					Graph gr;
					gr.COL = cols;
					gr.ROW = rows;
					clock_t begin2 = clock();
					gr.dijkstra_inicio(biomass, d.output_raster, centroid.x, centroid.y, demanda, heuristic);
					clock_t end2 = clock();
					double elapsed_secs2 = double(end2 - begin2) / CLOCKS_PER_SEC;
					cout << "Dijkstra time = " << elapsed_secs2 << " secs." << endl;

					/*clock_t begin3 = clock();
					di.matrix_to_tiff(gr.matrix_path, rows, cols);
					clock_t end3 = clock();
					double elapsed_secs3 = double(end3 - begin3) / CLOCKS_PER_SEC;
					cout << "Matrix to TIFF = " << elapsed_secs3 << endl;*/
					break;
				}

				case 'F': //Bellman-Ford search
				case 'f': {
					Bellmanford bf;
					bf.COL = cols;
					bf.ROW = rows;
					clock_t begin2 = clock();
					bf.bellford_start(biomass, d.output_raster, centroid.x, centroid.y, demanda, heuristic);
					clock_t end2 = clock();
					double elapsed_secs2 = double(end2 - begin2) / CLOCKS_PER_SEC;
					cout << "Bellman-Ford time = " << elapsed_secs2 << " secs." << endl;

					/*clock_t begin3 = clock();
					di.matrix_to_tiff(bf.matrix_path, rows, cols);
					clock_t end3 = clock();
					double elapsed_secs3 = double(end3 - begin3) / CLOCKS_PER_SEC;
					cout << "Matrix to TIFF = " << elapsed_secs3 << endl;*/
					break;
				}
				default:
					exit(0);
					break;
			}
			cont++;
		}

		cout << endl;

		map <float, points_export> :: reverse_iterator itr;

		int cont_puntos = 1;
		if(export_points == 1){
			cout << "*** Best point ***" << endl;
			coords.open("coords.txt");
			coords << bestY << " " << bestX;
			coords.close();
			di.reproject_coords(map_biomass);
			//cout << "Source = " << bestX << ", " << bestY << endl;
			bestInfo << bestX << ", " << bestY << ",";
			if(algorithm == 'a' || algorithm == 'A') {
				d.inicio_cost_distance(friction, bestX, bestY, biomass, di.intervals, bestxMin, bestxMax, bestyMin, bestyMax, di.projection);
				Explore e;
				string algName = "AStar";
				e.COL = cols;
				e.ROW = rows;
				e.inicio(biomass);
				clock_t begin2 = clock();
				e.explore(d.output_raster, bestX, bestY, demanda, bestInfo, heuristic);
				clock_t end2 = clock();
				double elapsed_secs2 = double(end2 - begin2) / CLOCKS_PER_SEC;
				cout << "A* Search time = " << elapsed_secs2 << " secs." << endl;
				clock_t begin3 = clock();
				di.write_image(e.matrix_path, rows, cols, hName, demanda, region, algName);
				di.matrix_to_tiff(e.matrix_path, rows, cols, hName, demanda, region, algName);
				clock_t end3 = clock();
				double elapsed_secs3 = double(end3 - begin3) / CLOCKS_PER_SEC;
				cout << "Creating the final route image took " << elapsed_secs3 << " secs." << endl;
				e.freeMem();
			}
			else if(algorithm == 'b' || algorithm == 'B') {
				d.inicio_cost_distance(friction, bestX, bestY, biomass, di.intervals, bestxMin, bestxMax, bestyMin, bestyMax, di.projection);
				string algName = "Breadth_First_Search";
				Tree rn;
				rn.COL = cols;
				rn.ROW = rows;
				clock_t begin2 = clock();
				rn.inicio_rutas(biomass, d.output_raster, bestX, bestY, demanda, bestInfo, heuristic);
				clock_t end2 = clock();
				double elapsed_secs2 = double(end2 - begin2) / CLOCKS_PER_SEC;
				cout << "Nodes time = " << elapsed_secs2 << " secs." << "\n\n";

				clock_t begin3 = clock();
				cout << "Image creation..." << endl;
				di.write_image(rn.matrix_path, rows, cols, hName, demanda, region, algName);
				di.matrix_to_tiff(rn.matrix_path, rows, cols, hName, demanda, region, algName);
				clock_t end3 = clock();
				double elapsed_secs3 = double(end3 - begin3) / CLOCKS_PER_SEC;
				cout << "Creating the final route image took " << elapsed_secs3 << " secs." << endl;
				rn.freeMem();
			}

		}else{
			for(itr = mapa_puntos.rbegin(); itr != mapa_puntos.rend(); itr++){
				if(cont_puntos <= export_points){
					cout << "Best points *** " << cont_puntos << "***" << endl;
					coords.open("coords.txt");
					coords << itr->second.y << " " << itr->second.x;
					coords.close();
					di.reproject_coords(map_biomass);
					//cout << "Source = " << itr->second.x << ", " << itr->second.y << endl;
					bestInfo << itr->second.x << ", " << itr->second.y << ",";
					if(algorithm == 'a' || algorithm == 'A') {
						d.inicio_cost_distance(friction, itr->second.x, itr->second.y, biomass, di.intervals, itr->second.xMin, itr->second.xMax, itr->second.yMin, itr->second.yMax, di.projection);
						Explore e;
						stringstream s;
						s << cont_puntos << "_A_Star";
						string algName = s.str();
						e.COL = cols;
						e.ROW = rows;
						e.inicio(biomass);
						clock_t begin2 = clock();
						e.explore(d.output_raster, itr->second.x, itr->second.y, demanda, bestInfo, heuristic);
						clock_t end2 = clock();
						double elapsed_secs2 = double(end2 - begin2) / CLOCKS_PER_SEC;
						cout << "A* Search time = " << elapsed_secs2 << " secs." << endl;
						clock_t begin3 = clock();
						di.write_image(e.matrix_path, rows, cols, hName, demanda, region, algName);
						di.matrix_to_tiff(e.matrix_path, rows, cols, hName, demanda, region, algName);
						clock_t end3 = clock();
						double elapsed_secs3 = double(end3 - begin3) / CLOCKS_PER_SEC;
						cout << "Creating the final route image took " << elapsed_secs3 << " secs." << endl;
						e.freeMem();
					}
					else if(algorithm == 'b' || algorithm == 'B') {
						d.inicio_cost_distance(friction, itr->second.x, itr->second.y, biomass, di.intervals, itr->second.xMin, itr->second.xMax, itr->second.yMin, itr->second.yMax, di.projection);
						stringstream s;
						s << cont_puntos << "_Breadth_First_Search";					string algName = s.str();
						Tree rn;
						rn.COL = cols;
						rn.ROW = rows;
						clock_t begin2 = clock();
						rn.inicio_rutas(biomass, d.output_raster, itr->second.x, itr->second.y, demanda, bestInfo, heuristic);
						clock_t end2 = clock();
                    double elapsed_secs2 = double(end2 - begin2) / CLOCKS_PER_SEC;
						cout << "Nodes time = " << elapsed_secs2 << " secs." << "\n\n";

						clock_t begin3 = clock();
						cout << "Image creation..." << endl;
						di.write_image(rn.matrix_path, rows, cols, hName, demanda, region, algName);
						di.matrix_to_tiff(rn.matrix_path, rows, cols, hName, demanda, region, algName);
						clock_t end3 = clock();
						double elapsed_secs3 = double(end3 - begin3) / CLOCKS_PER_SEC;
						cout << "Creating the final route image took " << elapsed_secs3 << " secs." << endl;
						rn.freeMem();
					}

					cont_puntos++;
				}else{
					break;
				}
			}
		}
	}
	info.close();
	bestInfo.close();
}
