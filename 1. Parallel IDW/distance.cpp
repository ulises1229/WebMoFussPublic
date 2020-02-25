/*
 * costDistance.cpp
 *
 *  Created on: 04/10/2017
 *      Author: lanase
 */

#include "distance.h"


	// Utility method for comparing two cells
	static bool operator<(const costDistance::cell& a, const costDistance::cell& b){
		if (a.distance == b.distance){
			if (a.x != b.x)
				return (a.x < b.x);
			else
				return (a.y < b.y);
		}
		return (a.distance < b.distance);
	}

	bool costDistance::isInsideGrid(int i, int j){
		return (i >= 0 && i < ROW && j >= 0 && j < COL && cost_raster[i][j] >= 0);
	}


	set<costDistance::cell>costDistance::vecinos(int origen_x, int origen_y){
		set<cell>distancias;
		float** dis;
		float sqtwo = (float)sqrt(2);
		dis = new float*[this->ROW];
		for(int i = 0; i< ROW; ++i)
			dis[i] = new float[COL];
		// initializing distance array by INT_MAX
		for (int i = 0; i < ROW; i++)
			for (int j = 0; j < COL; j++){
				dis[i][j] = INT_MAX;
			}

		// direction arrays for simplification of getting
		// neighbour
		int dx[] = { -1, -1, 0, 1, 1, 1, 0,-1 };
		int dy[] = {  0,  1, 1, 1, 0, -1, -1,-1 };

		set<cell> st;

		// insert (0, 0) cell with 0 distance
		st.insert(cell(origen_x, origen_y, 0));

		// initialize distance of (0, 0) with its grid value
		dis[origen_x][origen_y] = cost_raster[origen_x][origen_y];

		// get the cell with minimum distance and delete
		// it from the set
		cell k = *st.begin();
		st.erase(st.begin());
		int x, y;
		// looping through all neighbours
		for (int i = 0; i < 8; i++){
			x = k.x + dx[i];
			y = k.y + dy[i];
			//cout << "x = " << x << " y = " << y << endl;
			// if not inside boundry, ignore them
			if (!isInsideGrid(x, y)){
				//cout << "no inside grid" << endl;
				continue;
			}

			if(i % 2 == 0){//par = no es diagonal
				//cout << "N" << endl;
				dis[x][y] = (dis[k.x][k.y] + cost_raster[x][y])/2;
			}
			else{
				//cout << "D" << endl;
				dis[x][y] = sqtwo*(dis[k.x][k.y] + cost_raster[x][y])/2;
			}
			distancias.insert(cell(x, y, dis[x][y]));
			active_raster[x][y] = true;
			if( dis[x][y] < output_raster[x][y] ) //Se valida que el borde sea menor al valor existente en output_raster
				output_raster[x][y] = dis[x][y];
		}

		output_raster[origen_x][origen_y] = 0;

		for(int m=0; m<this->ROW; m++)
			delete[] dis[m];

		st.clear();
		return distancias;
	}


	void costDistance::acumulados(set<cell> active_costs, int origen_x, int origen_y, double max, double min, float** grid, bool diagonals_cost_more){
			int dx[] = { -1, -1, 0, 1, 1,  1,  0, -1 };
			int dy[] = {  0,  1, 1, 1, 0, -1, -1, -1 };
			int cont = 1;
			float sqtwo = (float)sqrt(2);
			while(!active_costs.empty()){
				cell k = *active_costs.begin();//inicio con la distancia menor
				active_costs.erase(active_costs.begin());

				//cout << "Estoy en: " << k.x << " - " << k.y << " - " << k.distance << endl;
//				if(k.x >= xMin && k.x <= xMax && k.y >= yMin && k.y <= yMax){
					cont++;
					int x, y;
					for (int i = 0; i < 8; i++){
						x = k.x + dx[i];
						y = k.y + dy[i];
						//cout << x << ", " << y << " " << active_raster[x][y]<< " - " << biomass[x][y] << endl;
						set<float>distancias;
						if(isInsideGrid(x,y)) {
							if(i % 2 != 0){// si es movimiento diagonal
								//cout << "Diagonal" << " "<< x << ", " << y << " " << cost_raster[x][y] << endl;
								if((x != origen_x || y != origen_y) && active_raster[x][y] == 0){
									float dist=0;
									int cont_x = k.x; int cont_y = k.y;
									int cont_x2 = k.x; int cont_y2 = k.y;
									int caminos = 1;
									//set<float>distancias;
									for(int z = 0; z < 8; z=z+2){

										int mov_x = rutas[i][z];
										int mov_y = rutas[i][z+1];

										dist += (cost_raster[cont_x+mov_x][cont_y+mov_y] + cost_raster[cont_x2][cont_y2])/2;
										//cout << mov_x << " " << mov_y << " " << "cost_raster[cont_x+mov_x][cont_y+mov_y] = "<< cost_raster[cont_x+mov_x][cont_y+mov_y]
										//<< " cost_raster[cont_x2][cont_y2] = "<< cost_raster[cont_x2][cont_y2] << " /2 =  " << dist << endl;
										if(caminos == 2 || caminos == 4) {
											//cout << dist+output_raster[k.x][k.y] << endl;
											distancias.insert(dist+output_raster[k.x][k.y]);
											cont_x = cont_x2 = k.x;
											cont_y = cont_y2 = k.y;
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
										distancias.insert( sqtwo * (cost_raster[k.x][k.y]+cost_raster[x][y]) / 2 + output_raster[k.x][k.y] );
									else
										distancias.insert( ((cost_raster[k.x][k.y]+cost_raster[x][y]) / 2) + output_raster[k.x][k.y] );

									float minimo = *distancias.begin();
									//cout << minimo << endl;
									if(minimo < output_raster[x][y] && minimo >= 0){
										//cout << "Inserté " << x << ", " << y << " - " << minimo << endl;
											bool validate_min = validate_limit(max, min, grid[x][y]);
											if(validate_min)
											{
//												cout << "es tru" << endl;
												output_raster[x][y] = minimo;
												active_costs.insert(cell(x,y,minimo));
											}
									}
								}
							}
							else{
								//cout << "Normal" << " "<< x << ", " << y << " " << cost_raster[x][y] << endl;
								float dist = 0;
								if((x != origen_x || y != origen_y) && active_raster[x][y] == 0 && isInsideGrid(x,y)){
									dist = ((cost_raster[k.x][k.y] + cost_raster[x][y])/2) + output_raster[k.x][k.y];
									if(dist < output_raster[x][y] && dist >= 0){
										//cout << "Inserté " << x << ", " << y << " - " << dist << endl;
										bool validate_dist = validate_limit(max, min, grid[x][y]);
										if(validate_dist)
										{
//											cout << "es tru" << endl;
											output_raster[x][y] = dist;
											active_costs.insert(cell(x,y,dist));
										}
									}
								}

							}
						}
					}
					active_raster[k.x][k.y] = true;
	//			}
				}


			active_costs.clear();
		}

	void costDistance::inicio_cost_distance(float** grid, int srcX, int srcY, double projection, bool is_relative, double max, double min, bool flag, bool diagonals_cost_more) {
		//cout << COL << " " << ROW << endl;
		this->cost_raster = new float*[this->ROW];
		this->output_raster = new float*[this->ROW];
		this->active_raster = new bool*[this->ROW];
		for(int i = 0; i< ROW; ++i) {
			this->cost_raster[i] = new float[COL];
			this->output_raster[i] = new float[COL];
			this->active_raster[i] = new bool[COL];
		}
		for(int x = 0; x < ROW; x++){
			for(int y = 0; y < COL; y++){
				/*if(grid[x][y] == 999999)
					grid[x][y] = -9999;*/

				if (is_relative)
					this->cost_raster[x][y] = grid[x][y] * projection;
				else
				 	this->cost_raster[x][y] = grid[x][y];

				this->active_raster[x][y] = false;

				if (flag)
					this->output_raster[x][y] = grid[x][y];
				else
					this->output_raster[x][y] = INT_MAX;
			}
		}

		for(int x = 0; x < ROW; x++){
			for(int y = 0; y < COL; y++){
				if(y%200==0)
					printf("[%f]\n", grid[x][y]);
			}
		}

		// the position of a cell that you want to display its neighbours
		cout << srcX << " - " << srcY << endl;
		active_raster[srcX][srcY] = 1;

		//se obtienen los vecinos proximos al origen y sus distancias calculadas. ordenas de menor a mayor
		set<cell> distancias = vecinos(srcX,srcY);

		/*set <cell> :: iterator itr;
		for (itr = distancias.begin(); itr != distancias.end(); ++itr){
			cout << (*itr).x << ", " << (*itr).y << " - dist =  " << (*itr).distance << endl;
		}*/
		//exit(0);
		acumulados(distancias, srcX, srcY, max, min, grid, diagonals_cost_more);//metodo para calcular demas vecinos.
		for(int r=0; r<ROW; r++){
			for(int c=0; c<COL; c++){
				if(output_raster[r][c] == INT_MAX)
					output_raster[r][c] = -9999;

			}
		}

	}

	void costDistance::freeMem() {
		for(int m=0; m < ROW; m++){
			delete[] output_raster[m];
			delete[] cost_raster[m];
			delete[] active_raster[m];
		}
	}

	bool costDistance::validate_limit(double max, double min, float friction) {
		bool res = false;
		float limit = (float)(max - min)*0.75;
		if(friction<=limit)
			res = true;
		return res;
	}
