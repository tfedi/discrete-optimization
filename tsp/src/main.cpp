#include <iostream>
#include <fstream>
#include <cstring>
#include "Graph/Graph.h"

Graph leer(string file, bool verbose){
    int n, m;
    ifstream f(file.c_str());
    f >> n >> m;

    Graph G(n, verbose);

    int i, j; float c;
    for (int k = 0; k < m; ++k) {
        f >> i >> j >> c;
        G.addEdge(i,j,c);
    }

    return G;
}

int main(int argc, char *argv[]) {
    string metodo = "NN";
    int K=0,cant_top=0,castigo=0,porcentaje=0,cant_iter=0;
    bool verbose = false, tabu = false, tabuSol = false, showCurrentOptimal = false;
    string file = "";
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-m") == 0) metodo = argv[i+1];
        if (strcmp(argv[i], "-v") == 0) verbose = true;
        if (strcmp(argv[i], "-o") == 0) showCurrentOptimal = true;
        if (strcmp(argv[i], "-f") == 0) file = argv[i+1];
        if (strcmp(argv[i], "-t") == 0){
            tabu = true;
            K = stoi(argv[i+1]);
            castigo = stoi(argv[i+2]);
            cant_top = stoi(argv[i+3]);
            porcentaje = stoi(argv[i+4]);
        }
        if (strcmp(argv[i], "-ts") == 0){
            tabuSol = true;
            cant_top = stoi(argv[i+1]);
            porcentaje = stoi(argv[i+2]);
            cant_iter = stoi(argv[i+3]);
        }
    }

    Graph G = leer(file, verbose);

    pair<cycle, weight> circuito;
    if (strcmp(metodo.c_str(), "NN") == 0) circuito = G.nearestNeighbor();
    if (strcmp(metodo.c_str(), "AC") == 0) circuito = G.shortestEdge();
    if (strcmp(metodo.c_str(), "heurAG") == 0) circuito = G.heurAG();
    if (strcmp(metodo.c_str(), "random") == 0) circuito = G.randomCycle();

    if (tabu)
        circuito.second = G.tabu_search(circuito.first, circuito.second, K, porcentaje, castigo, cant_top, showCurrentOptimal);
    else if (tabuSol)
        circuito.second = G.tabu_search_basado_en_ciclos(circuito.first, circuito.second, porcentaje, cant_top, cant_iter, showCurrentOptimal);

    cout << circuito.second << " 0" << endl;
    Graph::printCycle(circuito.first);

    return 0;
}


