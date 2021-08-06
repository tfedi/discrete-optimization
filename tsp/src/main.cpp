#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
using namespace std;

const double INF = 1e09;
const int NODO_INICIAL = 0;

struct Punto{
    double x,y;
};

struct Vecino{
    int index;
    double dist;
    Vecino(int i, double d): index(i), dist(d) {};
};

double d(Punto p1, Punto p2){
    return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
}

Vecino computar_vecindario(int i, vector<Punto> &v, vector<bool> &visitados, vector<double> &vecindario){
    Vecino mejor(i, INF);
    for (int j = 0; j < vecindario.size(); ++j) {
        vecindario[j] = visitados[j] ? INF : d(v[i], v[j]);
        if (vecindario[j] < mejor.dist){
            mejor.dist = vecindario[j];
            mejor.index = j;
        }
    }
    return mejor;
}

int main(int argc, char* argv[]) {
    // Levanto el archivo por parametro
    if (argc != 2) return 1;
    ifstream f(argv[1]);
    if (!f.is_open()) return 1;

    // Leo el archivo
    unsigned long n = 0;
    f >> n;
    vector<Punto> v(n);
    for (int i = 0; i < n; ++i) {
        f >> v[i].x >> v[i].y;
    }

    // Búsqueda local: heurística mejor vecino
    vector<bool> visitados(n,false);
    vector<double> vecindario(n, INF);
    vector<int> ciclo_hamiltoniano;
    int cant_visitados = 0, actual = NODO_INICIAL;
    double dist_total = 0;
    while (cant_visitados < n){
        visitados[actual]= true;
        ciclo_hamiltoniano.push_back(actual);
        if (cant_visitados < n-1){
            Vecino mejor_vecino = computar_vecindario(actual, v, visitados, vecindario);
            actual = mejor_vecino.index;
            dist_total += mejor_vecino.dist;
        }
        else {
            dist_total += d(v[ciclo_hamiltoniano[n-1]], v[NODO_INICIAL]);
        }

        cant_visitados++;
    }

    cout << dist_total << " 0" << endl;
    for (int i = 0; i < n; ++i) {
        cout << ciclo_hamiltoniano[i] << " ";
    }
    cout << endl;

    return 0;
}