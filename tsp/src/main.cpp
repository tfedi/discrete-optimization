#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <map>
#include <list>
#include <algorithm>
#include <thread>
using namespace std;

const double INF = 1e09;
const int NODO_INICIAL = 0;
const int CANT_ITER_TABU = 65;
const int MAX_ITER_TABU = 10000;
const int CANT_VECINOS_POR_HILO = 500;
struct Punto{
    double x,y;
};

struct Vecino{
    int index;
    double dist;
    Vecino(int i, double d): index(i), dist(d) {};
};

struct Parametros_Thread{
    vector<double> &distancias;
    vector<pair<int,int>> &vecindad;
    int from, to;
    vector<int> &ciclo;
    vector<Punto> &v;
    double dist;
    Parametros_Thread(vector<double> &distancias,
                      vector<pair<int,int>> &vecindad,
                      int from,
                      int to,
                      vector<int> &ciclo,
                      vector<Punto> &v,
                      double dist) : distancias(distancias), vecindad(vecindad), from(from), to(to), ciclo(ciclo), v(v), dist(dist) {};
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


double calcular_swap(vector<int> &ciclo, int i, int j, vector<Punto> &v, double dist){
    if (i == j) return dist;

    int prev_i = i == 0 ? ciclo.size()-1 : i-1;
    int sig_i = i+1;
    int prev_j = j-1;
    int sig_j = j == ciclo.size()-1 ? 0 : j+1;

    i = ciclo[i]; j = ciclo[j]; prev_i = ciclo[prev_i], sig_i = ciclo[sig_i], prev_j = ciclo[prev_j], sig_j = ciclo[sig_j];
    if (prev_i == j && sig_j == i){
        dist -= d(v[i], v[sig_i]);
        dist -= d(v[j], v[prev_j]);
        dist += d(v[j], v[sig_i]);
        dist += d(v[i], v[prev_j]);
    }
    else if (sig_i == j && prev_j == i){
        dist -= d(v[i], v[prev_i]);
        dist -= d(v[j], v[sig_j]);
        dist += d(v[j], v[prev_i]);
        dist += d(v[i], v[sig_j]);
    }
    else{
        dist -= d(v[i], v[prev_i]);
        dist -= d(v[i], v[sig_i]);
        dist -= d(v[j], v[prev_j]);
        dist -= d(v[j], v[sig_j]);
        dist += d(v[i], v[prev_j]);
        dist += d(v[i], v[sig_j]);
        dist += d(v[j], v[prev_i]);
        dist += d(v[j], v[sig_i]);
    }
    return dist;
}

void restar_iter_tabu(map<pair<int,int>, int> &tabu_list){
    list<pair<int,int>> a_borrar;
    for (auto &p: tabu_list) {
        if (p.second-- <= 0) a_borrar.push_back(p.first);
    }
    for (auto p: a_borrar){ tabu_list.erase(p); }
}

void calcular_distancia_swaps(Parametros_Thread params){
    for (int i = params.from; i < min((int)params.distancias.size(), params.to); ++i) {
        params.distancias[i] = calcular_swap(params.ciclo, params.vecindad[i].first, params.vecindad[i].second, params.v, params.dist);
    }
    //cout << "Calculadas las distancias desde " << params.from << " hasta " << min((int)params.distancias.size(), params.to) << endl;
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

    // Obtenemos un ciclo inicial

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

    // Búsqueda tabú
    int mejor_i, mejor_j, t = 0;
    double mejor_dist=dist_total, dist_swap, mejor_dist_vecina;
    vector<int> mejor_ciclo = ciclo_hamiltoniano;
    map<pair<int,int>, int> tabu_list;

    vector<pair<int,int>> vecindad;
    for (int i = 0; i < n-1; ++i) {
        for (int j = i+1; j < n; ++j) {
            vecindad.emplace_back(make_pair(i,j));
        }
    }

    vector<double> dist_por_swap(vecindad.size());
    vector<thread> threads;

    while (t < MAX_ITER_TABU){
        // Calculo nuevo vecindario
        for(int i = 0; i < vecindad.size(); i+=CANT_VECINOS_POR_HILO){
            std::thread th(calcular_distancia_swaps, Parametros_Thread(dist_por_swap, vecindad, i, i+CANT_VECINOS_POR_HILO, ciclo_hamiltoniano, v, dist_total));
            threads.push_back(move(th));
        }
        for (std::thread & th : threads)
        {
            // If thread Object is Joinable then Join that thread.
            if (th.joinable())
                th.join();
        }

        // Busco el mejor swap disponible
        mejor_i = 0;
        mejor_j = 0;
        mejor_dist_vecina = INF;
        for(int i = 0; i < dist_por_swap.size(); i++){
            if (dist_por_swap[i] < mejor_dist_vecina && (tabu_list.count(make_pair(vecindad[i].first, vecindad[i].second)) == 0 || dist_swap < mejor_dist)){
                mejor_dist_vecina = dist_por_swap[i];
                mejor_i = vecindad[i].first;
                mejor_j = vecindad[i].second;
            }
        }

        // Hago el mejor swap disponible
        swap(ciclo_hamiltoniano[mejor_i], ciclo_hamiltoniano[mejor_j]);
        dist_total = mejor_dist_vecina;

        // Si mejora el que teníamos, me lo guardo
        if (dist_total < mejor_dist){
            mejor_dist = dist_total;
            mejor_ciclo = ciclo_hamiltoniano;
        }

        // Actualizo lista tabú
        restar_iter_tabu(tabu_list);
        tabu_list[make_pair(mejor_i,mejor_j)] = CANT_ITER_TABU;
        t++;
    }

    cout << mejor_dist << " 0" << endl;
    for (int k = 0; k < n; ++k) {
        cout << mejor_ciclo[k] << " ";
    }
    cout << endl;

    return 0;
}