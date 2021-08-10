#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <map>
#include <list>
#include <algorithm>
#include <thread>
#include <random>
using namespace std;

const double INF = 1e09;
const int NODO_INICIAL = 0;
const int CANT_ITER_TABU = 10;
const int MAX_ITER_TABU = 100;
const int TAM_VECINDARIO_RANDOM = 1000;
const int ITER_TABU = 50;
const int CANT_HILOS = 8;

struct Punto{
    double x,y;
};

struct Vecino{
    int index;
    double dist;
    Vecino(int i, double d): index(i), dist(d) {};
};

struct Parametros_Thread{
    int from, to;
    vector<int> &ciclo;
    vector<Punto> &v;
    Parametros_Thread(int from,
                      int to,
                      vector<int> &ciclo,
                      vector<Punto> &v) : from(from), to(to), ciclo(ciclo), v(v) {};
};

double d(Punto p1, Punto p2){
    return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
}

double computar_distancia_recorrido(int from, int to, vector<int> &ciclo, vector<Punto> &v){
    double res = 0;
    for (int i = from; i < to-1; ++i) {
        res += d(v[ciclo[i]], v[ciclo[i+1]]);
    }
    return res;
}

double computar_distancia_ciclo(int from, int to, vector<int> &ciclo, vector<Punto> &v){
    return computar_distancia_recorrido(from, to, ciclo, v) + d(v[ciclo[to-1]], v[ciclo[from]]);
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

double calcular_swap(vector<int> &ciclo, tuple<int,int> s, vector<Punto> &v, double dist){
    int i = get<0>(s), j = get<1>(s);
    dist -= d(v[ciclo[i]], v[ciclo[i+1]]);
    dist -= d(v[ciclo[j]], v[ciclo[j+1]]);
    dist += d(v[ciclo[i]], v[ciclo[j]]);
    dist += d(v[ciclo[i+1]], v[ciclo[j+1]]);
    return dist;
}

vector<int> swap(vector<int> &ciclo, tuple<int,int> s){
    vector<int> nuevo;
    int i = get<0>(s), j = get<1>(s);
    for (int l = 0; l <= i; ++l) {nuevo.push_back(ciclo[l]);}
    for (int l = j; l >= i+1; --l) {nuevo.push_back(ciclo[l]);}
    for (int l = j+1; l < ciclo.size(); ++l) {nuevo.push_back(ciclo[l]);}
    return nuevo;
}

void two_opt(int from, int to, vector<int> &ciclo, vector<Punto> &v){
    double dist = computar_distancia_recorrido(from, to, ciclo, v);
    double mejor_dist = calcular_swap(ciclo, make_tuple(from, from+2), v, dist);
    tuple<int,int> mejor_swap = make_tuple(from, from+2);
    for (int i = from; i < to-4; ++i) {
        for (int j = i+2; j < to-2; ++j) {
            double nueva_dist = calcular_swap(ciclo, make_tuple(i,j),v,dist);
            if (nueva_dist < mejor_dist){
                mejor_dist = nueva_dist;
                mejor_swap = make_tuple(i,j);
            }
        }
    }
    if (mejor_dist < dist){
        vector<int> nuevo = swap(ciclo, mejor_swap);
        for(int i = from; i < to; i++) ciclo[i] = nuevo[i];
    }
}

void calcular_distancia_swaps(Parametros_Thread params){
    two_opt(params.from, params.to, params.ciclo, params.v);
}

void restar_iter_tabu(map<tuple<int,int>, int> &tabu_list){
    list<tuple<int,int>> a_borrar;
    for (auto &p: tabu_list) {
        if (p.second-- <= 0) a_borrar.push_back(p.first);
    }
    for (auto p: a_borrar){ tabu_list.erase(p); }
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
    double mejor_dist=dist_total;
    vector<int> mejor_ciclo = ciclo_hamiltoniano;
    map<tuple<int,int>, int> tabu_list;

    vector<thread> threads;
    const int CANT_VECINOS_POR_HILO = (int)ciclo_hamiltoniano.size()/CANT_HILOS;

    std::random_device rd;
    std::mt19937 gen(rd());
    for(int t = 0; t < MAX_ITER_TABU; t++){
        // Intensificación
        for(int i = 0; i < ciclo_hamiltoniano.size(); i+=CANT_VECINOS_POR_HILO){
            std::thread th(calcular_distancia_swaps, Parametros_Thread(i, min((int)ciclo_hamiltoniano.size(), i+CANT_VECINOS_POR_HILO), ciclo_hamiltoniano, v));
            threads.push_back(move(th));
        }
        for (std::thread & th : threads)
        {
            // If thread Object is Joinable then Join that thread.
            if (th.joinable())
                th.join();
        }

        dist_total = computar_distancia_ciclo(0, (int)ciclo_hamiltoniano.size(), ciclo_hamiltoniano, v);

        // Si mejora el que teníamos, me lo guardo
        if (dist_total < mejor_dist){
            mejor_dist = dist_total;
            mejor_ciclo = ciclo_hamiltoniano;
        }

        // Diversificación
        for (int k = 0; k < ITER_TABU; ++k) {
            // Vecindario random
            vector<tuple<int,int>> vecindario_random(TAM_VECINDARIO_RANDOM);
            for (int i = 0; i < TAM_VECINDARIO_RANDOM; ++i) {
                std::uniform_int_distribution<> distrib_i(0, n-4);
                int nuevo_i = distrib_i(gen);
                std::uniform_int_distribution<> distrib_j(nuevo_i+2, n-2);
                int nuevo_j = distrib_j(gen);
                vecindario_random[i] = make_tuple(nuevo_i, nuevo_j);
            }
            double mejor_dist_local = calcular_swap(ciclo_hamiltoniano, vecindario_random[0], v, dist_total);
            tuple<int,int> mejor_swap = vecindario_random[0];
            for (auto &vecino: vecindario_random) {
                double nueva_dist = calcular_swap(ciclo_hamiltoniano, vecino, v, dist_total);
                if (nueva_dist < mejor_dist_local && (tabu_list.count(vecino) == 0 || nueva_dist < mejor_dist)){
                    mejor_dist_local = nueva_dist;
                    mejor_swap = vecino;
                }
            }
            ciclo_hamiltoniano = swap(ciclo_hamiltoniano, mejor_swap);
            dist_total = mejor_dist_local;

            if (dist_total < mejor_dist){
                mejor_dist = dist_total;
                mejor_ciclo = ciclo_hamiltoniano;
            }

            restar_iter_tabu(tabu_list);
            tabu_list[mejor_swap] = CANT_ITER_TABU;
        }
    }

    cout << mejor_dist << " 0" << endl;
    for (int k = 0; k < n; ++k) {
        cout << mejor_ciclo[k] << " ";
    }
    cout << endl;

    return 0;
}