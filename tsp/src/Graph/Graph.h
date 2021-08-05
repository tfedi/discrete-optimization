#ifndef TP3_GRAPH_H
#define TP3_GRAPH_H
#include <vector>
#include <set>
#include "../UnionFind/UnionFind.h"

using namespace std;
typedef std::vector<int> cycle;
typedef std::vector<vector<float>> matrix;
typedef float weight;

class Graph {
public:
    explicit Graph(int n, bool verbose);
    void addEdge(int i, int j, float w);
    int size();
    pair<cycle, weight> nearestNeighbor();
    pair<cycle, weight> shortestEdge();
    pair<cycle, weight> heurAG();
    pair<cycle,weight> randomCycle();
    static void printCycle(cycle c);
    void DFS(matrix &E, vector<int> &pred, vector<int> &order);
    int minKey(vector<int> &key, vector<bool> &agmSet);
    void AGM_prim(matrix &g, vector<int> &parent);
    void AGM_a_mat(vector<int> &agm, matrix &E);
    bool checkCycle(cycle cycle, weight weight);
    void crear_circuito_hamiltoniano(vector<int> &order, weight &peso_total, cycle &H);
    weight calcular_weight_swap(cycle &S, vector<vector<int>> &tabu, int i, int j);
    weight calcular_weight_swap(cycle &S, int i, int j);
    void swap(cycle &S, int i, int j);
    weight tabu_search(cycle &ciclo, weight weight, int K, int porcentaje, int castigo, int cant_top, bool showCurrentOptimal);
    weight tabu_search_basado_en_ciclos(cycle &ciclo, weight weight, int porcentaje, int cant_top, int cant_iter, bool showCurrentOptimal);

private:
    matrix adjMatrix;
    bool verbose;
    const int INITIAL_NODE = 0;
    bool aspirationFunction(int currentOptimal, int tabuValue);
    int arg_min_not_visited(int v, vector<bool>& visited);
    pair<int,int> arg_min_shortest_edge(UnionFind& V, vector<int>& grados);
    vector<int> randSubVecindad(int porcentaje);

    bool cicloEsTabu(cycle &ciclo, matrix &tabu);
};


#endif //TP3_GRAPH_H




