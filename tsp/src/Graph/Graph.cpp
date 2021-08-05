#include <iostream>
#include "Graph.h"
#include <random>
#include <stack>
#include <list>
#include <map>
#include <algorithm>

const int INF = 10e9;
using namespace std;

Graph::Graph(int n, bool verbose) : adjMatrix(n, vector<float>(n,0)), verbose(verbose){}

int Graph::size() {
    return adjMatrix.size();
}

void Graph::addEdge(int i, int j, float w) {
    adjMatrix[i][j] = w;
    adjMatrix[j][i] = w;
}

int Graph::arg_min_not_visited(int v, vector<bool>& visited){
    int min = INF;
    int arg;
    for(int i = 0; i < this->size(); i++) {
        if (!visited[i] && adjMatrix[v][i] < min) {
            min = adjMatrix[v][i];
            arg = i;
        }
    }
    return arg;
}

// retornamos un ciclo hamiltoniano H y la suma de los pesos de las aristas.
pair<cycle, weight> Graph::nearestNeighbor(){ // O(n**2)
    cycle hamiltonianCycle;
    vector<bool> visited(this->size()); // vamos guardando los nodos ya visitados.

    int v = INITIAL_NODE; // Nodo donde empezamos a generar el cicrcuito Hamiltoniano
    weight weight = 0;

    hamiltonianCycle.push_back(v);
    visited[v] = true;

    while(hamiltonianCycle.size() < this->size()){
        int w = arg_min_not_visited(v, visited);
        hamiltonianCycle.push_back(w);
        visited[w] = true;
        weight += adjMatrix[v][w];
        v = w;
    }
    weight += adjMatrix[v][INITIAL_NODE];
    return make_pair(hamiltonianCycle, weight);
}

void Graph::printCycle(vector<int> c) {
    for(auto &v : c){
        cout << v << " ";
    }
    cout << endl;
}

bool Graph::checkCycle(cycle cycle, weight weight) {
    if (size() != cycle.size()) return false; // no es hamiltoniano

    set<int> visitedNodes;
    float sum = 0;
    for (int i = 0; i < cycle.size()-1; ++i) {
        if (visitedNodes.count(cycle[i]) > 0) return false; // hay un vertice repetido
        sum += adjMatrix[cycle[i]][cycle[i+1]];
        visitedNodes.insert(cycle[i]);
    }
    if (visitedNodes.count(cycle[size()-1]) > 0) return false;
    sum += adjMatrix[cycle[size()-1]][cycle[0]]; // agregamos el peso de la arista que une al ultimo nodo con el primero
    return sum == weight; // no hay vertices repetidos, pero la suma es diferente al peso estimado
}

// retorna un pair con los dos vertices que componen la arista con menor peso
// y que cumpla con que no forma ciclos y no genera un grado mayor igual a 3.
pair<int, int> Graph::arg_min_shortest_edge(UnionFind& V, vector<int>& grados){
    int min = INF;
    pair<int,int> e;
    for (int i = 0; i < this->size(); ++i) {
        for (int j = i+1; j < this->size(); ++j) {
            if ((grados[i] < 2) && (grados[j] < 2) && (V.find(i)!=V.find(j)) && (adjMatrix[i][j] < min)){
                e = make_pair(i,j);
                min = adjMatrix[i][j];
            }
        }
    }
    return e;
}

void Graph::DFS(matrix &E, vector<int>& pred, vector<int>& order)
{
    int i = 0;
    int next = 0;
    vector<bool> visited(this->size(), false);
    stack<int> stack;
    stack.push(i);
    visited[i] = true;

    while (!stack.empty())
    {
        i = stack.top();

        int j = 0;
        while (j < E.size() && (E[i][j] == 0 || visited[j])){j++;}

        if (j < E.size()){
            visited[j] = true;
            pred[j] = i;
            next += 1;
            order[j] = next;
            stack.push(j);
        }
        else{
            stack.pop();
        }
    }
}

// A utility function to find the vertex with
// minimum key value, from the set of vertices
// not yet included in MST
int Graph::minKey(vector<int> &key, vector<bool> &agmSet)
{
    // Initialize min value
    int min = INF, min_index;

    for (int v = 0; v < size(); v++)
        if (!agmSet[v] && key[v] < min)
            min = key[v], min_index = v;

    return min_index;
}


void Graph::AGM_prim(matrix &g, vector<int> &parent)
{
    // Valores 'Keys' usado para tomar la minima arista
    vector<int> key(size(), INF);

    // Para respresentar el conjunto de vertices incluidos en el AGM.
    vector<bool> agmSet(size());

    // Siempre incluimos primero al vertice 0 en el AGM.
    key[0] = 0;
    parent[0] = -1; // El primer nodo siempre es la raiz del AGM.

    // El AGM tendra n vertices.
    for (int count = 0; count < size() - 1; count++)
    {
        // Tomar el vertice que tenga minimo valor de 'key'.
        int u = minKey(key, agmSet);

        // Agregamos el vertice que elegimos al AGM
        agmSet[u] = true;

        // Entre todos los vetices del grafo:
        for (int v = 0; v < size(); v++)

            // Si:
            //      * u es adyacente a v (u!=v)
            //      * u no pertenece a el AGM
            //      * l(u,v) < key[v]
            if (g[u][v] && !agmSet[v] && g[u][v] < key[v])
                // actualizamos su arista y el valor de 'Key'
                parent[v] = u, key[v] = g[u][v];
    }

}

void Graph::AGM_a_mat(vector<int> &agm, matrix &E){
    for (int i = 1; i < size(); ++i) {
        E[i][agm[i]] = 1;
        E[agm[i]][i] = 1;
    }
}

// Dado un order, peso_total=0 y un vector vacio H
// retorna en H el circuito hamiltoniano determinado por order.
void Graph::crear_circuito_hamiltoniano(vector<int> &order, weight &peso_total, cycle &H){
    for (int i = 1; i < size(); ++i) {
        int j = 1;
        while (order[j] != i) {j++;}
        H[i] = j;
        peso_total += adjMatrix[H[i]][H[i-1]];
    }
    peso_total += adjMatrix[H[size()-1]][H[0]];
}

// Heuristica del arbol generador
pair<cycle, weight> Graph::heurAG(){
    // T <-- AGM(G)
    vector<int> T(size(), 0);
    AGM_prim(adjMatrix, T);

    // E <-- T como mat de adyacencia
    matrix E(size(), vector<float>(size(),0));
    AGM_a_mat(T, E);

    // recorrer E usando DFS
    vector<int> pred(this->size(), 0);
    vector<int> order(this->size(), 0);
    DFS(E, pred, order);

    // H <-- armar un circuito hamiltoniano siguiendo el orden dado por D
    vector<int> H(size(),0);
    weight peso = 0;
    crear_circuito_hamiltoniano(order, peso, H);

    return make_pair(H, peso);
}

// En este metodo vamos a crear el ciclo utilizando la estructura
// UnionFind para detectar ciclos.
pair<cycle, weight> Graph::shortestEdge(){
    matrix Xt(size(), vector<float>(size(),0));
    vector<int> grados(this->size(), 0);

    UnionFind V(size());

    int i = 0;
    while(i < size()-1){
        pair<int,int> e = arg_min_shortest_edge(V, grados);
        Xt[e.first][e.second] = 1;
        Xt[e.second][e.first] = 1;
        V.unionFind(e.first, e.second);
        grados[e.first]++;
        grados[e.second]++;
        i++;
    }
    // buscamos los vertices con exactamente grado 1
    pair<int,int> ultima_arista;
    bool primero = false;
    i = 0;
    while (i < size()){
        if(grados[i]==1){
            if(!primero){
                ultima_arista.first = i;
                primero = true;
            }else{
                ultima_arista.second = i;
            }
        }
        i++;
    }

    Xt[ultima_arista.first][ultima_arista.second] = 1;
    Xt[ultima_arista.second][ultima_arista.first] = 1;

    // Hacemos DFS con la matriz de adyacencia Xt
    vector<int> pred(this->size(), 0);
    vector<int> order(this->size(), 0);
    DFS(Xt, pred, order);

    cycle H(size(),0);
    weight peso = 0;
    crear_circuito_hamiltoniano(order, peso, H);

    return make_pair(H, peso);
}

weight Graph::calcular_weight_swap(cycle &S, vector<vector<int>> &tabu, int i, int j){
    vector<int> S_swapeado = S;
    swap(S_swapeado, i, j);
    weight w = 0;
    for (int k = 0; k < S_swapeado.size()-1; ++k) {
        w += adjMatrix[S_swapeado[k]][S_swapeado[k+1]];
    }
    w += tabu[j][i];
    return w;
}

weight Graph::calcular_weight_swap(cycle &S, int i, int j){
    vector<int> S_swapeado = S;
    swap(S_swapeado, i, j);
    weight w = 0;
    for (int k = 0; k < S_swapeado.size()-1; ++k) {
        w += adjMatrix[S_swapeado[k]][S_swapeado[k+1]];
    }
    return w;
}

void Graph::swap(cycle &S, int i, int j){
    int aux = S[i];
    S[i] = S[j];
    S[j] = aux;
}

bool Graph::aspirationFunction(int currentOptimal, int tabuValue) {
    return currentOptimal > tabuValue;
}

vector<int> Graph::randSubVecindad(int porcentaje){

    vector<int> vecindad;

    for (int i=0; i<size(); ++i) vecindad.push_back(i);

    std::shuffle ( vecindad.begin(), vecindad.end() , std::mt19937(std::random_device()()));

    vector<int> sub_vecindad(vecindad.begin(), vecindad.begin()+int(porcentaje*size()/100));

    return sub_vecindad;

}

pair<cycle,weight> Graph::randomCycle(){
    vector<int> ciclo = randSubVecindad(100);
    int peso = 0;
    for (int i = 0; i < ciclo.size()-1; ++i) {
        peso += adjMatrix[ciclo[i]][ciclo[i+1]];
    }
    peso += adjMatrix[ciclo[ciclo.size()-1]][ciclo[0]];
    return make_pair(ciclo, peso);
}

// dada una matriz vacia vecindad de size() x size(),
// el ciclo hayado por una funcion heuristica y una copia del peso del ciclo:
// modifica vecindad con los valores de cada cada swap en la posicion i, j.
weight Graph::tabu_search(cycle &ciclo, weight weight, int K = 20, int porcentaje = 20, int castigo = 1000, int cant_top = 10, bool showCurrentOptimal = false){
    vector<vector<int>> tabu(size(), vector<int>(size(),0));
    float currentOptimal = weight;
    float weightCastigado = 0;
    float currentOptimalCastigado = currentOptimal;
    cycle mejorCiclo = ciclo;
    list<pair<float,pair<int,int>>> candidatos; // top 5
    int t = 0;

    while(t < 1000){
        // Limpiamos el top 5 de la iteración anterior.
        candidatos.clear();

        // Restamos el turno a los swaps que ya hicimos y que están esperando para volver a ser utilizados.
        for (int i = 0; i < size(); ++i) {
            for (int j = i+1; j < size(); ++j) {
                if (tabu[i][j] > 0) tabu[i][j]--;
            }
        }

        // generamos una lista de tuplas (distintas) random para elegir una sub-vecindad.

        vector<int> sub_vecindad = randSubVecindad(porcentaje);

        // Calcular el top 5 de swaps posibles
        for (int i = 0; i < sub_vecindad.size()-1; ++i) {
            for (int j = i+1; j < sub_vecindad.size(); ++j) {
                int v1 = sub_vecindad[i];
                int v2 = sub_vecindad[j];
                float s = calcular_weight_swap(ciclo, tabu, v1, v2);
                candidatos.emplace_front(s,make_pair(v1,v2));
                candidatos.sort();
                if (candidatos.size() > cant_top) candidatos.pop_back();
            }
        }

        // Revisar si candidato esta en la lista tabu
        // En caso de no estarlo, se lo agrega a la lista tabu
        // y actualizamos el peso mejor actual (y si es mejor al mejor hasta el momento)
        // Si esta en la lista tabu, revisar si es > mejor hasta el momento
        for (auto const& candidato : candidatos) {
            if (tabu[candidato.second.first][candidato.second.second] == 0 || aspirationFunction(currentOptimal, candidato.first)){
                if (this->verbose){
                    cout << "IT: " << t << endl;
                    cout << "swapeando " << candidato.second.first << " " << candidato.second.second << endl;
                    cout << "peso viejo: " << weight << " nuevo peso: " << candidato.first << endl;
                }
                swap(ciclo, candidato.second.first, candidato.second.second);
                weight = candidato.first;
                weightCastigado = weight+tabu[candidato.second.second][candidato.second.first];
                tabu[candidato.second.first][candidato.second.second] = K+1;
                tabu[candidato.second.second][candidato.second.first] += castigo;
                break;
            }
        }
        if(currentOptimalCastigado > weightCastigado){
            currentOptimalCastigado = weightCastigado;
            currentOptimal = weight;
            mejorCiclo = ciclo;
        }
        t++;
        if (showCurrentOptimal) cout << weightCastigado << endl;
    }
    return currentOptimal;
}


weight Graph::tabu_search_basado_en_ciclos(cycle &ciclo, weight weight, int porcentaje, int cant_top, int cant_iter, bool showCurrentOptimal = false){
    map<cycle,int> tabuSoluciones;
    float currentOptimal = weight;
    cycle mejorCiclo = ciclo;
    list<pair<float,pair<int,int>>> candidatos; // top 5
    int t = 0;

    while(t < 1000){
        // Limpiamos el top 5 de la iteración anterior.
        candidatos.clear();

        list<map<cycle,int>::iterator> vencidos;
        for (auto& e:tabuSoluciones) {
            e.second--;
            if(e.second <= 0){
                vencidos.push_back(tabuSoluciones.find(e.first));
            }
        }
        for (auto it: vencidos){
            tabuSoluciones.erase(it);
        }

        // generamos una lista de tuplas (distintas) random para elegir una sub-vecindad.
        vector<int> sub_vecindad = randSubVecindad(porcentaje);

        // Calcular el top 5 de swaps posibles
        for (int i = 0; i < sub_vecindad.size()-1; ++i) {
            for (int j = i+1; j < sub_vecindad.size(); ++j) {
                int v1 = sub_vecindad[i];
                int v2 = sub_vecindad[j];
                float s = calcular_weight_swap(ciclo, v1, v2);
                candidatos.emplace_front(s,make_pair(v1,v2));
                candidatos.sort();
                if (candidatos.size() > cant_top) candidatos.pop_back();
            }
        }

        // Revisar si candidato esta en la lista tabu
        // En caso de no estarlo, se lo agrega a la lista tabu
        // y actualizamos el peso mejor actual (y si es mejor al mejor hasta el momento)
        // Si esta en la lista tabu, revisar si es > mejor hasta el momento
        for (auto const& candidato : candidatos) {
            swap(ciclo, candidato.second.first, candidato.second.second);

            if (tabuSoluciones.count(ciclo) == 0){
                if (this->verbose){
                    cout << "IT: " << t << endl;
                    cout << "swapeando " << candidato.second.first << " " << candidato.second.second << endl;
                    cout << "peso viejo: " << weight << " nuevo peso: " << candidato.first << endl;
                }

                weight = candidato.first;
                tabuSoluciones[ciclo] = cant_iter;
                break;
            }

            swap(ciclo, candidato.second.first, candidato.second.second);
        }
        if(currentOptimal > weight){
            currentOptimal = weight;
            mejorCiclo = ciclo;
        }
        t++;
        if (showCurrentOptimal) cout << weight << endl;
    }
    return currentOptimal;
}