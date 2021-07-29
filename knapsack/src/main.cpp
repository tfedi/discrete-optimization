#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cstring>

using namespace std;

struct item{
    int weight;
    int value;
    int order;
    bool operator<(item aux) const { return ((float)this->value/(float)this->weight  < (float)aux.value/(float)aux.weight);}
};

int n, K;
int** M;
vector<item> items;

const int UNDEFINED = -1;

int DP(int k, int j){
    if (j == 0) return 0;

    if (M[k][j] == UNDEFINED){
        if (items[j].weight <= k) M[k][j] = max(DP(k-items[j].weight, j-1) + items[j].value, DP(k, j-1));
        else M[k][j] = DP(k, j-1);
    }

    return M[k][j];
}

const int MINFTY = -10e6;
vector<bool> current_solution;
vector<bool> best_solution;
int best_value = 0;

float estimation(int k, int j){
    float res = 0;
    for (int i = j; i > 0; i--) {
        if (items[i].weight <= k){
            res += (float)items[i].value;
            k -= items[i].weight;
        }
        else{
            res += ((float)items[i].value/(float)items[i].weight)*(float)k;
            k = 0;
        }
    }
    return res;
}

int DFS(int k, int j, int v){
    if (k < 0) return MINFTY;
    if (j == 0){
        if (v > best_value){
            best_value = v;
            best_solution = current_solution;
        }
        return 0;
    }

    if ((float)v + estimation(k,j) <= (float)best_value) return MINFTY;

    current_solution[j] = true;
    int sol1 = DFS(k-items[j].weight, j-1, v+items[j].value) + items[j].value;
    current_solution[j] = false;
    int sol2 = DFS(k, j-1, v);

    return max(sol1,sol2);
}

int main(int argc, char* argv[]) {
    if (argc != 3) return 1;

    ifstream archivo(argv[2]);
    
    // Leer input
    archivo >> n >> K;
    items.resize(n+1);

    for (int i = 1; i <= n; ++i) {
        archivo >> items[i].value >> items[i].weight;
        items[i].order = i;
    }

    if (strcmp(argv[1], "DP") == 0){
        // Inicializar matriz de memoización
        M = new int*[K+1];
        for (int i = 0; i <= K; ++i) {
            M[i] = new int[n+1];
            for (int j = 0; j <= n; ++j) {
                M[i][j] = UNDEFINED;
            }
        }

        // Llamamos a la solución DP Bottom-Up
        for (int j = 0; j <= n; ++j) {
            for (int k = 0; k <= K; ++k) {
                DP(k,j);
            }
        }

        // Reconstruimos la solución
        int value = M[K][n];
        int* taken = new int[n+1];

        int i = K;
        int j = n;
        int leftValue = value;
        while (leftValue > 0){
            if (j > 0 and M[i][j] > M[i][j-1]){
                taken[j] = 1;
                leftValue = M[i][j] - items[j].value;
                while(i >= 0 and !(M[i][j-1] == leftValue and (i > 0 and M[i-1][j-1] < M[i][j-1])))
                    i -= 1;
            }
            j -= 1;
        }

        // Mostramos por stdout
        cout << value << " " <<  0 << endl;
        for (int k = 1; k <= n; ++k) {
            cout << taken[k] << " ";
        }
        cout << endl;
    }
    else if (strcmp(argv[1], "DFS") == 0){
        current_solution.resize(n+1, false);
        sort(items.begin()+1, items.end());
        int value = DFS(K, n, 0);

        bool* taken = new bool[n+1];
        for (int i = 0; i <= n; ++i) {
            taken[items[i].order] = best_solution[i];
        }

        cout << value << " " <<  0 << endl;
        for (int k = 1; k <= n; ++k) {
            cout << taken[k] << " ";
        }
        cout << endl;
    }


    
    return 0;
}
