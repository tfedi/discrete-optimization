#include "UnionFind.h"

 UnionFind::UnionFind(int n){
    for (int i = 0; i < n; i++){
        height.push_back(INITIAL_HEIGHT);
        parent.push_back(i);
    }
}

int UnionFind::find(int x){
    if(parent[x] != x){
        parent[x] = find(parent[x]);
    }
    return parent[x];
}

void UnionFind::unionFind(int x, int y){
    x = find(x);
    y = find(y);
    if(height[x] < height[y]){
        parent[x] = y;
    }
    else{
        parent[y] = x;
    }

    if(height[y] == height[x]) height[x]++;
}
