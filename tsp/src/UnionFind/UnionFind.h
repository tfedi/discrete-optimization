#ifndef TP3_UNIONFIND_H
#define TP3_UNIONFIND_H
#include <vector>

using namespace std;
class UnionFind{
public:
    UnionFind(int n);
    int find(int x);
    void unionFind(int x, int y);
private:
    vector<int> height;
    vector<int> parent;
    const int INITIAL_HEIGHT = 1;
};
#endif //TP3_UNIONFIND_H
