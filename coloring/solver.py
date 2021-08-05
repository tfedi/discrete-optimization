#!/usr/bin/python
# -*- coding: utf-8 -*-

from ortools.sat.python import cp_model
from igraph import *

def toGraph(n,m,E):
    graph = {k:[] for k in range(0,n)}
    for e in E:
        graph[e[0]].append(e[1])
        graph[e[1]].append(e[0])
    sorted_graph = sorted(graph.items(), key=lambda v: len(v[1]), reverse=True)
    return sorted_graph

def greedy(n,m,E):
    adj_list = toGraph(n,m,E)
    colors = [-1]*n
    for v in adj_list:
        color = 0
        j = 0
        while j < len(v[1]):
            if colors[v[1][j]] == color:
                color += 1
                j = 0
            else:
                j += 1
        colors[v[0]] = color
    return max(colors)+1, colors

def ColouringCPSat(n,m,E):
    upper_bound, greedy_solution = greedy(n,m,E)

    model = cp_model.CpModel()
    colors = [model.NewIntVar(0, upper_bound-1, 'c%i' %i) for i in range(n)]

    for e in E:
        model.Add(colors[e[0]] != colors[e[1]])   # c_i != c_j para todo (i,j) en E
    
    model.Add(colors[0] == 0)
    max_colors = [model.NewIntVar(0, upper_bound-1, 'c%i' %i) for i in range(n)]
    for i in range(1,n+1):
        model.AddMaxEquality(max_colors[i-1], colors[0:i])
        if(i < n): model.Add(colors[i] <= max_colors[i-1] + 1)
    
    max_color = max_colors[n-1]
    model.Minimize(max_color)

    '''g = Graph(n=n, edges=E)
    cliques = g.cliques(int(n/4), int(n/2))

    for clique in cliques:
        print(clique)
        colors_in_clique = [colors[i] for i in clique]
        model.AddAllDifferent(colors_in_clique)'''

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = 60.0
    status = solver.Solve(model)

    if status == cp_model.FEASIBLE or status == cp_model.OPTIMAL:
        res = []
        for color in colors:
            res.append(solver.Value(color))
        
        return solver.Value(max_color)+1, res
    else:
        return upper_bound, greedy_solution


def solve_it(input_data):
    # Modify this code to run your optimization algorithm

    # parse the input
    lines = input_data.split('\n')

    first_line = lines[0].split()
    node_count = int(first_line[0])
    edge_count = int(first_line[1])

    edges = []
    for i in range(1, edge_count + 1):
        line = lines[i]
        parts = line.split()
        edges.append((int(parts[0]), int(parts[1])))

    count,solution = ColouringCPSat(node_count, edge_count, edges)
    
    # prepare the solution in the specified output format
    output_data = str(count) + ' ' + str(0) + '\n'
    output_data += ' '.join(map(str, solution))

    return output_data


import sys

if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        file_location = sys.argv[1].strip()
        with open(file_location, 'r') as input_data_file:
            input_data = input_data_file.read()
        print(solve_it(input_data))
    else:
        print('This test requires an input file.  Please select one from the data directory. (i.e. python solver.py ./data/gc_4_1)')

