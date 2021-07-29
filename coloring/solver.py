#!/usr/bin/python
# -*- coding: utf-8 -*-

from ortools.sat.python import cp_model

def toGraph(n,m,E):
    graph = {k:[] for k in range(0,n)}
    for e in E:
        graph[e[0]].append(e[1])
        graph[e[1]].append(e[0])
    sorted_graph = sorted(graph.items(), key=lambda v: len(v[1]), reverse=True)
    return sorted_graph

def greedy(n,m,E):
    adj_list = toGraph(n,m,E)
    colours = [-1]*n
    for v in adj_list:
        color = 0
        j = 0
        while j < len(v[1]):
            if colours[v[1][j]] == color:
                color += 1
                j = 0
            else:
                j += 1
        colours[v[0]] = color
    return max(colours)+1, colours

def ColouringCPSat(n,m,E):
    upper_bound, greedy_solution = greedy(n,m,E)

    model = cp_model.CpModel()
    colours = [model.NewIntVar(0, upper_bound-1, 'c%i' %i) for i in range(n)]
    for e in E:
        model.Add(colours[e[0]] != colours[e[1]])   # c_i != c_j para todo (i,j) en E
    
    color_count = model.NewIntVar(0, upper_bound-1, 'color_count') # funcion objetivo minimax
    model.AddMaxEquality(color_count, colours)
    model.Minimize(color_count)

    model.AddDecisionStrategy(colours, cp_model.CHOOSE_FIRST, cp_model.SELECT_MIN_VALUE)

    solver = cp_model.CpSolver()
    solver.parameters.search_branching = cp_model.FIXED_SEARCH
    solver.parameters.max_time_in_seconds = 600.0
    status = solver.Solve(model)

    if status == cp_model.FEASIBLE or status == cp_model.OPTIMAL:
        res = []
        for color in colours:
            res.append(solver.Value(color))
        
        return solver.Value(color_count)+1, res
    else:
        return upper_bound, greedy_solution

def IntProgramming(n,m,E):
    model = cp_model.CpModel()
    w = [model.NewIntVar(0,1,'w%i'%i) for i in range(n)]
    x = [[model.NewIntVar(0,1,'x%i%i'%(i,j)) for j in range(n)] for i in range(n)]
    
    for i in range(n):
        model.Add(sum(x[i]) == 1)

    for j in range(n):
        for e in E:
            model.Add(x[e[0]][j] + x[e[1]][j] <= w[j])
        
        x_j = []
        for i in range(n):
            x_j.append(x[i][j])
        model.Add(w[j] <= sum(x_j))

        if (j < n-1):
            model.Add(w[j] >= w[j+1])

    model.Minimize(sum(w))

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = 60.0
    status = solver.Solve(model)
    
    if status == cp_model.FEASIBLE or status == cp_model.OPTIMAL:
        res = [0]*n
        for i in range(0,n):
            for j in range(0,n):
                if (solver.Value(x[i][j]) == 1):
                    res[i] = j

    return int(solver.ObjectiveValue()), res

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

    #count,solution = greedy(node_count, edge_count, edges)
    #count,solution = ColouringCPSat(node_count, edge_count, edges)
    count,solution = IntProgramming(node_count, edge_count, edges)

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

