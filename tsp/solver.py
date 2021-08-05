#!/usr/bin/python
# -*- coding: utf-8 -*-

import math
from collections import namedtuple

import os
from subprocess import Popen, PIPE

Point = namedtuple("Point", ['x', 'y'])

def length(point1, point2):
    return math.sqrt((point1.x - point2.x)**2 + (point1.y - point2.y)**2)

def solve_it(input_data):
    # parse the input
    lines = input_data.split('\n')

    nodeCount = int(lines[0])

    points = []
    for i in range(1, nodeCount+1):
        line = lines[i]
        parts = line.split()
        points.append(Point(float(parts[0]), float(parts[1])))

     # Writes the inputData to a temporay file

    tmp_file_name = 'tmp.data'
    tmp_file_path = './src/cmake-build-debug/' + tmp_file_name
    tmp_file = open(tmp_file_path, 'w')
    tmp_file.write("{} {}\n".format(nodeCount, int(((nodeCount-1)*nodeCount)/2)))
    for i in range (0, nodeCount-1):
        for j in range (i+1, nodeCount):
            tmp_file.write("{} {} {}\n".format(i,j,length(points[i], points[j])))
    tmp_file.close()

    process = Popen(['./src/cmake-build-debug/tsp', '-f', tmp_file_path, '-m', 'NN', '-t', '20', '0', '10', '50'], stdout=PIPE, universal_newlines=True)
    (stdout, stderr) = process.communicate()

    # removes the temporay file
    os.remove(tmp_file_path)

    return stdout.strip()



import sys

if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        file_location = sys.argv[1].strip()
        with open(file_location, 'r') as input_data_file:
            input_data = input_data_file.read()
        print(solve_it(input_data))
    else:
        print('This test requires an input file.  Please select one from the data directory. (i.e. python solver.py ./data/tsp_51_1)')

