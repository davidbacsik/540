'''Script to find the highest weight path in a graph.

To see usage:
    python process_graph.py --help

by David Bacsik'''

import os
import time
import argparse
import glob

def parseArgs():
    '''Takes in command line flags and passes them to variables.'''

    parser = argparse.ArgumentParser(
            description = 'Find the highest weight path in a graph. See '
            'template.csv for graph file format.')
    parser.add_argument('--graph', help='Graph file, in CSV format',
            required=True)
    parser.add_argument('--out', help='Prefix for output file (optional)',
            required=False, default='')

    return parser

def readGraph(graph_file):
    '''
    Reads graph, parses vertices and edges, and passes them each to a list.

    Args:
        graph_file (str): Name of file containing graph in CSV format.

    Returns:
        vertices (list): List of all vertex names; names must be unique
        edges (list): List of all edges. Format for each edge is:
                        edge_name, start_vertex, end_vertex, weight
        start_vertex (str): If supplied, name of contsraining start vertex
        end_vertex (str): If supplied, name of constraining end vertex
    '''

    # Initilize vertex and edge lists
    vertices = []
    edges = []
    start_vertex = None
    end_vertex = None

    # Open graph file
    with open(graph_file, 'r') as graph:
        # Iterate through each item in graph file
        for line in graph:
            line = line.strip('\n')
            line_split = line.split(',')

            # Vertices
            if line_split[0] == 'V':
                current_vertex = line_split[1]
                assert(len(current_vertex) <= 10),"Vertex name must be 10 characters or less."
                assert (current_vertex not in vertices),"Duplicate vertex name."
                vertices.append(current_vertex)

                # Parse START/END flags
                if 'START' in line_split:
                    start_vertex = line_split[1]
                elif 'END' in line_split:
                    end_vertex = line_split[1]

            # Edges
            elif line_split[0] == 'E':
                current_edge = [line_split[1], line_split[2], line_split[3], \
                    float(line_split[4])]
                edges.append(current_edge)

            # Ignore other lines
            else:
                pass

    return vertices, edges, start_vertex, end_vertex

def writeOutput(out_prefix, run_time, path_score, first_vertex, last_vertex, edge_path):
    '''
    Writes output text file.

    Args:
        out_prefix (str): Name of output file, from command line
        run_time (float): Total number of seconds the script took to run
        path_score (float): The total weight of the path; also the weight
                of the final vertex
        path (list): Ordered list of vertices, from end of path to beginning
        edge_path (list): Ordered list of edges, from beginning of path to end
        start_vertex (str): Name of the first vertex in the path
        end_vertex (str): Name of the last vertex in the path

    Returns:
        out_prefix_output.txt (txt): Text file containing output

    '''

    # Header
    if out_prefix != '':
        out_prefix = out_prefix+'_'
    with open ('{0}output.txt'.format(out_prefix), 'w') as out_file:
        # Header
        out_file.write('Assignment: GS 540 HW2\n')
        out_file.write('Name: David Bacsik\n')
        out_file.write('Email: dbacsik@uw.edu\n')
        out_file.write('Language: Python\n')
        out_file.write('Runtime: '+str(run_time)+' sec\n\n')

        # Path Summary
        out_file.write('Part X\n')
        out_file.write('Score: '+str(path_score)+'\n')
        out_file.write('Begin: '+str(first_vertex)+'\n')
        out_file.write('End: '+str(last_vertex)+'\n')
        out_file.write('Path: '+''.join(map(str, edge_path)))

def calculateWeights(vertices, edges, start_index, start_vertex):
    '''
    Calculates the weight (the maximum weight) of a path ending at all vertices
        in a graph, starting from designated vertex, if specified.

    Args:
        vertices (list): List of vertices in the graph
        edges (list of lists): List of all edges in the graph.
                Indices/Format:
                (0) edge_name, (1) start_vertex, (2) end_vertex, (3) weight
        start_index (int): Index of vertex to start from
        start_vertex (str): Name of the first vertex in the path

    Returns:
        weights (dict): Dictionary of vertex (max) weights, keyed by vertex name
        parent_dict (dict): Dictionary of parent giving max path weight to each vertex.
                Keyed by vertex, value is parent.
    '''

    weights = {}
    parent_dict = {}

    # Iterate through vertices
    for i in range(start_index,len(vertices)):

        # Initialize list of all weights for paths ending at current vertex
        current_vertex_weights = {vertices[i]: 0}

        # Iterate through all edges ending at current vertex
        # Set start index weight to 0
        if i == start_index:
            pass
        else:
            for e in edges:
                if e[2] == vertices[i] and weights.get(e[1]) != None:
                    parent_weight = weights[e[1]]
                    edge_weight = e[3]
                    total_weight = parent_weight + edge_weight
                    current_vertex_weights[e[1]] = total_weight

                    # If start vertex constrained, and a vertex is descended from constrained start vertex, must pick that as parent
                    if start_vertex != None and e[1] == vertices[start_index]:
                        current_vertex_weights = {e[1]: total_weight}

            # If there is a constrained start vertex, and the current vertex is not connected to it remove self-path as an option
            if start_vertex != None and len(current_vertex_weights) > 1 and current_vertex_weights.get(vertices[i]) != None:
                del current_vertex_weights[vertices[i]]

        # Pick the best graph of available option
        weight = max(current_vertex_weights.values())
        best_parent = max(current_vertex_weights.iterkeys(), key=lambda k: current_vertex_weights[k])
        parent_dict[vertices[i]] = best_parent
        weights[vertices[i]] = weight

    return weights, parent_dict

def highestParent(vertex, parent_dict):
    '''
    For a given vertex, returns the parent with the highest weight.

    Args:
        vertex (str): Name of vertex to be tested
        parent_dict (dict): Dictionary of parent giving max path weight to each vertex.
                Keyed by vertex, value is parent.

    Returns:
        parent (str): Name of parent with the highest weight
    '''

    parent = parent_dict[vertex]

    return parent

def trace(vertex, vertices, edges, weights, parent_dict, start_vertex):
    '''
    Traces the highest weight path from a given vertex until it reaches a
        vertex with weight 0.

    Args:
        vertex (str): Name of vertex to start traceback
        vertices (list): List of vertices in the graph
        edges (list of lists): List of all edges in the graph.
                Indices/Format:
                (0) edge_name, (1) start_vertex, (2) end_vertex, (3) weight
        weights (dict): Dictionary of vertex (max) weights, keyed by vertex name
        parent_dict (dict): Dictionary of parent giving max path weight to each vertex.
                Keyed by vertex, value is parent.
        start_vertex (str): Name of the first vertex in the path

    Returns:
        path_score (float): The total weight of the path; also the weight
                of the final vertex
        path (list): Ordered list of vertices, from end of path to beginning
        edge_path (list): Ordered list of edges, from beginning of path to end
        start_vertex (str): Name of the first vertex in the path
        end_vertex (str): Name of the last vertex in the path


    '''

    # Get total score for path
    path_score = weights[vertex]
    # Populate path
    path = [vertex]
    check_beginning = False
    path_counter = 0
    while not check_beginning:
        parent = highestParent(path[path_counter], parent_dict)
        path.append(parent)
        path_counter += 1
        if start_vertex == None and weights[parent] == 0:
            check_beginning = True
        elif start_vertex == parent:
            check_beginning = True

    # Parse path
    first_vertex = path[-1]
    last_vertex = path[0]

    # Translate to edges
    edge_path = []
    for i in range(len(path)):
        j = i + 1
        for e in edges:
            if e[2] == path[-j] and e[1] == path[-i]:
                edge_path.append(e[0])

    return path_score, path, edge_path, first_vertex, last_vertex

def main():
    '''Main body of script.'''
    # Start timer
    start_time = time.time()

    # Parse command line arguments
    parser = parseArgs()
    print("\nExecuting {0} at {1}".format(parser.prog, time.asctime()))
    args = parser.parse_args()
    print("\nParsed the following arguments:\n\t{0}".format(
            '\n\t'.join(['{0} = {1}'.format(arg, val) for (arg, val) in
            vars(args).items()])))

    # Read in graph
    print('\nReading in graph file: {0}'.format(args.graph))
    vertices, edges, start_vertex, end_vertex = readGraph(args.graph)
    print('File succesfully read.')

    if start_vertex == None:
        print('\nNo start vertex constraint.')
    else:
        print('\nStart vertex constrained.')
    if end_vertex == None:
        print('No end vertex constraint.')
    else:
        print('End vertex constrained.')

    # Calculate weights
    # If no start_vertex constraint:
    if start_vertex == None:
        start_index = 0
        print('\nCalculating weights of all vertices.')

    # If start vertex constraint:
    else:
        start_index = vertices.index(start_vertex)
        print('\nCalculating weights, starting at designated start vertex: {0}'.format(start_vertex))

    # Print results
    weights, parent_dict = calculateWeights(vertices, edges, start_index, start_vertex)
    print('Weights calculated.')

    # Find highest weight path from end vertex
    print('\nFinding highest weight path.')
    max_vertex = max(weights.iterkeys(), key=lambda k: weights[k])

    # If no end vetex constraint:
    if end_vertex == None:
        print('Starting at highest weight point: {0}'.format(max_vertex))
        path_score, path, edge_path, first_vertex, last_vertex = trace(max_vertex, vertices, edges, weights, parent_dict, start_vertex)

    # If end vertex constraint:
    else:
        print('Starting at designated end vertex: {0}'.format(end_vertex))
        path_score, path, edge_path, first_vertex, last_vertex = trace(end_vertex, vertices, edges, weights, parent_dict, start_vertex)

    # Print results
    print('\nPath found:\n{0}'.format(path[::-1]))
    print('Edges used:\n{0}'.format(edge_path))
    print('Start vertex: {0}'.format(first_vertex))
    print('End vertex: {0}'.format(last_vertex))
    print('Path score: {0}'.format(path_score))

    # End timer
    end_time = time.time()
    runtime = end_time-start_time

    # Write output
    writeOutput(args.out, runtime, path_score, first_vertex, last_vertex, edge_path)

if __name__ == '__main__':
    main()
