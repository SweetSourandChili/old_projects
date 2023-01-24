# function calculates manhattan distance between 2 2D points
def manhattan_distance(point1, point2):
    return abs(point1[0] - point2[0]) + abs(point1[1] - point2[1])


# function is used to understand if a customer found
def letter_check(nodes, letter, index):
    return nodes[index]['name'].startswith(letter)


# returns a list with a dictionary where each element is a node in the graph
def construct_nodes(environment):
    length = len(environment)
    nodes = []
    line_num = 0
    customer_count = 0
    final = None
    for line in environment:
        for i in range(length):
            if line[i] == 'C':
                customer_name = 'C' + str(customer_count)
                customer_count += 1
                nodes.append({'name': customer_name, 'index': [line_num, i]})
            elif line[i] == 'S':
                nodes.insert(0, {'name': 'S', 'index': [line_num, i]})
            elif line[i] == 'F':
                final = {'name': 'F', 'index': [line_num, i]}
        line_num += 1
    nodes.append(final)
    return nodes


# function constructs adjacency matrix
def construct_matrix(nodes):
    length = len(nodes)
    matrix = [[0 for col in range(length)] for row in range(length)]
    for i in range(length):
        for j in range(length):
            if i == j:
                continue
            index_i = nodes[i]['index']
            index_j = nodes[j]['index']
            matrix[i][j] = manhattan_distance(index_i, index_j)
    return matrix


def dfs(nodes, matrix, start_node_index, min_req):
    done = False
    length = len(nodes)
    stack = [start_node_index]                      # stack data structure for dfs - append start node
    visited = [False] * length                      # start all nodes as not visited
    result = []
    while stack:
        elem = stack.pop()
        if not visited[elem]:
            result.append(nodes[elem]['index'])     # add current node to result
            visited[elem] = True                    # make the node visited
            if letter_check(nodes, 'C', elem):      # check if the node is a customer
                min_req -= 1
                if min_req <= 0:
                    done = True
                    break
        for j in range(length - 1):                 # add neighbors to da stack
            if matrix[elem][j] != 0 and not visited[j]:
                stack.append(j)

    result.append(nodes[length - 1]['index'])       # add final to the result
    if done:
        return result
    else:
        return None


def bfs(nodes, matrix, start_node_index, min_req):
    done = False
    length = len(nodes)
    queue = [start_node_index]                      # queue data structure for bfs - append start node
    visited = [False] * length                      # start all nodes as not visited
    result = []
    while queue:
        elem = queue.pop(0)                         # pop 0th index - FIFO
        result.append(nodes[elem]['index'])         # add current node to result
        visited[elem] = True                        # make the node visited
        if letter_check(nodes, 'C', elem):          # check if the node is a customer
            min_req -= 1
            if min_req <= 0:
                done = True
                break
        for j in range(length - 1):                 # add neighbors to da queue
            value = matrix[elem][j]
            if value != 0 and not visited[j]:       # if the node is not itself and not visited yet
                queue.append(j)
                visited[j] = True
    result.append(nodes[length - 1]['index'])
    if done:
        return result
    else:
        return None


def ucs(nodes, matrix, start_node_index, min_req):
    length = len(nodes)
    if length-2 < min_req:                          # there is no enough customers for required packages
        return None
    explored = []
    queue = [{'index': start_node_index, 'cost': 0, 'path': [start_node_index]}]    # nodes are kept with the information - currentPath, currentCost, nodeIndex
    while True:
        if not queue:
            break

        queue = sorted(queue, key=lambda dict: dict['cost'])        # sort the queue according to cost since we need priority queue
        curr = queue.pop(0)

        if len(curr['path']) >= min_req + 1:                        # check if the min_requirenment is met and add to list
            if not any(elem['index'] == i for elem in explored):
                explored.append(curr)
            continue

        for i in range(length-1):
            index = curr['index']
            if matrix[index][i] != 0:                               # if it is not the node itself
                cost = matrix[index][i] + curr['cost']
                path = curr['path']
                # index is not in the queue or explored
                if not any(elem['index'] == i for elem in queue) and not any(elem['index'] == i for elem in explored):
                    if i not in path:
                        queue.append({'index': i, 'cost': cost, 'path': path + [i]})        # add connected nodes with their costs
                elif any(elem['index'] == i for elem in queue):
                    ind = next((x for (x, d) in enumerate(queue) if d["index"] == i), None)
                    if (cost < queue[ind]['cost'] or len(queue[ind]['path']) < min_req + 1) and not any(elem['index'] == i for elem in explored):   # append if cost is smaller
                        if i not in path:
                            queue.append({'index': i, 'cost': cost, 'path': path + [i]})

    min_dist_final = float('inf')       # add final distances and select the path with the minimum cost
    min_dist_index = 0
    for i in range(len(explored)):
        distance = explored[i]['cost'] + matrix[explored[i]['index']][length-1]
        if distance < min_dist_final:
            min_dist_final = distance
            min_dist_index = i

    path = explored[min_dist_index]['path'] + [length-1]
    path_list = []
    for ind in path:
        path_list.append(nodes[ind]['index'])
    return path_list


# the functions can also work if the start node is not specified with S
# with minor changes the DFS, BFS, UCS in this file then it also work in a case that if all the nodes were not connected
def UnInformedSearch(method_name, problem_file_name):
    file = open(problem_file_name, "r")
    contents = file.read()
    dict = eval(contents)
    file.close()
    min_requirenment = dict['min']
    environment = dict['env']
    nodes = construct_nodes(environment)        # construct nodes and matrix then pass them to required algorithm
    matrix = construct_matrix(nodes)
    result = None
    if method_name == "DFS":
        result = dfs(nodes, matrix, 0, min_requirenment)
    elif method_name == "BFS":
        result = bfs(nodes, matrix, 0, min_requirenment)
    elif method_name == "UCS":
        result = ucs(nodes, matrix, 0, min_requirenment)
    return result
