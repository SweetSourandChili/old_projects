def read_grid(filename):
    """
    Reads the txt file and return it as 2D array
    :param filename: Txt file that has grid values
    :return: 2D array that has grid values
    """
    f = open(filename, 'r')
    lines = [line.split() for line in f]
    return lines


def construct_nodes(lines):
    """
    Constructs every grid as nodes like a graph
    Node has properties:
        value: S, E, or . for UCS - S, E, or heuristic values for A*
    :param lines: 2D array that has grid values
    :return: created nodes 2D array, index of the start array
    """
    # node data structure = {'value': '# or . etc.','index': index, 'cost': cost, 'path': path}
    size = len(lines)
    nodes = []
    start = -1
    for i in range(size):
        innodes = []
        for j in range(size):
            value = lines[i][j]
            if value == 'S':
                start = [i, j]
            elem = {'value': value, 'index': [i, j], 'cost': 1, 'path': [[i, j]], 'total': None}
            innodes.append(elem)
        nodes.append(innodes)
    return nodes, start


def correct_path(path):
    new_path = [(p[1], p[0]) for p in path[::-1]]
    return new_path


def ucs(nodes, start):
    """
    Calculates the shortest path using uniformed cost search algorithm
    :param nodes: 2D array of nodes
    :param start: index of the start node
    :return: optimum path
    """
    length = len(nodes)
    explored = []                                               # array to keep visited nodes
    start_node = nodes[start[0]][start[1]]
    start_node['cost'] = 0
    start_node['total'] = 0
    queue = [nodes[start[0]][start[1]]]                         # add start node to the queue before entering the loop
    while True:
        if not queue:
            return None

        queue = sorted(queue, key=lambda dict: dict['cost'])    # make it a priority queue
        curr = queue.pop(0)                                     # get min_cost element
        index = curr['index']
        if curr['value'] == 'E':                                # final node is reached return the path
            return correct_path(curr['path'])

        explored.append(curr)
        # direction vectors
        # there is a connection in every node left-up-right-down direction
        # add direction index values to current index to get neighbor index
        # if the new indexes is not in the grid border then continue
        dr = (0, -1, 0, +1)
        dc = (-1, 0, +1, 0)
        for k in range(4):
            new_r = index[0] + dr[k]
            new_c = index[1] + dc[k]
            if (new_r < 0 or new_c < 0) or (new_r >= length or new_c >= length):
                continue
            child = nodes[new_r][new_c]
            if child['value'] == '#':   # hit the border
                continue
            cost = curr['cost'] + child['cost']     # update cost and path values
            path = curr['path'].copy()
            path.append([new_r, new_c])
            # the child not in queue or visited
            if not any(elem['index'] == [new_r, new_c] for elem in queue) and not any(elem['index'] == [new_r, new_c] for elem in explored):
                child['cost'] = cost
                child['path'] = path
                queue.append(child)
            # child in the queue but not visited
            elif any(elem['index'] == [new_r, new_c] for elem in queue):
                ind = next((x for (x, d) in enumerate(queue) if d["index"] == [new_r, new_c]), None)
                if cost < queue[ind].get('cost'):
                    queue[ind]['cost'] = cost
                    queue[ind]['path'] = path


def astar(nodes, start):
    """
    Logic is the same with ucs
    The main difference is to sort the priority queue according to the heuristic values(cost + given value in the grid)
    :param nodes: 2D array of nodes
    :param start: index of the start node
    :return: optimum path
    """

    length = len(nodes)
    explored = []
    start_node = nodes[start[0]][start[1]]
    start_node['cost'] = 0
    queue = [nodes[start[0]][start[1]]]  # nodes are kept with the information - currentPath, currentCost, nodeIndex

    while True:
        if not queue:
            return None

        queue = sorted(queue, key=lambda dict: dict['total'])       # priority queue
        curr = queue.pop(0)
        index = curr['index']
        if curr['value'] == 'E':                                    # final node reached - return
            return correct_path(curr['path'])
        explored.append(curr)
        # direction vectors
        # there is a connection in every node left-up-right-down direction
        # add direction index values to current index to get neighbor index
        # if the new indexes is not in the grid border then continue
        dr = (0, -1, 0, +1)
        dc = (-1, 0, +1, 0)
        for k in range(4):
            new_r = index[0] + dr[k]
            new_c = index[1] + dc[k]
            if (new_r < 0 or new_c < 0) or (new_r >= length or new_c >= length):
                continue
            child = nodes[new_r][new_c]
            if child['value'] == '#':       # hit the border
                continue
            # calculate heuristic + cost and add to total
            # update cost,total and path values
            cost = curr['cost'] + child['cost']
            total = None
            if curr['value'].isnumeric():
                total = cost + int(curr['value'])
            else:
                total = cost
            path = curr['path'].copy()
            path.append([new_r, new_c])

            # the child not in queue or visited
            if not any(elem['index'] == [new_r, new_c] for elem in queue) and not any(elem['index'] == [new_r, new_c] for elem in explored):
                child['cost'] = cost
                child['path'] = path
                child['total'] = total
                queue.append(child)
                # child in the queue but not visited
            elif any(elem['index'] == [new_r, new_c] for elem in queue):
                ind = next((x for (x, d) in enumerate(queue) if d["index"] == [new_r, new_c]), None)
                if total < queue[ind].get('total'):
                    queue[ind]['cost'] = cost
                    queue[ind]['path'] = path
                    queue[ind]['total'] = total


def InformedSearch(method_name, problem_file_name):
    solution = None
    lines = read_grid(problem_file_name)
    nodes, start = construct_nodes(lines)
    if method_name == 'UCS':
        solution = ucs(nodes, start)
    elif method_name == 'AStar':
        solution = astar(nodes, start)
    else:
        raise Exception('There is no method named - {}'.format(method_name))
    return solution


if __name__ == '__main__':
    """time_calc(InformedSearch, ('UCS', 'sampleUCS.txt'))
    time_calc(InformedSearch, ('AStar', 'sampleAstar.txt'))"""

    sol = InformedSearch('UCS', 'sampleUCS.txt')
    sol2 = InformedSearch('AStar', 'sampleAstar.txt')
    print(sol)
    print(sol2)

