import random


def readFile(filepath):
    values = dict()
    with open(filepath, encoding='utf-8') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("["):        # find section
                idx2 = line.index("]")
                section = line[1:idx2]      # get section name
                values[section] = []        # create an empty list for key: section_name
            else:
                try:
                    values[section].append(eval(line))      # if the line can be evaluated
                except:
                    values[section].append(line.replace('\n', ''))
    global table
    table = values["ProbabilityTable"]      # fill the global table variable
    return values


# global table variable to avoid passing it for recursions
table = None


# normalizing function
def normalize(arr):
    total = sum(arr)
    normalized = [round(x / total, 3) for x in arr]
    return normalized


# check if l1 and l2 have any common variable - for enumeration
def check_intersection(l1, l2):
    s1 = set(l1)
    s2 = set(l2)
    inter = s1.intersection(s2)
    return len(inter) != 0


def get_prob(value, parents, assignments, value_assign):
    if table is None:       # error in file reading
        raise Exception("Table value is None!!")
    for entry in table:
        if value == entry[0] and parents == entry[1]:
            if not parents:     # no parents no dependencies
                prob = next(iter(entry[2]))
                if not value_assign:    # get true probability on value: False
                    prob = 1 - prob
                return prob
            else:
                # make assignment variable as key for the dict
                if len(assignments) == 1:
                    assignments = assignments[0]
                else:
                    assignments = tuple(assignments)
                prob = entry[2][assignments]
                if not value_assign:    # get true probability on value: False
                    prob = 1 - prob
                return prob


# get the parents of given node
def get_parents(value, paths):
    for path in paths:
        if path[1] == value:
            return path[0]
    return []


# get the parents of given node(for gibbs ask)
def get_children(value, paths):
    children = []
    for path in paths:
        if value in path[0]:
            children.append(path[1])
    return children


# enumeration ask algorithm as in the pdf
def enumeration_ask(data):
    ret_arr = []
    X, e = data["Query"][0][0], data["Query"][0][1]
    paths = data["Paths"]
    for val in (True, False):       # give X all possible values then return the probabilites
        send_e = e.copy()
        send_vars = data["BayesNetNodes"].copy()
        send_e[X] = val
        ret_arr.append(enumerate_all(send_vars, send_e, paths))
    normalized = normalize(ret_arr)
    return tuple(normalized)


# enumerate all algorithm as in the pdf
def enumerate_all(vars, e, paths):
    if not vars:
        return 1
    var = vars.pop(0)
    parents = get_parents(var, paths)
    # check if a child come before its parent
    if check_intersection(parents, vars):  # add to the last position and pass
        vars.append(var)
        return enumerate_all(vars, e, paths)
    assignments = [e[parent] for parent in parents]
    if var in e:        # get the probability if it is assigned in e
        send_e = e.copy()
        send_vars = vars.copy()
        return get_prob(var, parents, assignments, send_e[var]) * enumerate_all(send_vars, send_e, paths)
    else:
        sum_values = []
        for val in (True, False):       # assign var in e then get the probabilities
            send_e = e.copy()
            send_vars = vars.copy()
            send_e[var] = val
            sum_values.append(
                get_prob(var, parents, assignments, send_e[var]) * enumerate_all(send_vars, send_e, paths))
        return sum(sum_values)


# calculate true-false probability of a value in its markov blanket(parents, children, parents of childs)
def calculate_blanket(value, path, e):
    arr = []
    parents = get_parents(value, path)
    children = get_children(value, path)
    assignments = [e[parent] for parent in parents]
    for val in (True, False):
        e[value] = val
        prob = get_prob(value, parents, assignments, val)
        for child in children:
            c_parents = get_parents(child, path)
            c_assignments = [e[parent] for parent in c_parents]
            prob *= get_prob(child, c_parents, c_assignments, e[child])
        arr.append(prob)
    return tuple(arr)


# gibbs ask algorithm as in the pdf
def gibbs_ask(data, N):
    if N == 0:
        print("N can not be 0!")
        return None
    random.seed(10)
    count = {True: 0, False: 0}
    X, e = data["Query"][0][0], data["Query"][0][1]
    z = [x if x not in e else -1 for x in data["BayesNetNodes"]]            # z = all nodes - e
    z = list(filter(lambda a: a != -1, z))
    x = e.copy()
    # assign nodes randomly
    for elem in z:
        x[elem] = random.choice((True, False))
    paths = data["Paths"]

    for i in range(N):
        # iterate through nodes than count X=True and X=False then return normalized counts as probability
        for zi in z:
            choice = random.choices(population=[True, False], weights=calculate_blanket(zi, paths, x), k=1)
            x[zi] = choice[0]
        count[x[X]] += 1
    normalized = normalize(list(count.values()))
    return tuple(normalized)


def DoInference(method_name, problem_file_name, num_iteration):
    data = readFile(problem_file_name)
    if method_name == "ENUMERATION":
        result = enumeration_ask(data)
    elif method_name == "GIBBS":
        result = gibbs_ask(data, num_iteration)
    else:
        print("Wrong method name!")
        return None
    return result


if __name__ == '__main__':
    # result = enumeration_ask(data)
    true, false = DoInference("GIBBS", "query2.txt", 0)
