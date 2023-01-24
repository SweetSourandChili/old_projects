# check if the game is finished - break recursion
def is_finished(nim):
    return all(elem == 0 for elem in nim)


# minimax recursive function
def minimax(nim, is_max, depth, prune=False, alpha=None, beta=None):
    if is_finished(nim):        # if finished return the score according to the type of the player
        if is_max:
            return 1, depth
        else:
            return -1, depth
    # I started with -inf and inf however when I give -1,1 I can reach the numbers in the pdf
    max_score = -1 if is_max else 1
    # iterate through each chunk and each possibility then call minimax
    for idx, chunk in enumerate(nim):
        for i in range(chunk, 0, -1):
            pass_nim = nim.copy()
            pass_nim[idx] -= i
            # call minimax with the current nim board, opposite of the last player, depth and prune variables
            score, depth = minimax(pass_nim, not is_max, depth + 1, prune, alpha, beta)
            # get the score according to the player_type
            max_score = max(max_score, score) if is_max else min(max_score, score)
            # prune
            if prune:
                if is_max:
                    alpha = max(max_score, alpha)
                else:
                    beta = min(beta, max_score)
                # if alpha is greater or equal than beta, stop the recursion and return
                if alpha >= beta:
                    return max_score, depth
    # if there is no pruning return the values
    return max_score, depth


players = {'MAX': False, 'MIN': True}


def make_move(nim, player_type, prune=False, alpha=None, beta=None):
    player = players.setdefault(player_type)
    if player is None:
        print(f'There exist no player as {player_type}')
        return None

    select_chunk = 0
    number_piles = 0
    max_score = -1
    actual_depth = 0
    # for each possible move calculate the depth and return the best move with its depth
    for idx, chunk in enumerate(nim):
        for i in range(chunk, 0, -1):
            pass_nim = nim.copy()
            pass_nim[idx] -= i
            if prune:
                score, depth = minimax(pass_nim, players[player_type], 0, True, alpha, beta)
            else:
                score, depth = minimax(pass_nim, players[player_type], 0)
            if score > max_score:
                max_score = score
                select_chunk = idx
                number_piles = i
                actual_depth = depth

    # update nim board and return
    nim[select_chunk] -= number_piles
    return [tuple(nim), actual_depth]


def SolveGame(method_name, problem_file_name, player_type):
    file = open(problem_file_name, "r")
    contents = file.read()
    nim = eval(contents)
    file.close()
    retval = []
    if method_name == 'Minimax':
        retval = make_move(nim, player_type)
    elif method_name == 'AlphaBeta':
        retval = make_move(nim, player_type, True, -1, 1)
    else:
        print('Wrong method name')
    return retval


if __name__ == '__main__':
    print(SolveGame('Minimax', "nim2.txt", "MAX"))
    print(SolveGame('AlphaBeta', "nim2.txt", "MAX"))
