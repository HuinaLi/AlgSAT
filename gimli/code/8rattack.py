import logging
from typing import Any
from sage.all import *
from sage.rings.polynomial.pbori.pbori import *
from sage.rings.polynomial.pbori import *
from sage.sat.boolean_polynomials import solve as solve_sat
from gimli import Gimli, get_logger

def attack_8round() -> Any:
    """attack intermediate 8 differential trail

    Returns:
        Any: the sat solution
    """
    ROUNDS = 8
    # get a logger
    logger = get_logger("attack")
    logger.info("start attack")
    # define ring variables
    diff_var = 'x'  # input differential
    input_var = 'x'  # input x
    auxiliary_var = 'u'
    # defined a polynomial ring, only s0,3 is active, the length of a word s0,3 = Gimli.z = 32
    # set a bit in s0,3 to 1 in order to get a non-trivial solution
    # so the length of diff_var is len(word) - 1
    R = declare_ring([Block(input_var,  (2*ROUNDS + 1) *Gimli.state + Gimli.z - 1)] + [auxiliary_var], globals())
    a_vars = [[R(x(Gimli.state * (2*r + 1) + i)) for i in range (Gimli.state)] for r in range(ROUNDS)]
    b_vars  = [[R(x(Gimli.state * (2*r + 2) + i)) for i in range (Gimli.state)] for r in range(ROUNDS )]
    # get a Gimli object on R
    gimli = Gimli(R)
    # difference values
    diff = [[0] * Gimli.state for i in range(ROUNDS + 1)]
    # initial diff[0]
    diff[0][gimli.index_start[0][3]:gimli.index_end[0][3] - 1] = [R(diff_var + "({})".format(i)) for i in range((2*ROUNDS + 1) *Gimli.state, (2*ROUNDS + 1) *Gimli.state + Gimli.z - 1)]
    # set a bit in diff[0] s0,3 to 1, we set the last bit of s0,3
    diff[0][gimli.index_end[0][3] - 1] = 1
    # ignore s0,3 s1,3 s2,3 in diff[1], diff[3]
    # ignore s0,2 s1,2 s2,2 in diff[5], diff[7]
    for i in range(Gimli.x):
        diff[1][gimli.index_start[i][3]:gimli.index_end[i][3]] \
        = diff[3][gimli.index_start[i][3]:gimli.index_end[i][3]] \
        = diff[5][gimli.index_start[i][2]:gimli.index_end[i][2]] \
        = diff[7][gimli.index_start[i][2]:gimli.index_end[i][2]] \
        = [-1] * Gimli.z
    # ignore s1,3 s2,3 in diff[2], s0,2 in diff[4]
    # ignore s1,2 s2,2 in diff[6], s0,3 in diff[8]
    diff[2][gimli.index_start[1][3]:gimli.index_end[1][3]] \
    = diff[2][gimli.index_start[2][3]:gimli.index_end[2][3]] \
    = diff[4][gimli.index_start[0][2]:gimli.index_end[0][2]] \
    = diff[6][gimli.index_start[1][2]:gimli.index_end[1][2]] \
    = diff[6][gimli.index_start[2][2]:gimli.index_end[2][2]] \
    = diff[8][gimli.index_start[0][3]:gimli.index_end[0][3]] \
    = [-1] * Gimli.z
    ##########################################################################################################################
    
    # initial input values and SAT claues
    # input X and diff
    X = [R(input_var + "({})".format(i)) for i in range(Gimli.state)]
    # set of SAT clauses
    Q = set()
    # add initial differential
    logger.info("add initial differential")
    for i in range(Gimli.state):
        X[i] += diff[0][i] * R(auxiliary_var)
    logger.info("start adding round difference to clauses")
    # start round function
    for r in range(ROUNDS):
        # set current round number
        current_round = r + 1
        # non-linear layer
        X = gimli.non_linear(X)
        # linear layer
        X = gimli.linear_mixing(X, 24 - current_round)
        # add const
        X = gimli.round_const(X, 24 - current_round)
        # variable subsitution - x = a * u + b
        for i in range(Gimli.state):
            a = X[i] / R(auxiliary_var)
            b = X[i] + a * R(auxiliary_var)
            # the r th round, the i th variable
            nva_r_i = a_vars[r][i]
            nvb_r_i = b_vars[r][i]
            Q.add(a + R(nva_r_i))
            Q.add(b + R(nvb_r_i))
            X[i] = R(nva_r_i) * R(auxiliary_var) + R(nvb_r_i)
        # add differential after a round
        logger.info("start adding round {}".format(r + 1))
        for i in range(Gimli.state):
            # if difference bit is 1, we add x / u + 1
            # else add x / u
            if diff[r + 1][i] == 1:
                d = X[i] / R(auxiliary_var)
                if d == 1:
                    # meaning it is True, no need to add to Q
                    pass
                elif d == 0:
                    logger.warning("Impossible")
                    exit(0)
                else:
                    Q.add(X[i] / R(auxiliary_var) + 1)
            elif diff[r + 1][i] == 0:
                d = X[i] / R(auxiliary_var)
                if d == 0:
                    # meaning it is True, no need to add to Q
                    pass
                elif d == 1:
                    logger.warning("Impossible")
                    exit(0)
                else:
                    Q.add(X[i] / R(auxiliary_var))
            else:
                # ignore -1
                pass
    for q in Q:
        print(q)
    
    logger.info("end adding ")
        
    

if __name__ == '__main__':
    attack_8round()