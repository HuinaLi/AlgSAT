from sage.all import *
from sage.rings.polynomial.pbori.pbori import *
from sage.rings.polynomial.pbori import *
from sage.sat.boolean_polynomials import solve as solve_sat
import logging
import argparse

from read_trails import read_trails

# create logger
logger = logging.getLogger("4rkeccak_w133")
logger.setLevel(logging.DEBUG)

# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)

# create formatter
formatter = logging.Formatter("c: %(asctime)s - %(name)s - %(levelname)s - %(message)s")

# add formatter to ch
ch.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)
 

"""
Keccak-f[1600]: a state-bit permutation
little endian: for example - 0x00000004 = [0, 0, 0, 0, 0, 1, 0, 0, ..., 0] i.e. 04000000
                                                        ^
                                                        |
                                        lowest byte at lowest significant
But hex2vector() function does not convert to little endian, and function's output is [0, 0, 1, 0, 0, 0, 0, 0, ..., 0],

state: 1600 bits = 5 * 5 * 64 
We use a list to denote a state/word and the element of the list is a bit
Keccak uses x,y,z to denote each bit, 0 <= x <= 4, 0 <= y <= 4, 0 <= z <= 63
index of each bit: 64 * (5 * y + x) + z
p0: 0-320
p1: 320-640
...
p4: 1280-1600
<<<: cycle left shift
<<: non-cycle left shift
""" 

def SinglePlane(X, dx, dz):
    P = []
    for i in range(5):
        for j in range(lane_z):
            P.append(X[lane_z * ((i - dx) % 5) + (j - dz) % lane_z])
    return P

def theta(X):
    E1 = [R(0) for i in range(5 * lane_z)]
    E2 = [R(0) for i in range(5 * lane_z)]
    P = []
    Y = []
    for i in range(5 * lane_z):
        P.append(X[i] + X[i + 5 * lane_z] + X[i + 10*lane_z] + X[i + 15 * lane_z] + X[i + 20 * lane_z])
    E1[0:320] = SinglePlane(P, 1, 0)
    E2[0:320] = SinglePlane(P, -1, 1)

    for j in range(5 * lane_z):
        Y.append(X[j] + E1[j] + E2[j])
        
    for j in range(5 * lane_z):
        Y.append(X[j + 5 * lane_z] + E1[j] + E2[j])

    for j in range(5 * lane_z):
        Y.append(X[j + 10*lane_z] + E1[j] + E2[j])

    for j in range(5 * lane_z):
        Y.append(X[j + 15 * lane_z] + E1[j] + E2[j])

    for j in range(5 * lane_z):
        Y.append(X[j + 20*lane_z] + E1[j] + E2[j])
    return Y

def SingleLane(X, dz):
    Y = []
    for i in range(lane_z):
        Y.append(X[(i - dz + lane_z) % lane_z])
    return Y
                                                                                    
def rhoPi(X):
    ## Rotation offsets
    r=[[0,    36,     3,    41,    18]    ,
        [1,    44,    10,    45,     2]    ,
        [62,    6,    43,    15,    61]    ,
        [28,   55,    25,    21,    56]    ,
        [27,   20,    39,     8,    14]    ]
    for y in range(5):
        for x in range(5):
            X[(5 * lane_z * y + lane_z * x):(5 * lane_z * y + lane_z * x + lane_z)] = SingleLane(X[(5 * lane_z * y + lane_z * x):(5 * lane_z * y + lane_z * x + lane_z)], r[x][y])

    Y = [R(0) for i in range(state)]
    for y in range(5):
        for x in range(5):
            for z in range(lane_z):
                Y[(5 * lane_z * ((2*x+3*y)%5) + lane_z * y) +  z] = X[lane_z * (5 * y + x) + z]
    return(Y)

def SingleSbox(x0, x1, x2, x3, x4):
    y0 = x0 + (1 + x1) * x2
    y1 = x1 + (1 + x2) * x3
    y2 = x2 + (1 + x3) * x4
    y3 = x3 + (1 + x4) * x0
    y4 = x4 + (1 + x0) * x1
    return y0, y1, y2, y3, y4


def sbox(A):
    """state bits sbox

    Args:
        Y (list): state bits input

    Returns:
        list: state bits output
    """
    B = [R(0) for i in range(state)]
    # 5 bits as a block, each block uses a 5-bits sbox
    for z in range(lane_z):
        for y in range(5):
            B[0 + 5 * lane_z * y + z], B[lane_z + 5 * lane_z * y + z], B[2*lane_z + 5 * lane_z * y + z], B[3*lane_z +5 * lane_z * y + z] , B[4*lane_z + 5 * lane_z * y + z] = SingleSbox(A[0 + 5 * lane_z * y + z], A[lane_z + 5 * lane_z * y + z], A[2 * lane_z + 5 * lane_z * y + z], A[3*lane_z + 5 * lane_z * y + z], A[4*lane_z+ 5 * lane_z* y + z])
    return B

def addConst ( X, r ):
    constant = [0x0000000000000001,
        0x0000000000008082,
        0x800000000000808A,
        0x8000000080008000,
        0x000000000000808B,
        0x0000000080000001,
        0x8000000080008081,
        0x8000000000008009,
        0x000000000000008A,
        0x0000000000000088,
        0x0000000080008009,
        0x000000008000000A,
        0x000000008000808B,
        0x800000000000008B,
        0x8000000000008089,
        0x8000000000008003,
        0x8000000000008002,
        0x8000000000000080,
        0x000000000000800A,
        0x800000008000000A,
        0x8000000080008081,
        0x8000000000008080,
        0x0000000080000001,
        0x8000000080008008]
    for i in range(lane_z):
        if constant[r] >> i  & 0x1:
            X[i] += 1
    return X

def round(X,r):
    for i in range(r):
        X = theta(X)
        X = rhoPi(X)
        X = sbox(X)
        X = addConst(X, i)
    return X

def checkround(X):
    # for keccak-f[1600]
    X = round(X, 24)
    for i in range(25):
        print('\n')
        for j in range(lane_z):
            print(X[lane_z*i + j], end=' ')
    return X


if __name__ == '__main__':
    parse = argparse.ArgumentParser(description="check trails")
    parse.add_argument("-f", "--file", type=str, help="trail file path")
    parse.add_argument("-r", "--rounds", type=int, default=4, help="file path")
    args = parse.parse_args()
    filepath = args.file
    ROUNDS = args.rounds
    state = 1600
    lane_z = state // 25
    trails = read_trails(filepath, ROUNDS)
    trail = trails[0]

    # number of vars
    # X: state                                 v 1 - 1600
    # a_vars[0]:                              v 1601 - 3200
    # b_vars[0]:                              v 3201 - 4800
    # a_vars[1]:                              v 4801 - 6400
    # b_vars[1]:                              v 5401 - 8000
    # a_vars[2]:                              v 8001 - 9600
    # b_vars[2]:                              v 9601 - 11200
    # 4r: b0-->b1-->b2-->b3
    R = declare_ring([Block('x', (2*ROUNDS - 1)*state ), 'u'], globals())
    X = [R(x(i)) for i in range(state)]
    a_vars = [[R(x(state*(2*r + 1) + i)) for i in range(state)] for r in range(ROUNDS-1)]
    b_vars = [[R(x(state*(2*r + 2) + i)) for i in range (state)] for r in range(ROUNDS-1)]       
    diff = [[0] * state  for i in range(ROUNDS+1)]
    ######### diff pre #############
    for r in range(ROUNDS+1):
        # b_{r}_add_start(equal to diff[r])
        for i in range(5):
            for j in range(5):
                for k in range(lane_z):
                    if (i,j,k) in trail[r]:
                        diff[r][lane_z*(i + 5 * j) + k] = 1
                    else:
                        diff[r][lane_z*(i + 5 * j) + k] = 0
    ######## Initialization ########
    for i in range(state):
        X[i] += diff[0][i] * R(u)
    Q = set()
    ######## Start Add #############
    for r in range(1, ROUNDS):
        X = sbox(X)
        for i in range( state ):
            a = X[i] / R(u)
            b = X[i] + a * R(u)
            # the r th round, the i th variable
            # x = a * u + b
            Q.add(a + a_vars[r-1][i])
            Q.add(b + b_vars[r-1][i])
            X[i] = a_vars[r-1][i] * R(u) + b_vars[r-1][i]
        X = addConst(X, r)
        X = theta(X)
        X = rhoPi(X)
        for i in range(state):
            if diff[r][i] == 1:
                d = X[i] / R(u) 
                if d == 1:
                    pass
                elif d == 0:
                    print ( diff[r][i], d )
                    print( "Impossible" )
                    exit(0)
                else:
                    Q.add(X[i]/R(u) + 1) 
            else:
                d = X[i] / R(u) 
                if d == 0:
                    pass
                elif d == 1:
                    print ( diff[r][i], d )
                    print( "Impossible" )
                    exit(0)
                else:
                    Q.add(X[i]/R(u) ) 
   
    """
    for q in Q:
        print(q)
    """
    logger.info( " start solve " )
    s = solve_sat ( list ( Q ))
    logger.info( s )
    logger.info("finished")