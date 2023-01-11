from __future__ import division
from sage.all import *
from copy import copy, deepcopy
from sage.rings.polynomial.pbori.pbori import *
from sage.rings.polynomial.pbori import *
from random import randint
from sage.sat.boolean_polynomials import solve as solve_sat
from sage.sat.converters.polybori import CNFEncoder
from sage.sat.solvers.dimacs import DIMACS
import sys
import logging

# create logger
logger = logging.getLogger("4rkeccak_1600")
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
state = 1600
lane_z = state // 25
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
    ROUNDS = 4
    state = 1600
    lane_z = state // 25
    R = declare_ring( [ Block( 'x', (4*ROUNDS-4)*state ),'u' ], globals() )
    a_vars  = [[R(x(state*r + i)) for i in range(state)] for r in range(2*ROUNDS -2)]
    b_vars  = [[R(x(state*(r + 2*ROUNDS -2) + i)) for i in range (state)] for r in range(2*ROUNDS -2)]      
    # numvars
    # a_vars[0]      v 1    -  1600
    # a_vars[1]      v 1601 -  3200
    # a_vars[2]      v 3201 -  4800
    # a_vars[3]      v 4801 -  6400
    # a_vars[4]      v 6401 -  8000
    # a_vars[5]      v 8001 -  9600
    # b_vars[0]      v 9601 -  11200
    # b_vars[1]      v 11201 - 12800
    # b_vars[2]      v 12801 - 14400
    # b_vars[3]      v 14401 - 16000
    # b_vars[4]      v 16001 - 17600
    # b_vars[5]      v 17601 - 19200
    diff = [[0] * state  for i in range(ROUNDS)]
    ######### diff pre #############
    
    # b_{0}_add_start(equal to diff[0])
    for i in range(5):
        for j in range(5):
            for k in range(64):
                if (i,j,k) in [(0,0,2), (0,1,2), (0,2,2), (0,3,2), (0,4,2), 
                                (2,3,60), (2,4,60), (2,4,1), 
                                (3,0,1), (3,1,1), (3,2,1), (3,3,1), 
                                (4,0,3), (4,1,3), (4,2,3), (4,4,3)]:
                    diff[0][64*(i + 5 * j) + k] = 1
                else:
                    diff[0][64*(i + 5 * j) + k] = 0
    
    # b_{1}_add_start(equal to diff[1])
    for i in range(5):
        for j in range(5):
            for k in range(64):
                if (i,j,k) in[(2,3,11), (3,2,11), (3,3,11), (4,1,57), (4,3,57)]:   
                    diff[1][64 * (i + 5 * j) + k] = 1
                else:
                    diff[1][64 * (i + 5 * j) + k] = 0
    
    # b_{2}_add_start(equal to diff[2])
    for i in range(5):
        for j in range(5):
            for k in range(64):
                if (i,j,k) in [(1,1,13), (2,2,36), (3,0,32), (3,2,1)]:    
                    diff[2][64 * (i + 5 * j) + k] = 1
                else:
                    diff[2][64 * (i + 5 * j) + k] = 0
    
    # b_{3}_add_start(equal to diff[3])
    for i in range(5):
        for j in range(5):
            for k in range(64):
                if (i,j,k) in [(0, 0, 14), (1, 0, 57), (1, 0, 17), (2, 0, 56), (2, 0, 45), (2, 0, 12), (2, 0, 15), (3, 0, 57), (4, 0, 46), (4, 0, 15), (0, 1, 60), (0, 1, 0), (1, 1, 52), (1, 1, 21), (2, 1, 17), (3, 1, 18), (4, 1, 63), (4, 1, 30), (4, 1, 10), (0, 2, 38), (1, 2, 39), (1, 2, 19), (1, 2, 8), (2, 2, 61), (2, 2, 26), (3, 2, 40), (3, 2, 9), (4, 2, 32), (0, 3, 59), (0, 3, 28), (1, 3, 50), (2, 3, 47), (3, 3, 48), (3, 3, 28), (3, 3, 17), (4, 3, 28), (0, 4, 31), (0, 4, 11), (0, 4, 0), (1, 4, 27), (2, 4, 40), (2, 4, 7), (3, 4, 55), (4, 4, 39)]:
                    diff[3][64 * (i + 5 * j) + k] = 1
                else:
                    diff[3][64 * (i + 5 * j) + k] = 0
    Q = set()
    ################################ Start Get ANFs ###################################################
    r = 0
    for i in range(state):
        Q.add(a_vars[2*r][i] + b_vars[2*r][i] + diff[2*r][i])
    for r in range(1,ROUNDS):
        # ignore non-linear operation a/b_vars[even] & a/b_vars[odd]
        X = a_vars[2*r-1].copy() 
        Y = b_vars[2*r-1].copy()
        #a/b_1 --> a/b_2
        X = addConst(X, r)
        X = theta(X)
        X = rhoPi(X)
        Y = addConst(Y, r)
        Y = theta(Y)
        Y = rhoPi(Y)
        for i in range( state ):
            Q.add(X[i] + Y[i] + diff[r][i]) 
        # a/b_2
        if r < ROUNDS -1 :
            for i in range(state):
                Q.add(X[i] + a_vars[2*r][i] ) 
                Q.add(Y[i] + b_vars[2*r][i] )      
    
    solver = DIMACS(filename = "keccak_linear")
    e = CNFEncoder(solver, R)
    e(list(Q))
    solver.write()
    
    
    