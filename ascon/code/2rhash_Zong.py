from __future__ import division
from sage.all import *
from copy import copy, deepcopy
from sage.rings.polynomial.pbori.pbori import *
from sage.rings.polynomial.pbori import *
from random import randint
import sys
from sage.sat.boolean_polynomials import solve as solve_sat
import logging

# create logger
logger = logging.getLogger('2rhash_Zong')
logger.setLevel(logging.DEBUG)

# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)

# create formatter
formatter = logging.Formatter('c:%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# add formatter to ch
ch.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)

def SingleMatrix( X, r0, r1 ):
    Y = []
    for i in range(64):
        Y.append( X[i] + X[( i + (64 - r0) ) % 64] + X[ ( i + (64 - r1) ) % 64 ] )
    return Y
        
def Matrix( X ):
    X[0  :64]  = SingleMatrix(X[0  : 64], 19, 28 )
    X[64 :128] = SingleMatrix(X[64 :128], 61, 39 )
    X[128:192] = SingleMatrix(X[128:192], 1,  6 )
    X[192:256] = SingleMatrix(X[192:256], 10, 17 )
    X[256:320] = SingleMatrix(X[256:320], 7, 41 )
    return X

def SingleSbox( y0, y1, y2, y3, y4 ):
    x0 = y4*y1 + y3 + y2*y1 + y2 + y1*y0 + y1 + y0
    x1 = y4 + y3*y2 + y3*y1 + y3 + y2*y1 + y2 + y1 + y0
    x2 = y4*y3 + y4 + y2 + y1 + 1
    x3 = y4*y0 + y4 + y3*y0 + y3 + y2 + y1 + y0
    x4 = y4*y1 + y4 + y3 + y1*y0 + y1
    return x0, x1, x2, x3, x4

def Sbox( Y ):
    Z = [ R(0) for i in range(320)]
    for j in range(64):
        Z[0 + j], Z[64 + j], Z[128 + j], Z[192 + j] , Z[256 + j] = SingleSbox( Y[0 + j], Y[64 + j], Y[128 + j], Y[192 + j], Y[256+j] )
    return Z

def addConst ( X, r ):
    constant = [ 0xf0, 0xe1, 0xd2, 0xc3, 0xb4, 0xa5, 0x96, 0x87, 0x78, 0x69,
            0x5a, 0x4b ]
    base = 184
    for i in range(8):
        if constant[r] >> ( 7 - i ) & 0x1:
            X[ base + i] += 1
    return X

if __name__ == '__main__':
    ROUNDS = 2
    R = declare_ring( [ Block( 'x', (2*ROUNDS + 1)*320 ),'u' ], globals() )
    X = [R(x(i)) for i in range(320)]
    a_vars  = [[R(x(320*(2*r + 1) + i)) for i in range(320)] for r in range(ROUNDS)]
    b_vars  = [[R(x(320*(2*r + 2) + i)) for i in range (320)] for r in range(ROUNDS)]       
    diff = [[0] * 320 for i in range(2*ROUNDS)]
    ####### diff_pre ############
    d0 = 0xe6765f2bfb737f78
    d3 = 0xd255739452530b86
    for i in range(64):
        diff[0][i] = d0 >> ( 63 - i ) & 0x1
        diff[3][i] = d3 >> ( 63 - i ) & 0x1
    w = [ 0 for i in range(5) ]
    w[0] = 0x00144000c0404000
    w[1] = 0xe6765f2bfb737f78
    w[2] = 0x0000000000000000
    w[3] = 0x0400000008101000
    w[4] = 0xe6621f2b3b333f78
    for i in range(5):
        for j in range(64):
            diff[1][64*i+j] = w[i] >> ( 63 - j ) & 0x1
    a = [ 0 for i in range(5) ]
    a[0] = 0x0c10400249045804
    a[1] = 0x8232408ad1246801
    a[2] = 0x0000000000000000
    a[3] = 0x0c0102000812100c
    a[4] = 0x8233428ad1366809
    for i in range(5):
        for j in range(64):
            diff[2][64*i+j] = a[i] >> ( 63 - j ) & 0x1  
    ### Initialization ######
    Q = set()
    for i in range(320):
        X[i] += diff[0][i] * R(u)
    ###########Start Add ##################
    for r in range(ROUNDS): 
        X = addConst(X,r)
        X = Sbox(X)
        for i in range( 320 ):
            a = X[i] / R(u)
            b = X[i] + a * R(u)
            # the r th round, the i th variable
            # x = a * u + b
            Q.add(a + a_vars[r][i])
            Q.add(b + b_vars[r][i])
            X[i] = a_vars[r][i] * R(u) + b_vars[r][i]
        X = Matrix( X )
        for i in range(320):
            if diff[r+2][i] == 1:
                d = X[i] / R(u) 
                if d == 1:
                    pass
                elif d == 0:
                    print ( diff[r+2][i], d )
                    print( "Impossible" )
                    exit(0)
                else:
                    Q.add(X[i]/R(u) + 1) 
            else:
                d = X[i] / R(u) 
                if d == 0:
                    pass
                elif d == 1:
                    print ( diff[r+2][i], d )
                    print( "Impossible" )
                    exit(0)
                else:
                    Q.add(X[i]/R(u) ) 
    
    for q in Q:
        print (q) 
    """
    logger.info( " start solve " )
    s = solve_sat ( list ( Q ))
    logger.info( s ) 
    logger.info("finished")
    """
