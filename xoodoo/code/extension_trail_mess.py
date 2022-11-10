"""
    Xoodoo: a 384-bit permutation
    lane: 32 bits
    state: 384 bits = 3 * 4 words
    0 <= x <= 3, 0 <= y <= 2, 0 <= z <= 31
    index of a bit [x,y,z]: 32 * (x + 4 * y) + z
"""
import logging
import argparse

from sage.all import *
from sage.rings.polynomial.pbori.pbori import *
from sage.rings.polynomial.pbori import *
from sage.sat.boolean_polynomials import solve as solve_sat

from read_dc import read_sol_i

# create logger
logger = logging.getLogger('xoodoo_mess')
logger.setLevel(logging.DEBUG)
# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
# create formatter
formatter = logging.Formatter('c: %(asctime)s - %(name)s - line:%(lineno)d - %(levelname)s - %(message)s')
# add formatter to ch
ch.setFormatter(formatter)
# add ch to logger
logger.addHandler(ch)

def SingleMatrix(X, r0, r1):
  Y = []
  for i in range(4):
      for j in range(32):
          Y.append(X[32 * ((i - r0 + 4) % 4) + (j - r1 + 32) % 32])
  return Y

def theta(X):
    E1 = [R(0) for i in range(128)]
    E2 = [R(0) for i in range(128)]
    P = []
    Y = []
    for i in range(128):
        P.append(X[i] + X[i + 128] + X[i + 256])
    E1[0:128] = SingleMatrix(P, 1, 5)
    E2[0:128] = SingleMatrix(P, 1, 14)

    for j in range(128):
        Y.append(X[j] + E1[j] + E2[j])
        
    for j in range(128):
        Y.append(X[j + 128] + E1[j] + E2[j])

    for j in range(128):
        Y.append(X[j + 256] + E1[j] + E2[j])
    return Y

def rhowest(X):
    X[0:128] = SingleMatrix(X[0: 128], 0, 0)
    X[128:256] = SingleMatrix(X[128:256], 1, 0)
    X[256:384] = SingleMatrix(X[256:384], 0, 11)

    return (X)

def addConst ( X, r ):
    constant = [ 0x0058, 0x0038, 0x03c0, 0x00d0, 0x0120, 0x0014, 0x0060, 0x002c, 0x0380, 0x00f0,
            0x01a0, 0x0012 ]
    for i in range(16):
        if constant[r] >> i  & 0x1:
            X[i] += 1
    return X

def SingleSbox(x0, x1, x2):
    y0 = x0 + (1 + x1) * x2
    y1 = x1 + (1 + x2) * x0
    y2 = x2 + (1 + x0) * x1
    return y0, y1, y2

def chi(A):
    B = [R(0) for i in range(384)]
    for j in range(128):
        B[0 + j], B[128 + j], B[256 + j] = SingleSbox(A[0 + j], A[128 + j], A[256 + j])
    return B

def rhoeast(X):
    X[0:128] = SingleMatrix(X[0: 128], 0, 0)
    X[128:256] = SingleMatrix(X[128:256], 0, 1)
    X[256:384] = SingleMatrix(X[256:384], 2, 8)
    return (X)

def round(X):
    for i in range(2):
       X = theta(X)
       X = rhowest(X)
       X = addConst(X, i)
       X = chi(X)
       X = rhoeast(X)
    return X

def xoodoo_mess(diff: list, ROUNDS: int) -> None:
    """verify a differential trail

    Args:
        diff (list): a differential trail
        ROUNDS (int): number of rounds
    """
    # number of vars
    # X: 384                                  v 1 - 384
    # a_vars[0]:                              v 385 - 768
    # b_vars[0]:                              v 769 - 1152
    # a_vars[1]:                              v 1153 - 1536
    # b_vars[1]:                              v 1537 - 1920
    # a_vars[2]:                              v 1921 - 2304
    # b_vars[2]:                              v 2305 - 2688
    # a_vars[3]:                              v 2689 - 3072
    # b_vars[3]:                              v 3073 - 3456
    # a_vars[4]:                              v 3457 - 3840
    # b_vars[4]:                              v 3841 - 4224   
    # verify_3r: b0-->a1-->b1-->a2-->b2-->a3 weight = w(a1 + a2 + b2)
    # verify_4r: b0-->a1-->b1-->a2-->b2-->a3-->b3-->a4    weight = w(a1 + a2 + a3 + b3)
    # define a ring
    R = declare_ring([Block('x', 2*(ROUNDS-1)*384 + 384), 'u'], globals())
    # message input x
    X = [R(x(i)) for i in range(384)]
    # auxiliary vars
    a_vars = [[R(x(384*(2*r + 1) + i)) for i in range(384)] for r in range(ROUNDS - 1)]
    b_vars = [[R(x(384*(2*r + 2) + i)) for i in range (384)] for r in range(ROUNDS - 1)]
    Q = set()
    ################################ Start Get ANFs ##########################################
    # 3r: b0-->a1-->b1-->a2-->b2-->a3              weight = w(a1 + a2 + b2)
    # 4r: b0-->a1-->b1-->a2-->b2-->a3-->b3-->a4    weight = w(a1 + a2 + a3 + b3)
    for i in range(384):
        X[i] += diff[0][i] * R(u)
    for r in range(ROUNDS):
        # b_{r}-->a_{r+1}
        X = chi(X)
        # add diff[r+1]
        for i in range(384):
            if diff[r+1][i] == 1:
                d = X[i] / R(u) 
                if d == 1:
                    pass
                elif d == 0:
                    logger.info("Impossible")
                    exit(0)
                else:
                    Q.add(X[i]/R(u) + 1) 
            else:
                d = X[i] / R(u) 
                if d == 0:
                    pass
                elif d == 1:
                    logger.info("Impossible")
                    exit(0)
                else:
                    Q.add(X[i]/R(u))
        
        # variable substitution after chi and add diff
        ## compute a_{r+1} = a_vars[r-1][i]
        if r < ROUNDS -1:
            for i in range(384):
                a = X[i] / R(u)
                b = X[i] + a * R(u)
                # the r th round, the i th variable
                # x = a * u + b
                Q.add(a + a_vars[r][i])
                Q.add(b + b_vars[r][i])
                X[i] = a_vars[r][i] * R(u) + b_vars[r][i]
        
            # a_{r}-->b_{r}
            X = rhoeast(X)
            X = theta(X)
            X = rhowest(X)
            X = addConst(X, r + 1)
    
    for q in Q:
        print(q)
    
    logger.info("end adding")

if __name__ == '__main__':
    parse = argparse.ArgumentParser(description="verify the i-th solution in solution file")
    parse.add_argument("-f", "--file", type=str, 
                    default="/home/user/anonymous/bosphorus/xoodoo/3rxoodoo/3rxoomess/dc.txt", help="solution file")
    parse.add_argument("-r", "--rounds", type=int, default=3, help="rounds of the solution")
    parse.add_argument("-i", "--ith", type=int, default=0, help="the i-th solution to verify")
    args = parse.parse_args()
    sol_file = args.file
    ROUNDS = args.rounds
    i = args.ith

    # read diff from solution file
    diffi = read_sol_i(sol_file, i, ROUNDS, ROUNDS+1)
    if ROUNDS+1 != len(diffi):
        raise ValueError("reading {}-th solution wrong".format(i))
    xoodoo_mess(diffi, ROUNDS)