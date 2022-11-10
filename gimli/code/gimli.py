import logging
import argparse
from typing import Any

from sage.all import *
from sage.rings.polynomial.pbori.pbori import *
from sage.rings.polynomial.pbori import *
from sage.sat.boolean_polynomials import solve as solve_sat

def get_logger(msg: str ="example") -> logging.Logger:
    """get a format logger

    Args:
        msg (str, optional): description of the logger. Defaults to "example".

    Returns:
        logging.Logger: the format logger
    """
    # create logger, set message
    logger = logging.getLogger(msg)
    # clear existed handler
    logger.handlers.clear()
    # set level
    logger.setLevel(logging.DEBUG)
    # create console handler and set level to debug
    ch = logging.StreamHandler()
    # set level
    ch.setLevel(logging.DEBUG)
    # create formatter 
    formatter = logging.Formatter("c: %(asctime)s - %(name)s - %(levelname)s - %(filename)s[line:%(lineno)d]: %(message)s")
    # add formatter to ch
    ch.setFormatter(formatter)
    # add ch to logger
    logger.addHandler(ch)
    return logger

def hex2vector(n: int, length: int = 32) -> vector:
    """convert a hex number to a fixed length sage vector

    Args:
        n (int): input hex number
        length(int, optional): the length of the output vector, default is 32

    Returns:
        vector: output sage vector
    """
    # length should be the multiple of 8
    if length % 8 != 0:
        logger = get_logger()
        logger.error("length should be the multiple of 8")
        exit(1)
    # convert hex to binary list, now is big endian
    res = [int(i) for i in bin(n)[2:]]
    if len(res) > length:
        logger = get_logger()
        logger.error("number length overflow!!!")
        exit(1)
    # reverse
    res.reverse()
    # padding
    res += [0] * (length - len(res))
    # convert binary list to binary vector
    res = vector(GF(2), res)
    return res

class Gimli:
    """Gimli: a 384-bit permutation
    word: 32 bits
    little endian: for example - 0x00000004 = [0, 0, 0, 0, 0, 1, 0, 0, ..., 0] i.e. 04000000
                                                          ^
                                                          |
                                            lowest byte at lowest significant
    But hex2vector() function does not convert to little endian, and function's output is [0, 0, 1, 0, 0, 0, 0, 0, ..., 0],
    so that cycle_lshift() and noncycle_lshift() function is simpler
    state: 384 bits = 3 * 4 words
    We use a list to denote a state/word and the element of the list is a bit
    Gimli uses si,j to denote each word, 0 <= i <= 2, 0 <= j <= 3
    index of si,j: i * 128 + j * 32
    s0,0: 0-31
    s0,1: 32-63
    s0,2: 64-95
    s0,3: 96-127
    s1,0: 128-159
    s1,1: 160-191
    s1,2: 192-223
    s1,3: 224-255
    s2,0: 256-287
    s2,1: 288-319
    s2,2: 320-351
    s2,3: 352-383
    <<<: cycle left shift
    <<: non-cycle left shift
    """
    # define gimli permutation size
    x = 3
    y = 4
    z = 32  # a word length
    state = x * y * z
    def __init__(self, R: Any) -> None:
        """initial a gimli permutation object

        Args:
            R (Any): the ring where the permutation is
        """
        # define gimli parameters
        # convert 0x9e377900 to binary vector
        self.const = hex2vector(0x9e377900, self.z)
        # indexes of each word
        self.index_start = [[i * self.z * self.y + j * self.z for j in range(self.y)] for i in range(self.x)]
        self.index_end = [[self.index_start[i][j] + self.z for j in range(self.y)] for i in range(self.x)]
        self.ring = R
    
    def cycle_lshift(self, word: list, offset: int) -> list:
        """cycle left shift offset bits

        Args:
            word (list): the word(32 bits) to shift
            offset (int): shift bits

        Returns:
            list: output word
        """
        # check word's length
        if len(word) != self.z:
            logger = get_logger()
            logger.error("word length not equal to {}".format(self.z))
            exit(1)
        # cycle left shift
        return word[self.z - offset:self.z] + word[0:self.z - offset]

    def noncycle_lshift(self, word: list, offset: int) -> list:
        """noncycle left shift offset bits

        Args:
            word (list): the word(32 bits) to shift
            offset (int): shift bits

        Returns:
            list: output word
        """
        # check word's length
        if len(word) != self.z:
            logger = get_logger()
            logger.error("word length not equal to {}".format(self.z))
            exit(1)
        # non-cycle left shift
        word = [0] * offset + word[0:self.z - offset]
        return word

    def vector_and(self, word1: vector, word2: vector) -> list:
        """compute word1 & word2(bitwise and)

        Args:
            word1 (vector): input word(32 bits)
            word2 (vector): input word(32 bits)

        Returns:
            list: output word
        """
        # check word's length
        if len(word1) != self.z or len(word2) != self.z:
            logger = get_logger()
            logger.error("word length not equal to {}".format(self.z))
            exit(1)
        # and computation: a & b
        res = [word1[i] * word2[i] for i in range(self.z)]
        return res

    def vector_or(self, word1: vector, word2: vector) -> list:
        """compute word1 | word2(bitwise or)

        Args:
            word1 (vector): input word(32 bits)
            word2 (vector): input word(32 bits)

        Returns:
            list: output word
        """
        # check word's length
        if len(word1) != self.z or len(word2) != self.z:
            logger = get_logger()
            logger.error("word length not equal to {}".format(self.z))
            exit(1)
        # or computation: (a & b) ^ (a ^ b)
        res = [word1[i] * word2[i] + word1[i] + word2[i] for i in range(self.z)]
        return res

    def sp_box(self, word1: list, word2: list, word3: list) -> tuple:
        """96-bit SP-box in non-linear layer, applied to each column

        Args:
            word1 (list): first word of input column
            word2 (list): second word of input column
            word3 (list): third word of input column

        Returns:
            tuple: 3 words of output column
        """
        # check word's length
        if len(word1) != self.z or len(word2) != self.z or len(word3) != self.z:
            logger = get_logger()
            logger.error("word length not equal to {}".format(self.z))
            exit(1)
        # 3 steps
        # first x <<< 24, y <<< 9
        tmpx = vector(self.ring, self.cycle_lshift(word1, 24))
        tmpy = vector(self.ring, self.cycle_lshift(word2, 9))
        tmpz = vector(self.ring, word3)
        # second and third step
        y_and_z = self.vector_and(tmpy, tmpz)
        x_or_z = self.vector_or(tmpx, tmpz)
        x_and_y = self.vector_and(tmpx, tmpy)
        resz = tmpx + vector(self.ring, self.noncycle_lshift(word3, 1)) + vector(self.ring, self.noncycle_lshift(y_and_z, 2))
        resy = tmpy + tmpx + vector(self.ring, self.noncycle_lshift(x_or_z, 1))
        resx = tmpz + tmpy + vector(self.ring, self.noncycle_lshift(x_and_y, 3))
        return resx.list(), resy.list(), resz.list()

    def non_linear(self, X: list) -> list:
        """The non-linear layer of gimli, 
        containing 3 96-bit SP-box applied to each column, 
        and a column is 96 bits(3 words)

        Args:
            X (list): 384 bits input state

        Returns:
            list: 384 bits output state
        """
        for j in range(self.y):
            # spbox a column(3 words)
            X[self.index_start[0][j]:self.index_end[0][j]], X[self.index_start[1][j]:self.index_end[1][j]], X[self.index_start[2][j]:self.index_end[2][j]] = \
                self.sp_box(X[self.index_start[0][j]:self.index_end[0][j]], X[self.index_start[1][j]:self.index_end[1][j]], X[self.index_start[2][j]:self.index_end[2][j]])
        return X

    def linear_mixing(self, X: list, r: int) -> list:
        """The linear mixing layer of gimli containing two swap operations. 
        The linear mixing layer is used in every second round.

        Args:
            X (list): 384 bits input state
            r (int): the round number indicating which round we are running

        Returns:
            list: 384 bits output state
        """
        # small swap -- r mod 4 == 0
        if r % 4 == 0:
            X[self.index_start[0][0]:self.index_end[0][0]], X[self.index_start[0][1]:self.index_end[0][1]] = \
                X[self.index_start[0][1]:self.index_end[0][1]], X[self.index_start[0][0]:self.index_end[0][0]]
            X[self.index_start[0][2]:self.index_end[0][2]], X[self.index_start[0][3]:self.index_end[0][3]] = \
                X[self.index_start[0][3]:self.index_end[0][3]], X[self.index_start[0][2]:self.index_end[0][2]]
        # big swap -- r mod 4 == 2
        elif r % 4 == 2:
            X[self.index_start[0][0]:self.index_end[0][0]], X[self.index_start[0][2]:self.index_end[0][2]] = \
                X[self.index_start[0][2]:self.index_end[0][2]], X[self.index_start[0][0]:self.index_end[0][0]]
            X[self.index_start[0][1]:self.index_end[0][1]], X[self.index_start[0][3]:self.index_end[0][3]] = \
                X[self.index_start[0][3]:self.index_end[0][3]], X[self.index_start[0][1]:self.index_end[0][1]]
        return X

    def round_const(self, X: list, r: int) -> list:
        """The round consts of gimli. s0,0 ^ 0x9e377900 ^ r.
        Gimli add consts in every fourth round.

        Args:
            X (list): 384 bits input state
            r (int): the round number indicating which round we are running

        Returns:
            list: 384 bits output state
        """
        # every fourth round
        if r % 4 == 0:
            # convert s0_0 to sage binary vector
            s0_0 = vector(self.ring, X[self.index_start[0][0]:self.index_end[0][0]])
            # convert r to sage binary vector
            r = hex2vector(r, self.z)
            # add const
            s0_0 = s0_0 + self.const + r
            # convert to list and return
            X[self.index_start[0][0]:self.index_end[0][0]] = s0_0.list()
        return X

    def round(self, X: list, rounds: int) -> list:
        """round function of gimli

        Args:
            X (list): 384 bits input state
            rounds (int): the number of rounds
        
        Returns:
            list: 384 bits output state
        """
        # a round contains 3 parts
        for r in range(rounds):
            # non-linear layer
            X = self.non_linear(X)
            # linear layer
            X = self.linear_mixing(X, 24 - r)
            # add const
            X = self.round_const(X, 24 - r)
        return X

    def round_without_const(self, X: list, rounds: int) -> list:
        """round function without add const

        Args:
            X (list): 384 bits input state
            rounds (int): the number of rounds
        
        Returns:
            list: 384 bits output state
        """
        # a round contains 2 parts
        for r in range(rounds):
            # non-linear layer
            X = self.non_linear(X)
            # linear layer
            X = self.linear_mixing(X, 24 - r)
        return X

def check_round() -> bool:
    """check if a round function is valid

    Returns:
        bool: whether round function is valid
    """
    # get a logger
    logger = get_logger("check_round")
    logger.info("start")
    # get a Gimli object on GF(2)
    gimli = Gimli(GF(2))
    # use default 24 rounds and the input is all zeroes
    X = [0] * Gimli.state
    X = gimli.round(X, 24)
    # set the truth result(run ./test to see the result)
    ground_truth = "0 0 1 0 0 0 1 1 0 0 0 1 1 0 1 1 1 1 1 0 0 1 1 0 0 0 1 0 0 1 1 0 \
1 1 0 1 1 1 0 0 0 0 0 1 1 1 1 1 0 0 1 1 1 0 1 1 1 1 1 0 0 0 0 0 \
0 0 1 0 1 0 1 1 0 0 0 0 1 1 0 1 1 1 0 1 0 0 0 0 1 1 0 1 1 1 0 0 \
0 0 1 1 0 0 1 0 0 1 1 0 1 1 0 0 1 0 0 0 0 1 0 0 1 1 0 1 1 0 0 0 \
0 0 1 1 1 0 1 1 1 0 0 0 1 1 0 0 0 0 1 0 1 1 0 0 0 0 0 1 0 0 0 0 \
0 1 1 1 0 0 0 1 0 1 1 1 1 1 0 1 1 1 0 1 1 1 1 1 0 1 1 1 0 0 0 0 \
0 0 1 0 0 0 0 1 0 0 0 1 0 1 1 1 0 0 1 0 1 0 1 0 0 0 0 0 0 0 0 0 \
1 0 1 0 1 0 1 0 1 0 0 1 1 0 1 1 1 1 0 1 0 0 0 1 0 0 1 0 0 1 1 0 \
0 1 1 1 0 1 0 0 0 0 1 0 1 1 0 1 1 0 1 1 1 0 1 0 0 1 0 1 0 0 1 0 \
1 1 0 1 0 0 1 1 1 0 0 0 0 0 1 0 0 1 1 0 0 0 0 0 0 1 0 1 0 0 1 1 \
0 1 0 0 0 0 1 1 0 1 0 0 1 0 1 1 1 1 0 0 1 1 1 0 0 1 1 0 0 0 0 1 \
1 0 0 1 0 0 0 0 0 0 0 1 1 0 1 1 0 0 0 0 1 1 0 0 0 1 1 1 0 1 0 0"
    # convert the truth result to list
    ground_truth = [int(i) for i in ground_truth.split()]  # default split by space ' '
    # check if round output is equal to truth result
    res = (X == ground_truth)
    logger.info("check result is " + str(res))
    return res

def check_differential(rounds: int, MAX_ROUNDS: int = 6) -> Any:
    """check if the differential trail is valid(rounds <= MAX_ROUNDS)

    Args:
        rounds (int): the number of rounds
        MAX_ROUNDS (int): the max number of rounds. Defaults to 6
    
    Returns:
        Any: the sat solution
    """
    # get a logger
    logger = get_logger("attack")
    logger.info("start attack")
    # check rounds value
    if rounds > MAX_ROUNDS:
        logger.error("only check rounds <= {}!".format(MAX_ROUNDS))
        exit(1)
    # define ring variables
    # x = a * u + b
    a_vars = ['NVa' + str(i) for i in range(MAX_ROUNDS)]
    b_vars = ['NVb' + str(i) for i in range(MAX_ROUNDS)]
    input_var = 'x'
    block_vars = [input_var] + a_vars + b_vars  # a block size is a gimli state
    auxiliary_var = 'u'
    # defined a polynomial ring
    R = declare_ring([Block(v, Gimli.state) for v in block_vars] + [auxiliary_var], globals())
    # usage of R: R('u') or R(u), input can be a string, we use string here
    # get a Gimli object on R
    gimli = Gimli(R)
    # difference value parameters
    logger.info("start converting hex difference values to binary list")
    # consts in the paper
    c = [0xff898081, 0x80618880, 0x81ff8980, 0x42668080, 0xc0400000, 0x00011100, 
        0x80010080, 0x00402000, 0x80400080, 0x00000080, 0x00400000, 0x80000000]
    # convert hex values to binary list
    for i in range(len(c)):
        c[i] = hex2vector(c[i], Gimli.z).list()
    logger.info("start setting round difference values")
    # all rounds differential trails, each difference value is a state, initial is all zeroes
    diff = [[0] * Gimli.state for i in range(MAX_ROUNDS + 1)]
    # set difference values
    diff[0][gimli.index_start[0][1]:gimli.index_end[0][1]] = diff[0][gimli.index_start[0][3]:gimli.index_end[0][3]] = c[0]
    diff[1][gimli.index_start[1][1]:gimli.index_end[1][1]] = diff[1][gimli.index_start[1][3]:gimli.index_end[1][3]] = c[1]
    diff[1][gimli.index_start[2][1]:gimli.index_end[2][1]] = diff[1][gimli.index_start[2][3]:gimli.index_end[2][3]] = c[2]
    diff[2][gimli.index_start[0][1]:gimli.index_end[0][1]] = diff[2][gimli.index_start[0][3]:gimli.index_end[0][3]] = c[3]
    diff[2][gimli.index_start[1][1]:gimli.index_end[1][1]] = diff[2][gimli.index_start[1][3]:gimli.index_end[1][3]] = c[4]
    diff[2][gimli.index_start[2][1]:gimli.index_end[2][1]] = diff[2][gimli.index_start[2][3]:gimli.index_end[2][3]] = c[5]
    diff[3][gimli.index_start[0][1]:gimli.index_end[0][1]] = diff[3][gimli.index_start[0][3]:gimli.index_end[0][3]] = c[6]
    diff[3][gimli.index_start[1][1]:gimli.index_end[1][1]] = diff[3][gimli.index_start[1][3]:gimli.index_end[1][3]] = c[7]
    diff[3][gimli.index_start[2][1]:gimli.index_end[2][1]] = diff[3][gimli.index_start[2][3]:gimli.index_end[2][3]] = c[8]
    diff[4][gimli.index_start[0][1]:gimli.index_end[0][1]] = diff[4][gimli.index_start[0][3]:gimli.index_end[0][3]] = c[9]
    diff[4][gimli.index_start[1][1]:gimli.index_end[1][1]] = diff[4][gimli.index_start[1][3]:gimli.index_end[1][3]] = c[10]
    diff[4][gimli.index_start[2][1]:gimli.index_end[2][1]] = diff[4][gimli.index_start[2][3]:gimli.index_end[2][3]] = c[11]
    diff[5][gimli.index_start[2][1]:gimli.index_end[2][1]] = diff[5][gimli.index_start[2][3]:gimli.index_end[2][3]] = c[11]
    diff[6][gimli.index_start[0][1]:gimli.index_end[0][1]] = diff[6][gimli.index_start[0][3]:gimli.index_end[0][3]] = c[11]
    # input X
    X = [R(input_var + "({})".format(i)) for i in range(Gimli.state)]
    # set of SAT clauses
    Q = set()
    # add initial differential
    logger.info("add initial differential")
    for i in range(Gimli.state):
        X[i] += diff[0][i] * R(auxiliary_var)
    ##########################################################################################################################
    logger.info("start adding round difference to clauses")
    for r in range(rounds):
        # a round
        # non-linear layer
        X = gimli.non_linear(X)
        # linear layer
        X = gimli.linear_mixing(X, 24 - r)
        # add const
        X = gimli.round_const(X, 24 - r)
        # variable subsitution - x = a * u + b
        for i in range(Gimli.state):
            a = X[i] / R(auxiliary_var)
            b = X[i] + a * R(auxiliary_var)
            # the r th round, the i th variable
            nva_r_i = a_vars[r] + "({})".format(i)
            nvb_r_i = b_vars[r] + "({})".format(i)
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
            else:
                d = X[i] / R(auxiliary_var)
                if d == 0:
                    # meaning it is True, no need to add to Q
                    pass
                elif d == 1:
                    logger.warning("Impossible")
                    exit(0)
                else:
                    Q.add(X[i] / R(auxiliary_var))
        logger.info("end adding round {}".format(r + 1))
    ##########################################################################################################################
    logger.info("start solving")
    res = solve_sat(list(Q))
    logger.info("end solving, result is: ")
    logger.info(res)
    return res

if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(prog="Gimli")
    arg_parser.description = "Check if the differential trail is valid."
    arg_parser.add_argument('-r', '--rounds', type=int, default=6,
                            help="int, please input the number of differential rounds")
    args = arg_parser.parse_args()
    # check_round()
    check_differential(args.rounds)