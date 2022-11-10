import argparse
from read_dc import read_sol

state_x = 4
state_y = 3
state_z = 32

def index(x: int, y: int, z: int) -> int:
    """compute the index of (x,y,z)

    Args:
        x (int): x
        y (int): y
        z (int): z

    Returns:
        int: index
    """
    return state_z*state_x*y + state_z*x + z

def chi(x: int, xbits: int) -> int:
    """chi transformation

    Args:
        x (int): input
        xbits (int): bits length of x

    Returns:
        int: output
    """
    x_bin = []
    y = 0
    # chi transform
    for i in range(xbits):
        x_bin.append((x>>i) & 1)
    for i in range(xbits):
        y += (x_bin[i] ^ ((x_bin[(i+1)%xbits]^1) & x_bin[(i+2)%xbits])) * pow(2, i)
    return y

def get_ddt(xbits: int) -> dict:
    """generate a ddt table of chi

    Args:
        xbits (int): bits length of chi

    Returns:
        dict: ddt table
    """
    x_ = pow(2, xbits)
    ddt = {}
    # delta_a is the input difference
    for delta_a in range(x_):
        tmp = []
        for ax in range(x_):
            # input message ax and ay
            ay = ax ^ delta_a
            # delta_b is the output difference
            delta_b = chi(ax, xbits) ^ chi(ay, xbits)
            # only log the output difference that has not been logged before
            if delta_b not in tmp:
                tmp.append(delta_b)
        ddt[delta_a] = sorted(tmp)
    return ddt

def print_state(X: list) -> None:
    """print a state in column form

    Args:
        X (list): input state
    """
    # print 4*32 columns
    for x in range(state_x):
        lane_print = ""
        for z in range(state_z):
            # now start convert binary column to int
            # get the binary column
            col = 0
            for y in range(state_y):
                col += X[index(x,y,z)] * pow(2, y)
            lane_print += str(col) if col else "."
        print(lane_print)
    print("")


class Extend:
    """extend a trail forward and backward through sbox using dfs
    """
    def __init__(self, trail: list, max_weight: int) -> None:
        # get the ddt table
        self.ddt = get_ddt(state_y)
        # init the trail to extend
        self.trail = trail
        # max weight to extend
        self.max_weight = max_weight
        # set the active columns of the trail's a1 and br-1
        self.set_active_cols()

    def set_active_cols(self) -> None:
        """set the active columns of the trail's a1 and br-1

        """
        # a_1
        a1 = self.trail[0]
        # b_{r-1}
        br_1 = self.trail[-1]

        # get all active columns
        a1_active = []
        br_1_active = []
        for x in range(state_x):
            for z in range(state_z):
                # now start convert binary column to int
                # get the binary column
                col_a = col_b = 0
                for y in range(state_y):
                    col_a += a1[index(x,y,z)] * pow(2, y)
                    col_b += br_1[index(x,y,z)] * pow(2, y)
                if col_a:
                    # use (x,z,value) to represent a column
                    a1_active.append((x,z,col_a))
                if col_b:
                    br_1_active.append((x,z,col_b))
        # set the active columns
        self.active_cols = a1_active + br_1_active
        # set the number of active columns of a1
        self.a1_weight = len(a1_active)

    def dfs(self, i: int, cur_col: list) -> None:
        """dfs self.active_cols

        Args:
            i (int): dfs current index
            cur_col (list): dfs current output
        """
        if i == self.max_weight:
            self.print_extend(cur_col)
            return
        # use (x,z,value) to represent a column
        for delta in self.ddt[self.active_cols[i][2]]:
            tmp_col = (self.active_cols[i][0],self.active_cols[i][1],delta)
            cur_col.append(tmp_col)
            self.dfs(i+1, cur_col)
            cur_col.pop(i)
    
    def extend_by_col(self, state: list, cols: list) -> list:
        """change state with the correspond column changes

        Args:
            state (list): the state to change
            cols (list): the correspond column changes

        Returns:
            list: _description_
        """
        res = state.copy()
        for col in cols:
            # use (x,z,value) to represent a column
            x = col[0]
            z = col[1]
            value = col[2]
            # set column(x,z) to the value of the correspond column changes
            for y in range(state_y):
                res[index(x,y,z)] = (value>>y) & 0x1
        return res
    
    def print_extend(self, col: list) -> None:
        """print the extended trail

        Args:
            col (list): backward and forward column changes (through sbox)
        """
        # print b0, a1, a2, a3
        # if max weight is less than the weight of a1
        if self.a1_weight >= self.max_weight:
            # col0 is the backward column changes
            col0 = col + self.active_cols[self.max_weight:self.a1_weight]
            # colr is the forward column changes
            colr = self.active_cols[self.a1_weight:]
        else:
            col0 = col[:self.a1_weight]
            colr = col[self.a1_weight:self.max_weight] + self.active_cols[self.max_weight:]
        
        print("find one:")
        # backward to b0
        b0 = self.extend_by_col(self.trail[0], col0)
        print("b0:")
        print_state(b0)
        for i in range(len(self.trail) - 1):
            print("a{}:".format(i + 1))
            print_state(self.trail[i])
        # forward to ar
        ar = self.extend_by_col(self.trail[-1], colr)
        print("a3:")
        print_state(ar)


def extend_trail(trails: list, max_weight: int) -> None:
    """extend trail by backward and forward (through sbox)

    Args:
        trails (list): trails to extend
        max_weight (int): max weight to extend, i.e. get 4^(max-weight) extended trails
    """
    for trail in trails:
        ext = Extend(trail, max_weight)
        cur_col = []
        ext.dfs(0, cur_col)


if __name__ == '__main__':
    parse = argparse.ArgumentParser(description="extend solutions")
    parse.add_argument("-f", "--file", type=str, 
                    default="/home/user/anonymous/bosphorus/xoodoo/3rxoodoo/3rxoomess/dc.txt", help="solution file")
    parse.add_argument("-r", "--rounds", type=int, default=3, help="rounds of the solution")
    parse.add_argument("-m", "--max-weight", type=int, default=7, help="max weights to extend, i.e. get 4^(max-weight) extended trails")
    args = parse.parse_args()

    diffs = read_sol(args.file, args.rounds)
    extend_trail(diffs, args.max_weight)