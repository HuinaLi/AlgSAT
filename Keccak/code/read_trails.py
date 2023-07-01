import argparse

start_sign = "β0"
zero_sign = "-"
z_len = 64
x_len = y_len = 5

def active_hex(h: str) -> list:
    """
    return active bits of a hex(4 bits) value
    """
    # convert to int(base 16)
    h = int(h, 16)
    ret = []
    for i in range(4):
        if (h>>i) & 1:
            ret.append(i)
    return ret

def read_trails(path: str, ROUNDS: int) -> list:
    """
    return the active bits of the trails

    trail starts with "β0"
    """
    f = open(path, "r")
    # read in lines
    lines = f.readlines()
    trails = []
    i = 0
    while i < len(lines):
        if lines[i].find(start_sign) == -1:
            i += 1
            continue
        # find the start of a trail
        else:
            i += 1
            trail = []
            # round
            r = 0
            ### read a state
            while i < len(lines) and r < (ROUNDS+1):
                # y axis
                y = 0
                # only return trail active bits
                state_active_bits = []
                ### read a plane
                while i < len(lines) and y < y_len:
                    # a line is a plane(5 lanes)
                    # exclude '\n'
                    line = lines[i].rstrip("\n")
                    if len(line) < z_len/4 * x_len:
                        i += 1
                        continue
                    # split by "|"
                    plane = line.split("|")
                    if len(plane) != x_len:
                        raise ValueError("a plane at line {} should contains {} lanes, not {} lanes".\
                                            format(i+1, x_len, len(plane)))
                    for x in range(x_len):
                        # read each lane in plane
                        lane = plane[x].replace(" ", "")
                        if len(lane) != z_len//4:
                            raise ValueError("a lane at line {} should contains {} bits, not {} bits".\
                                                format(i+1, z_len, len(lane)*4))
                        # read lane in hex(4bits)
                        for h in range(z_len//4):
                            hbits = lane[h]
                            # not zero (i.e. active)
                            if hbits != zero_sign:
                                # get active bits index
                                for tmpz in active_hex(hbits):
                                    # transform to z index
                                    z = tmpz + (15-h)*4
                                    state_active_bits.append((x, y, z))
                    # next plane
                    i += 1
                    y += 1
                    if y == y_len:
                        trail.append(state_active_bits)
                r += 1
                if r == ROUNDS+1:
                    trails.append(trail)
    return trails

if __name__ == "__main__":
    parse = argparse.ArgumentParser(description="read trails")
    parse.add_argument("-f", "--file", type=str, default='/home/user/lhn/bosphorus/keccak/trails_800.txt', help="file path")
    parse.add_argument("-r", "--rounds", type=int, default=4, help="file path")
    args = parse.parse_args()
    filepath = args.file
    ROUNDS = args.rounds

    # print(read_trails(filepath, ROUNDS))