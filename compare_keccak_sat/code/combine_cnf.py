n=4
# open file
f1 = open("/home/n2107349e/SAT/keccak/sim_keccak_linear.cnf", "r")
f2 = open("/home/n2107349e/SAT/keccak/firstmessage_chi.cnf", "r")
f3 = open("/home/n2107349e/SAT/keccak/secondmessage_chi.cnf", "r")
# file to write
f4 = open("/home/n2107349e/SAT/keccak/keccak_bos.cnf", "w")
# buf to write
buf = ""
# read file
buf = f1.read()
buf += f2.read()
buf += f3.read()
# find the end of the first line
index_1st = buf.find("\n")
# modify first line
line_1st = buf[:index_1st + 1].split()
line_1st[3] = str(int(line_1st[3]) + (n - 1) * 2 * 320 * 29)
buf = " ".join(line_1st) + "\n" + buf[index_1st + 1:]
f4.write(buf)
# close
f1.close()
f2.close()
f3.close()
f4.close()