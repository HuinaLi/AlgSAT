# AlgSAT
The  source codes abd results are used to help verify the results in our paper.


code files:
1. Keccak
code: the code to describe model initialization and generate the final SAT model in indirect encoding way for verifing differential trails of Keccak.
result: output all ANF and CNF files of each verified trails as well as print the final feasible solution (i.e. a right message pair) and run time.

2. ascon
code: the code to describe model initialization and generate the final SAT model in indirect encoding way for verifing differential trails of Ascon.
result: output all ANF and CNF files of each verified trails as well as print the final feasible solution (i.e. a right message pair) and run time.

3. compare_keccak_sat
code: the code to describe model initialization and generate the final SAT model in direct encoding way for verifing differential trails of Keccak.
result: output all CNF files of each verified trails as well as print the final feasible solution (i.e. a right message pair) and run time.

4. gimli
code: the code to describe model initialization and generate the final SAT model in indirect encoding way for verifing differential trails of Gimli.
result: output all ANF and CNF files of each verified trails as well as print the final feasible solution (i.e. a right message pair) and run time.

Note: A brief user's guide with instructions on how to use Algsat is available in "USER_GUIDE.md" file.

