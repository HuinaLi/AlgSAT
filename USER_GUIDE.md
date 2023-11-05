# AlgSAT_a-SAT-Method-for-Search-and-Verification
The  source codes and results are used to help verify the results in our paper.

## Step 1: Set up environment

# conda setup
1. install miniconda: https://docs.conda.io/en/latest/miniconda.html
```bash Miniconda3-latest-Linux-x86_64.sh```
Note: if you have error "UnicodeDecodeError: ‘ascii‘ codec can‘t decode byte", just type command:
```export LC_ALL=C.UTF-8```
```source ~/.bashrc```
# delete previous environment: </br>
```conda env remove -n name```

# [sagemath install](https://doc.sagemath.org/html/en/installation/conda.html)
1. set conda-forge: </br>```conda config --add channels conda-forge```
2. change to restrict mode: </br>```conda config --set channel_priority strict```
3. create a conda environment named `sage`: </br>```conda create -n sage sage python=3.9```
4. note: you can change the environment name as you like</br>```conda create -n {name} sage python=3.9```


# Install Bosphorus 
The open-source tool is available at https://github.com/meelgroup/bosphorus
1. ```git clone https://github.com/meelgroup/bosphorus.git```

# Install SAT Solver

Choose 1: Install cryptominisat
1. activate your sage environment: </br>```conda activate sage```
2. install cryptominisat: </br>``` conda install cryptominisat``` 
3. To deactivate an active environment, use
```conda deactivate```

cryptominisat usage in linux
for example:
 ```cryptominisat5 6rgimli.cnf -t 20 > 6rgimli_cms20.log 2>&1 &```

Choose 2: Install cadical
1. git clone https://github.com/arminbiere/cadical.git
2. ./configure && make

cadical usage in linux
for example:
 ```cadical  6rgimli.cnf > 6rgimli_cadical.log 2>&1 &```

## Step 2: Generate ANF file (SageMath & Python)
1. Run .py file to generate the original .anf file.
   
We take 6rgimli.anf for example:
```python 6rattack.py > 6rgimli.anf 2>&1 &```

## Step 3: Generate CNF file using Bosphorus
To get final CNFs:
1. Basic command
``` ./build/bosphorus --anfread 6rgimli.anf --anfwrite .6rgimli_out.anf --cnfwrite 6rgimli.cnf ```


## Step 4: Solve our SAT model
 To get a feasible solution

 cryptominisat usage:
for example:
 ```cryptominisat5 6rgimli.cnf -t 20 > 6rgimli_cms20.log 2>&1 &```

## Note
1. We suggest using a script file to automatically finish the above four steps, eg., please refer to our solve.sh file.
2. We  suggest running CryptoMiniSat or CadiCaL sat solvers for at least 5 different, randomly generated, similarly hard problems to get an average time to solve.
