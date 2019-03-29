# MOwNiT lab 2
Operations on matricies



## Contents
Zestaw funkcji i "bibliotek" implementujących:
* rozwiązywanie układu równań metodą Gaussa-Jordana z pełnym pivotingiem
* faktoryzację LU metodą Choleskiego
* obliczanie prądów w gałęziach przy pomocy metody potencjałów węzłowych
* generowanie obrazów obwodów z skalą RGB oznaczającą poziomy natężeń
  w gałęziach i napięć na węzłach - korzysta z pakietu graphviz


## Prerequisites
Graphviz suite - should be already installed on your distro. If not:
```
$ sudo apt install graphviz
```


## Installing
```
$ cd /path/to/folder/containing/this/readme
$ cmake .
$ make
```


## Running
### GAUSS-JORDAN SYSTEM OF EQUATION SOLVER
For this program graph to be solved is hardcoded. Modify values or use
one of the funciton form libmatrix to initialize matrix with random values
```
$ ./bin/zad1
```

### LU DECOMPOSITION
Here values are also hardcoded. Make necessary adjustments incode yourself.
```
$ ./bin/zad2
```

### BRANCH CURRENT SOLVER
Here you can modify ./src/circuit.h by adding '#define DEBUG 1' bellow libraries includes.
This will enable different color coding scheme
Specify file with information on circuit connections, nodes, resistances and voltage applied to given nodes

file format (nodes are counted from 0!):
```
(int)node_count (int)branch_count
(int)ground_node_id (int)source_node_id (float)voltage_applied
(int)start_node_id (int)end_node_id (float)branch_resistance
(int)start_node_id (int)end_node_id (float)branch_resistance
(int)start_node_id (int)end_node_id (float)branch_resistance
itd...
```
examples of these files are in ./connections
 
code bellow will create circuit image in file_name's directory
```
$ ./bin/zad3 file_name
```

## Contributing
File matrix.cpp is straight copy of file matrix.c
I did it becouse cmake Makefile created by cmake couldn't link c library to cpp file(?)
If anybody came across this problem before please let me know


## Authors
* **Krzysztof Kolasa** - *Man needs to do his homework*
