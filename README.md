# Simplex Algorithms
This repository contains several C++ programs utilizing the Simplex algorithm to solve different types of linear optimization problems.

## Simplex
The simplex algorithm solves linear optimization problems of the forms:
```
minimize cTx subject to Ax = b, x >= 0
minimize cTx subject to Ax <= b, x >= 0
minimize cTx subject to Ax >= b, x >= 0
```
```
maximize cTx subject to Ax = b, x >= 0
maximize cTx subject to Ax <= b, x >= 0
maximize cTx subject to Ax >= b, x >= 0
```
For example, a problem of this form is:
```
minimize 3w + x + 3y + z subject to
w + 2x - y + z = 3,
5w - 4x + 2y + 3z = 10,
w - 2x + 2y + 3z = 5,
w >= 0, x >= 0, y >= 0, z >= 0
```
This program can be solved with a minimum of 5.5 at w = 1.3333, x = 0.1667, y = 0, z = 1.3333.

### Running Simplex
The Simplex program can be compiled with the following command:

__Unix__
```
g++ Simplex.cpp -o Simplex
```

__Windows__
```
cl Simplex.cpp /o Simplex
```
It can then be run with the command:

__Unix__
```
./Simplex
```

__Windows__
```
.\Simplex
```
Simplex reads the input from the "in.txt" file in the same directory. The format of "in.txt" is as follows:
```
number of constraint equations
number of variables
min/max
objective function
constraint equations
```
The example program above has an "in.txt" file specified as such:
```
3
4
min
3 1 3 1
1 2 -1 1 = 3
5 -4 2 3 = 10
1 -2 2 3 = 5
```
The program assumes all variables are greater than or equal to 0, so this does not need to be specified. Simplex will compute the desired minimum or maximum and report the result in "out.txt" in the same directory. If the problem is unbounded or infeasible, it will also be reported in "out.txt".

## Revised Simplex
The revised simplex algorithm solves linear optimization problems of the forms:
```
minimize cTx subject to Ax = b, x >= 0
minimize cTx subject to Ax <= b, x >= 0
minimize cTx subject to Ax >= b, x >= 0
```
```
maximize cTx subject to Ax = b, x >= 0
maximize cTx subject to Ax <= b, x >= 0
maximize cTx subject to Ax >= b, x >= 0
```
For example, a problem of this form is:
```
minimize 3w + x + 3y + z subject to
w + 2x - y + z = 3,
5w - 4x + 2y + 3z = 10,
w - 2x + 2y + 3z = 5,
w >= 0, x >= 0, y >= 0, z >= 0
```
This program can be solved with a minimum of 5.5 at w = 1.3333, x = 0.1667, y = 0, z = 1.3333.

### Running Revised Simplex
The Revised Simplex program can be compiled with the following command:

__Unix__
```
g++ RevisedSimplex.cpp -o RevisedSimplex
```

__Windows__
```
cl RevisedSimplex.cpp /o RevisedSimplex
```
It can then be run with the command:

__Unix__
```
./RevisedSimplex
```

__Windows__
```
.\RevisedSimplex
```
RevisedSimplex reads the input from the "in.txt" file in the same directory. The format of "in.txt" is as follows:
```
number of constraint equations
number of variables
min/max
objective function
constraint equations
```
The example program above has an "in.txt" file specified as such:
```
3
4
min
3 1 3 1
1 2 -1 1 = 3
5 -4 2 3 = 10
1 -2 2 3 = 5
```
The program assumes all variables are greater than or equal to 0, so this does not need to be specified. RevisedSimplex will compute the desired minimum or maximum and report the result in "out.txt" in the same directory. If the problem is unbounded or infeasible, it will also be reported in "out.txt".

## Dual Simplex
The dual simplex algorithm solves linear optimization problems of the forms:
```
minimize cTx subject to Ax = b, x >= 0
minimize cTx subject to Ax <= b, x >= 0
minimize cTx subject to Ax >= b, x >= 0
```
```
maximize cTx subject to Ax = b, x >= 0
maximize cTx subject to Ax <= b, x >= 0
maximize cTx subject to Ax >= b, x >= 0
```
For example, a problem of this form is:
```
minimize -7u + 7v - 2w - x - 6y subject to
3u - v + w - 2x = -3,
2u + v + x + y = 4,
-u + 3v - 3x + z = 12,
u >= 0, v >= 0, w >= 0, x >= 0, y >= 0, z >= 0
```
This program can be solved with a maximum of -16.5 at u = 0, v = 0, w = 0, x = 1.5, y = 2.5, z = 16.5.
__Note:__ For the dual simplex algorithm, the problem entered must already have a feasible solution available from the beginning (positive or negative identity matrix).

### Running Dual Simplex
The Dual Simplex program can be compiled with the following command:

__Unix__
```
g++ DualSimplex.cpp -o DualSimplex
```

__Windows__
```
cl DualSimplex.cpp /o DualSimplex
```
It can then be run with the command:

__Unix__
```
./DualSimplex
```

__Windows__
```
.\DualSimplex
```
Dual Simplex reads the input from the "in.txt" file in the same directory. The format of "in.txt" is as follows:
```
number of constraint equations
number of variables
min/max
objective function
constraint equations
```
The example program above has an "in.txt" file specified as such:
```
3
6
min
-7 7 -2 -1 -6 0
3 -1 1 -2 0 0 = -3
2 1 0 1 1 0 = 4
-1 3 0 -3 0 1 = 12
```
The program assumes all variables are greater than or equal to 0, so this does not need to be specified. Dual Simplex will compute the desired minimum or maximum and report the result in "out.txt" in the same directory. If the problem is unbounded or infeasible, it will also be reported in "out.txt".

## Primal-Dual Simplex
The primal-dual simplex algorithm solves linear optimization problems of the forms:
```
minimize cTx subject to Ax = b, x >= 0
maximize cTx subject to Ax = b, x >= 0
```

For example, a problem of this form is:
```
minimize 2x + y + 4z subject to
x + y + 2z = 3,
2x + y + 3z = 5,
x >= 0, y >= 0, z >= 0
```
This program can be solved with a minimum of 5 at x = 2, y = 1, z = 0.
__Note:__ For the primal-dual simplex algorithm, the problem entered must have an initial dual feasible solution (lambda) specified in "in.txt".

### Running Primal-Dual Simplex
The Primal-Dual Simplex program can be compiled with the following command:

__Unix__
```
g++ PrimalDualSimplex.cpp -o PrimalDualSimplex
```

__Windows__
```
cl PrimalDualSimplex.cpp /o PrimalDualSimplex
```
It can then be run with the command:

__Unix__
```
./PrimalDualSimplex
```

__Windows__
```
.\PrimalDualSimplex
```
Primal-Dual Simplex reads the input from the "in.txt" file in the same directory. The format of "in.txt" is as follows:
```
number of constraint equations
number of variables
initial lambda
min/max
objective function
constraint equations
```
The example program above has an "in.txt" file specified as such:
```
2
3
0 0
min
2 1 4
1 1 2 = 3
2 1 3 = 5
```
The program assumes all variables are greater than or equal to 0, so this does not need to be specified. Primal-Dual Simplex will compute the desired minimum or maximum and report the result in "out.txt" in the same directory. If the problem is unbounded or infeasible, it will also be reported in "out.txt".
