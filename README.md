# Simplex Algorithms
This repository contains several C++ programs utilizing the Simplex algorithm to solve different types of liner optimization problems.

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
minimize -3w + x + 3y - z subject to
w + 2x - y + z = 0,
2w - 2x + 3y + 3z = 9,
w - x + 2y - z = 6,
w >= 0, x >= 0, y >= 0, z >= 0
```
This program can be solved with a minimum of 7 at w = 1, x = 1, y = 3, z = 0.
### Running Simplex
The Simplex program can be compiled with the following command:
```
g++ Simplex.cpp -o Simplex
```
It can then be run with the command
```
./Simplex
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
-3 1 3 -1
1 2 -1 1 = 0
2 -2 3 3 = 9
1 -1 2 -1 = 6
```
The program assumes all variables are greater than or eqaul to 0, so this does not need to be specified. Simplex will compute the desired minimum or maximum and report the result in "out.txt" in the same directory. If the problem is unbounded, it will also be reported in "out.txt".
