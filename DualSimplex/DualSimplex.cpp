/*! \file DualSimplex.cpp
    \Author Brian Burns
    \brief Performs the Dual Simplex Algorithm
*/
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>
using namespace std;

//! Outputs the result of the Dual Simplex algorithm to "out.txt"
/*!
  \param a the A matrix
  \param b the b matrix
  \param c the c matrix
  \param x the optimal x* vector
  \param opt the optimal value found
  \param numVars the number of variables in the program
  \sa outputUnbounded(), outputInfeasible()
*/
void output(vector<vector<float> > *a, vector<float> *b, vector<float> *c, vector<float> x, float opt, int numVars)
{
  // open out.txt
  ofstream outfile;
  outfile.open("out.txt");
  // make sure it's open
  if (outfile.is_open())
  {
    // output the final simplex tableau seen
    outfile << "Final Tableau:" << endl;
    for (int i=0; i<(*a).size(); i++)
    {
      for (int j=0; j<(*a)[0].size(); j++)
      {
        // output the A matrix
        outfile << (*a)[i][j] << " ";
      }
      // output the b matrix
      outfile << "| " << (*b)[i] << endl;
    }
    // output the c matrix
    for (int i=0; i<(*c).size(); i++)
    {
      outfile << (*c)[i] << " ";
    }
    // output the optimal value found
    outfile << "| " << opt << endl;
    outfile << endl;

    outfile << "z* = " << opt << endl;

    // output the optimal x
    outfile << "x* = (";
    for (int i=0; i<x.size(); i++)
    {
      if (i == x.size()-1)
        outfile << x[i] << ")" << endl;
      else
        outfile << x[i] << ", ";
    }

    // close out.txt
    outfile.close();
  }
  // error
  else
    cout << "Could not open output file." << endl;
}

//! Outputs an unbounded solution to "out.txt"
/*!
  \param a the A matrix
  \param b the b matrix
  \param c the c matrix
  \param minimize true if min, false if max
  \sa output(), outputInfeasible()
*/
void outputUnbounded(vector<vector<float> > *a, vector<float> *b, vector<float> *c, bool minimize)
{
  // open out.txt
  ofstream outfile;
  outfile.open("out.txt");
  // make sure it's open
  if (outfile.is_open())
  {
    // output the final simplex tableau seen
    outfile << "Final Tableau:" << endl;
    for (int i=0; i<(*a).size(); i++)
    {
      for (int j=0; j<(*a)[0].size(); j++)
      {
        // output the A matrix
        outfile << (*a)[i][j] << " ";
      }
      // output the b matrix
      outfile << "| " << (*b)[i] << endl;
    }
    // output the c matrix
    for (int i=0; i<(*c).size(); i++)
    {
      outfile << (*c)[i] << " ";
    }
    outfile << endl;

    // optimal is unbounded
    if (minimize)
      outfile << "z* = -infinity" << endl;
    else
      outfile << "z* = infinity" << endl;
    outfile << "The program is unbounded." << endl;

    // close out.txt
    outfile.close();
  }
  // error
  else
    cout << "Could not open output file." << endl;
}

//! Outputs an infeasible solution to "out.txt"
/*!
  \param a the A matrix
  \param b the b matrix
  \param c the c matrix
  \sa output(), outputUnbounded()
*/
void outputInfeasible(vector<vector<float> > *a, vector<float> *b, vector<float> *c)
{
  // open out.txt
  ofstream outfile;
  outfile.open("out.txt");
  // make sure it's open
  if (outfile.is_open())
  {
    // output the final simplex tableau seen
    outfile << "Final Tableau:" << endl;
    for (int i=0; i<(*a).size(); i++)
    {
      for (int j=0; j<(*a)[0].size(); j++)
      {
        // output the A matrix
        outfile << (*a)[i][j] << " ";
      }
      // output the b matrix
      outfile << "| " << (*b)[i] << endl;
    }
    // output the c matrix
    for (int i=0; i<(*c).size(); i++)
    {
      outfile << (*c)[i] << " ";
    }
    outfile << endl;

    // problem is infeasible
    outfile << "The program is infeasible." << endl;

    // close out.txt
    outfile.close();
  }
  // error
  else
    cout << "Could not open output file." << endl;
}

//! Performs row operations to pivot on an element and reduce a tableau
/*!
  \param a the A matrix
  \param b the b matrix
  \param c the c matrix
  \param opt the optimal value in the current tableau
  \param row the pivot row
  \param col the pivot column
  \sa reduceC()
*/
void reduce(vector<vector<float> > *a, vector<float> *b, vector<float> *c, float *opt, int row, int col)
{
  // get the value needed to make the pivot element 1
  float pivdiv = 1 / (*a)[row][col];
  // reduce the pivot row with this value
  for (int i=0; i<(*a)[0].size(); i++)
  {
    (*a)[row][i] = (*a)[row][i] * pivdiv;
  }
  (*b)[row] = (*b)[row] * pivdiv;

  // for every other row
  for (int i=0; i<(*a).size(); i++)
  {
    if (i == row)
      continue;

    // get the value needed to make the pivot column element 0
    float coeff = (*a)[i][col] * (*a)[row][col];
    // reduce the row with this value
    for (int j=0; j<(*a)[0].size(); j++)
    {
      (*a)[i][j] = (*a)[i][j] - ((*a)[row][j] * coeff);
    }
    (*b)[i] = (*b)[i] - ((*b)[row] * coeff);
  }
  // get the value needed to make the c element 0 in the pivot column
  float coeff = (*c)[col] * (*a)[row][col];
  // reduce the c row with this value
  for (int i=0; i<(*c).size(); i++)
  {
    (*c)[i] = (*c)[i] - ((*a)[row][i] * coeff);
  }
  (*opt) = (*opt) - ((*b)[row] * coeff);
}

//! Performs row operations to make c elements in basic columns 0
/*!
  \param a the A matrix
  \param b the b matrix
  \param c the c matrix
  \param opt the optimal value in the current tableau
  \sa reduce()
*/
void reduceC(vector<vector<float> > *a, vector<float> *b, vector<float> *c, float *opt)
{
  // for each row
  for (int i=0; i<(*a).size(); i++)
  {
    // default the row operation coefficient to 0
    float coeff = 0;
    for (int j=0; j<(*a)[0].size(); j++)
    {
      // found a 1, check for a basic column
      if ((*a)[i][j] == 1)
      {
        bool basic = true;
        for (int k=0; k<(*a).size(); k++)
        {
          // not 0 or 1, not basic
          if ((*a)[k][j] != 0 && k != i)
            basic = false;
        }
        // basic, make sure the c element becomes 0
        if (basic)
          coeff = (*c)[j] * (*a)[i][j];
      }
    }
    // reduce the c vector
    for (int k=0; k<(*a)[0].size(); k++)
    {
      (*c)[k] = (*c)[k] - ((*a)[i][k] * coeff);
    }
    // apply the reduction to the optimal value
    (*opt) = (*opt) - ((*b)[i] * coeff);
  }
}

//! Performs the Dual Simplex algorithm on a given tableau
/*!
  \param a the A matrix
  \param b the b matrix
  \param c the c matrix
  \param initialC the original c matrix (for calculating optimal)
  \param numVars the number of variables in the program
  \param minimize true if min, false if max
*/
void dualSimplex(vector<vector<float> > *a, vector<float> *b, vector<float> *c, vector<float> initialC, int numVars, bool minimize)
{
  float opt = 0;

  // reduce the cost vector for basic variables
  reduceC(a, b, c, &opt);

  // run the Dual Simplex algorithm
  bool stop = false;
  while (!stop)
  {
    // get the pivot row
    float min_b = 0;
    int piv_row = 0;
    for (int i=0; i<(*b).size(); i++)
    {
      // check for the smallest b < 0
      if ((*b)[i] < min_b)
      {
        min_b = (*b)[i];
        piv_row = i;
      }
    }

    // no b < 0, we're done
    if (min_b >= 0)
    {
      stop = true;

      // get x* from the matrix
      vector<float> x;
      // for each column
      for (int i=0; i<numVars; i++)
      {
        // check if x_i is basic
        bool basic = true;
        int basic_i = 0;
        for (int j=0; j<(*a).size(); j++)
        {
          // not 0 or 1, not basic
          if ((*a)[j][i] != 0 && (*a)[j][i] != 1)
          {
            basic = false;
          }
          // get the basic index
          if ((*a)[j][i] == 1)
            basic_i = j;
        }

        // basic x, get the b value
        if (basic)
          x.push_back((*b)[basic_i]);
        // nonbasic, 0
        else
          x.push_back(0.0);
      }

      // get optimal
      opt = 0;
      // multiply x values with cost function, sum up
      for (int i=0; i<x.size(); i++){
        opt += x[i] * initialC[i];
      }
      // min function, just output
      if (minimize)
        output(a, b, c, x, opt, numVars);
      // max function, output with negative optimal
      else
        output(a, b, c, x, -1*opt, numVars);
      continue;
    }

    // get the pivot column
    vector<float> ratios;
    float ratio = 0;
    for (int j=0; j<(*a)[piv_row].size(); j++)
    {
      // don't divide by 0
      if ((*a)[piv_row][j] == 0)
        ratios.push_back(0.0);
      // r_j / -y_j
      else
        ratios.push_back((*c)[j] / (-1 * (*a)[piv_row][j]));
    }

    // find the minimum positive ratio
    float min_ratio = ratios[0];
    int piv_col = 0;
    for (int i=1; i<ratios.size(); i++)
    {
      if (ratios[i] > 0 && (ratios[i] < min_ratio || min_ratio <= 0))
      {
        min_ratio = ratios[i];
        piv_col = i;
      }
    }

    // no positive ratios, problem is infeasible
    if (min_ratio <= 0)
    {
      stop = true;
      outputInfeasible(a, b, c);
      continue;
    }

    // reduce the tableau on the pivot row and column
    reduce(a, b, c, &opt, piv_row, piv_col);

    // return to step 1 of the Simplex algorithm
  }
}

//! Reads input and starts the Dual Simplex procedure
/*!
  \return 0 for a successful run, variable otherwise
*/
int main()
{
  // open "in.txt"
  string line;
  ifstream infile;
  infile.open("in.txt");
  // make sure it's open
  if (infile.is_open())
  {
    // get the number of constraints
    int numConstraints;
    if (getline(infile, line))
      numConstraints = atoi(line.c_str());

    // get the number of variables
    int numVars;
    if (getline(infile, line))
      numVars = atoi(line.c_str());

    // get whether this problem is min or max
    bool minimize;
    if (getline(infile, line))
    {
      if (line.compare(0, 3, "min") == 0)
        minimize = true;
      else
        minimize = false;
    }

    // get the c row
    vector<float> vec_c;
    if (getline(infile, line))
    {
      // read in one at a time
      int i = 0;
      char *token = strtok((char *)line.c_str(), " ");
      while (token)
      {
        if (minimize)
          vec_c.push_back(atof(token));
        else
          // maximize function, negate c
          vec_c.push_back(-1*atof(token));
        i++;
        token = strtok(NULL, " ");
      }
    }

    // get the a and b matrices
    vector<vector<float> > vec_a(numConstraints);
    vector<float> vec_b;
    int row = 0;
    // keep track of greater than constraints to potentially flip constraints
    int numGT = 0;
    while (getline(infile, line))
    {
      // A or b side of the constraint
      bool past = false;
      // read in one at a time
      char *token = strtok((char *)line.c_str(), " ");
      while (token)
      {
        // check for = | <= | >=
        if (strcmp(token, "=") == 0)
          past = true;
        // less than, add an identity row
        else if (strcmp(token, "<=") == 0)
        {
          for (int i=0; i<numConstraints; i++)
          {
            if (row == i)
              vec_a[row].push_back(1.0);
            else
              vec_a[row].push_back(0.0);
          }
          past = true;
        }
        // greater than, add a negative identity row
        else if (strcmp(token, ">=") == 0)
        {
          for (int i=0; i<numConstraints; i++)
          {
            if (row == i)
              vec_a[row].push_back(-1.0);
            else
              vec_a[row].push_back(0.0);
          }
          past = true;
          // keep track of greater than constraints
          numGT++;
        }
        // number
        else
        {
          // add to b matrix
          if (past)
            vec_b.push_back(atof(token));
          // add to A matrix
          else
            vec_a[row].push_back(atof(token));
        }
        token = strtok(NULL, " ");
      }
      // go to the next row
      row++;
    }

    // close in.txt
    infile.close();

    // convert the c array, add 0s for added columns
    for (int i=numVars; i<vec_a[0].size(); i++)
      vec_c.push_back(0.0);

    // check if we need to flip constraints
    if (numGT == numConstraints)
    {
      for (int i=0; i<vec_a.size(); i++)
      {
        for (int j=0; j<vec_a[0].size(); j++)
        {
          // -0 looks tacky
          if (vec_a[i][j] != 0)
            // flip A matrix
            vec_a[i][j] = -1 * vec_a[i][j];
        }
        // flip b matrix
        vec_b[i] = -1 * vec_b[i];
      }
    }

    // run the Dual Simplex algorithm
    dualSimplex(&vec_a, &vec_b, &vec_c, vec_c, numVars, minimize);

  }
  // error
  else
    cout << "Could not open input file." << endl;

  return 0;
}
