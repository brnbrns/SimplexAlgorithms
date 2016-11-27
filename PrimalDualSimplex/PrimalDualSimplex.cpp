/*! \file PrimalDualSimplex.cpp
    \Author Brian Burns
    \brief Performs the Primal-Dual Simplex Algorithm
*/
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>
using namespace std;

//! Outputs the result of the Primal-Dual Simplex algorithm to "out.txt"
/*!
  \param a the A matrix
  \param b the b matrix
  \param r the r matrix
  \param cA the c - lambda*A matrix
  \param x the optimal x* vector
  \param opt the optimal value found
  \sa outputUnbounded(), outputInfeasible()
*/
void output(vector<vector<float> > *a, vector<float> *b, vector<float> r, vector<float> cA, vector<float> x, float opt)
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
    // output the r matrix
    for (int i=0; i<r.size(); i++)
    {
      outfile << r[i] << " ";
    }
    // output the optimal value found
    outfile << "| " << 0 << " (r_ARP)" << endl;
    // output the c - lambda * A matrix
    for (int i=0; i<cA.size(); i++)
      outfile << cA[i] << " ";
    outfile << "(c_T - lambda*A)" << endl;
    outfile << endl;

    // output the optimal value found
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
  \param r the r matrix
  \param cA the c - lambda*A matrix
  \param minimize true if min, false if max
  \sa output(), outputInfeasible()
*/
void outputUnbounded(vector<vector<float> > *a, vector<float> *b, vector<float> r, vector<float> cA, bool minimize)
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
    // output the r matrix
    for (int i=0; i<r.size(); i++)
    {
      outfile << r[i] << " ";
    }
    // output the optimal value found
    outfile << " (r_ARP)" << endl;
    // output the c - lambda * A matrix
    for (int i=0; i<cA.size(); i++)
      outfile << cA[i] << " ";
    outfile << "(c_T - lambda*A)" << endl;

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
  \param r the r matrix
  \param cA the c - lambda*A matrix
  \sa output(), outputUnbounded()
*/
void outputInfeasible(vector<vector<float> > *a, vector<float> *b, vector<float> r, vector<float> cA)
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
    // output the r matrix
    for (int i=0; i<r.size(); i++)
    {
      outfile << r[i] << " ";
    }
    // output the optimal value found
    outfile << " (r_ARP)" << endl;
    // output the c - lambda * A matrix
    for (int i=0; i<cA.size(); i++)
      outfile << cA[i] << " ";
    outfile << "(c_T - lambda*A)" << endl;

    // problem is infeasible
    outfile << "The program is infeasible." << endl;

    // close out.txt
    outfile.close();
  }
  // error
  else
    cout << "Could not open output file." << endl;
}

//! Calculates lambda * A
/*!
  \param lambda the lambda vector
  \param a the A matrix
  \param numVars the number of variables in the program
  \return the value lambda * A
*/
vector<float> lambdaA(vector<float> lambda, vector<vector<float> > a, int numVars)
{
  vector<float> result;

  // multiply the given lambda by the A matrix
  for (int i=0; i<numVars; i++)
  {
    // multiply each column
    float sum = 0;
    for (int j=0; j<a.size(); j++)
    {
      sum += lambda[j] * a[j][i];
    }
    result.push_back(sum);
  }

  return result;
}

//! Calculates c - lambda * A
/*!
  \param v1 the c vector
  \param v2 the lambda * A vector
  \return the value c - lambda * A
*/
vector<float> cMinusA(vector<float> v1, vector<float> v2)
{
  vector<float> result;

  // subtract v2 from v1 (c - lambda*A)
  for (int i=0; i<v2.size(); i++)
  {
    float diff = v1[i] - v2[i];
    result.push_back(diff);
  }

  return result;
}

//! Performs row operations to pivot on an element and reduce a tableau
/*!
  \param a the A matrix
  \param b the b matrix
  \param c the c matrix
  \param opt the optimal value in the current tableau
  \param row the pivot row
  \param col the pivot column
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
  // reduce the optimal value
  (*opt) = (*opt) - ((*b)[row] * coeff);
}

//! Performs the Primal-Dual Simplex algorithm on a given tableau
/*!
  \param a the A matrix
  \param b the b matrix
  \param initialC the initial value of the c matrix
  \param lambda the initial lambda value
  \param numVars the number of variables in the program
  \param minimize true if min, false if max
*/
void primalDualSimplex(vector<vector<float> > *a, vector<float> *b, vector<float> initialC, vector<float> *lambda, int numVars, bool minimize)
{
  // start by reducing the tableau to get r and opt
  float opt = 0;
  for (int i=0; i<(*b).size(); i++)
    opt -= (*b)[i];

  // get the r vector
  vector<float> r;
  for (int i=0; i<(*a)[0].size(); i++)
  {
    // nonbasic x, subtract the A values
    if (i<numVars)
    {
      float diff = 0;
      for (int j=0; j<(*a).size(); j++)
        diff -= (*a)[j][i];
      r.push_back(diff);
    }
    // basic y, 0
    else
      r.push_back(0.0);
  }

  // get lambda * A
  vector<float> lamA = lambdaA((*lambda), (*a), numVars);
  // get c - lambda * A
  vector<float> cA = cMinusA(initialC, lamA);

  // run the Primal-Dual Simplex algorithm
  bool stop = false;
  while (!stop)
  {
    // get epsilon
    vector<float> ratios;
    for (int i=0; i<cA.size(); i++)
      // calculate each ratio (c - lambda*A) / -r
      if (r[i] < 0)
        ratios.push_back(cA[i] / (-1 * r[i]));
      // r >= 0, unwanted ratio
      else
      {
        // ensure the ratio won't be picked
        if (i == 0)
          ratios.push_back(1000000);
        else
          ratios.push_back(i*100000);
      }
    // get the minimum ratio
    float epsilon = ratios[0];
    int piv_col = 0;
    for (int i=1; i<ratios.size(); i++)
    {
      if (ratios[i] < epsilon)
      {
        epsilon = ratios[i];
        // save the pivot column
        piv_col = i;
      }
    }

    // reduce ARP
    if (epsilon <= 0)
    {
      // get piv_row
      float col_ratio = (*b)[0] / (*a)[0][piv_col];
      int piv_row = 0;
      for (int i=1; i<(*b).size(); i++)
      {
        // get the minimum positive ratio
        if ((*a)[i][piv_col] > 0 && (*b)[i] / (*a)[i][piv_col] < col_ratio)
        {
          col_ratio = (*b)[i] / (*a)[i][piv_col];
          // save the pivot row
          piv_row = i;
        }
      }

      // reduce
      reduce(a, b, &r, &opt, piv_row, piv_col);
    }
    else
    {
      // check if y = 0
      bool ynonbasic = true;
      // for each column
      for (int i=numVars; i<(*a)[0].size(); i++)
      {
        // check if y_i is basic
        bool basic = true;
        for (int j=0; j<(*a).size(); j++)
        {
          // not 0 or 1, not basic
          if ((*a)[j][i] != 0 && (*a)[j][i] != 1)
          {
            basic = false;
          }
        }
        // basic y, y != 0
        if (basic)
          ynonbasic = false;
      }

      // y = 0, we're done
      if (ynonbasic)
      {
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
            {
              basic_i = j;
              // make sure other column entries are 0
              for (int k=0; k<(*a).size(); k++)
              {
                if (j == k)
                  continue;
                if ((*a)[k][i] != 0)
                  basic = false;
              }
            }
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

        // output
        output(a, b, r, cA, x, opt);
        stop = true;
        continue;
      }
      else
      {
        // check for infeasibility
        bool feasible = false;
        for (int i=0; i<r.size(); i++)
        {
          // any r[i] < 0, feasible
          if (r[i] < 0)
            feasible = true;
        }

        if (!feasible)
        {
          outputInfeasible(a, b, r, cA);
          stop = true;
          continue;
        }

        // get new c - lambda * A
        for (int i=0; i<r.size(); i++)
          cA[i] = cA[i] + (epsilon * r[i]);
      }
    }
  }
}

//! Reads input and starts the Primal-Dual Simplex procedure
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

    // get the initial lambda
    vector<float> vec_lambda;
    if (getline(infile, line))
    {
      // read in one at a time
      char *token = strtok((char *)line.c_str(), " ");
      while (token)
      {
        vec_lambda.push_back(atof(token));
        token = strtok(NULL, " ");
      }
    }

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
      char *token = strtok((char *)line.c_str(), " ");
      while (token)
      {
        if (minimize)
          vec_c.push_back(atof(token));
        else
          // maximize function, negate c
          vec_c.push_back(-1*atof(token));
        token = strtok(NULL, " ");
      }
    }

    // get the a and b matrices
    vector<vector<float> > vec_a(numConstraints);
    vector<float> vec_b;
    int row = 0;
    while (getline(infile, line))
    {
      // A or b side of the constraint
      bool past = false;
      // read in one at a time
      char *token = strtok((char *)line.c_str(), " ");
      while (token)
      {
        // check for =
        if (strcmp(token, "=") == 0)
        {
          past = true;
          // add an identity row for dual problem
          for (int i=0; i<numConstraints; i++)
          {
            if (row == i)
              vec_a[row].push_back(1.0);
            else
              vec_a[row].push_back(0.0);
          }
          past = true;
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

    // run the Dual Simplex algorithm
    primalDualSimplex(&vec_a, &vec_b, vec_c, &vec_lambda, numVars, minimize);

  }
  // error
  else
    cout << "Could not open input file." << endl;

  return 0;
}
