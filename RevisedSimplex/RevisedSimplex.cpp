/*! \file RevisedSimplex.cpp
    \Author Brian Burns
    \brief Performs the 2-Phase Revised Simplex Algorithm
*/
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <cstring>
#include <cstdlib>
using namespace std;

//! Outputs the result of the Revised Simplex algorithm to "out.txt"
/*!
  \param a the A matrix
  \param b the b matrix
  \param x the optimal x* vector
  \param bvi the list of basic variable indices
  \param initialBVI the original list of basic variable indices
  \param opt the optimal value found
  \sa outputUnbounded(), outputInfeasible()
*/
void output(vector<vector<float> > *a, vector<float> *b, vector<float> x, vector<int> *bvi, vector<int> initialBVI, float opt)
{
  // open out.txt
  ofstream outfile;
  outfile.open("out.txt");
  // make sure it's open
  if (outfile.is_open())
  {
    // output the final simplex tableau seen
    outfile << "Final Tableau:" << endl;
    outfile << "BVI | B^-1 | b" << endl;
    for (int i=0; i<(*a).size(); i++)
    {
      outfile << (*bvi)[i] << " | ";
      int coli = 0;
      int col = initialBVI[coli];
      for (int j=0; j<(*a)[0].size(); j++)
      {
        if (j == col)
          {
            outfile << (*a)[i][j] << " ";
            coli++;
            col = initialBVI[coli];
          }
      }
      outfile << "| " << (*b)[i] << endl;
    }
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
  \param bvi the list of basic variable indices
  \param initialBVI the original list of basic variable indices
  \param minimize true if min, false if max
  \sa output(), outputInfeasible()
*/
void outputUnbounded(vector<vector<float> > *a, vector<float> *b, vector<int> *bvi, vector<int> initialBVI, bool minimize)
{
  // open out.txt
  ofstream outfile;
  outfile.open("out.txt");
  // make sure it's open
  if (outfile.is_open())
  {
    // output the final simplex tableau seen
    outfile << "Final Tableau:" << endl;
    outfile << "BVI | B^-1 | b" << endl;
    for (int i=0; i<(*a).size(); i++)
    {
      outfile << (*bvi)[i] << " | ";
      int coli = 0;
      int col = initialBVI[coli];
      for (int j=0; j<(*a)[0].size(); j++)
      {
        if (j == col)
          {
            outfile << (*a)[i][j] << " ";
            coli++;
            col = initialBVI[coli];
          }
      }
      outfile << "| " << (*b)[i] << endl;
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
  \param bvi the list of basic variable indices
  \param initialBVI the original list of basic variable indices
  \sa output(), outputUnbounded()
*/
void outputInfeasible(vector<vector<float> > *a, vector<float> *b, vector<int> *bvi, vector<int> initialBVI)
{
  // open out.txt
  ofstream outfile;
  outfile.open("out.txt");
  // make sure it's open
  if (outfile.is_open())
  {
    // output the final simplex tableau seen
    outfile << "Final Tableau:" << endl;
    outfile << "BVI | B^-1 | b" << endl;
    for (int i=0; i<(*a).size(); i++)
    {
      outfile << (*bvi)[i] << " | ";
      int coli = 0;
      int col = initialBVI[coli];
      for (int j=0; j<(*a)[0].size(); j++)
      {
        if (j == col)
          {
            outfile << (*a)[i][j] << " ";
            coli++;
            col = initialBVI[coli];
          }
      }
      outfile << "| " << (*b)[i] << endl;
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
  \param y the y_i vector
  \param bvi the list of basic variable indices
  \param row the pivot row
*/
void reduce(vector<vector<float> > *a, vector<float> *b, vector<float> y, vector<int> bvi, int row)
{
  int coli = 0;
  int col = bvi[coli];
  // get the value needed to make the pivot element 1
  float pivdiv = 1 / y[row];
  // reduce the pivot row with this value
  for (int i=0; i<(*a)[0].size(); i++)
  {
    if (i == col)
    {
      (*a)[row][i] = (*a)[row][i] * pivdiv;
      coli++;
      col = bvi[coli];
    }
  }
  // reduce the right hand side
  (*b)[row] = (*b)[row] * pivdiv;

  coli = 0;
  col = bvi[coli];
  // for every other row
  for (int i=0; i<y.size(); i++)
  {
    if (i == row)
      continue;

    // get the value needed to make the pivot column element 0
    float coeff = y[i];
    // reduce the row with this value
    for (int j=0; j<(*a)[0].size(); j++)
    {
      if (j == col)
      {
        (*a)[i][j] = (*a)[i][j] - ((*a)[row][j] * coeff);
        coli++;
        col = bvi[coli];
      }
    }
    // reduce the right hand side
    (*b)[i] = (*b)[i] - ((*b)[row] * coeff);
    coli = 0;
    col = bvi[coli];
  }
}

//! Calculates c_B * B^(-1)
/*!
  \param c the c matrix
  \param bvi the list of basic variable indices
  \param a the A matrix
  \param initialBVI the original list of basic variable indices
  \return the lambda value lambda = c_B * B^(-1)
*/
vector<float> cBasicB(vector<float> c, vector<int> bvi, vector<vector<float> > a, vector<int> initialBVI)
{
  vector<float> result;
  vector<float> cbasic;

  // get the basic costs
  int col = 0;
  for (int i=bvi[col]; col<bvi.size(); i=bvi[col])
  {
    cbasic.push_back(c[i]);
    col++;
  }

  // multiply the basic cost vector with B^-1
  col = 0;
  for (int i=initialBVI[col]; col<initialBVI.size(); i=initialBVI[col])
  {
    float sum = 0;
    for (int j=0; j<a.size(); j++)
      sum += cbasic[j] * a[j][i];
    result.push_back(sum);
    col++;
  }

  return result;
}

//! Calculates lambda * D
/*!
  \param lambda the lambda vector
  \param a the A matrix
  \param dvi the list of nonbasic variable indices
  \return the value lambda * D
*/
vector<float> lambdaD(vector<float> lambda, vector<vector<float> > a, vector<int> dvi)
{
  vector<float> result;

  // multiply the given lambda by the D matrix
  int col = 0;
  for (int i=dvi[col]; col<dvi.size(); i=dvi[col])
  {
    float sum = 0;
    for (int j=0; j<a.size(); j++)
    {
      sum += lambda[j] * a[j][i];
    }
    result.push_back(sum);
    col++;
  }

  return result;
}

//! Calculates c_D - lambda * D
/*!
  \param v1 the c_D vector
  \param v2 the lambda * D vector
  \param indices the indices of the D matrix within A
  \return the r_D value r_D = c_D - lambda * D
*/
vector<float> cMinusD(vector<float> v1, vector<float> v2, vector<int> indices)
{
  vector<float> result;

  // subtract v2 from v1 (c_D - lambda(D))
  int col = 0;
  int colv2 = 0;
  for (int i=indices[col]; col<indices.size(); i=indices[col])
  {
    float diff = v1[i] - v2[colv2];
    result.push_back(diff);
    col++;
    colv2++;
  }

  return result;
}

//! Calculates y_i = B^(-1) * a_i
/*!
  \param a the A matrix
  \param initialA the original A matrix
  \param bvi the list of basic variable indices
  \param col which column of A to use
  \return the y value y_i = B^(-1) * a_i
*/
vector<float> bInverseA(vector<vector<float> > a, vector<vector<float> > initialA, vector<int> bvi, int col)
{
  // get the column from A that we're using
  vector<float> acol;
  for (int i=0; i<initialA.size(); i++)
    acol.push_back(initialA[i][col]);

  vector<float> result;
  // multiply B^-1 by the A column
  int acoli = 0;
  int bcoli = 0;
  int bcol = bvi[bcoli];
  for (int i=0; i<a.size(); i++)
  {
    float sum = 0;
    for (int j=0; j<a[0].size(); j++)
    {
      // make sure we're using B^-1
      if (j == bcol)
      {
        sum += acol[acoli] * a[i][j];
        acoli++;
        bcoli++;
        bcol = bvi[bcoli];
      }
    }
    result.push_back(sum);
    acoli = 0;
    bcoli = 0;
    bcol = bvi[bcoli];
  }

  return result;
}

//! Performs the Revised Simplex algorithm on a given tableau
/*!
  \param a the A matrix
  \param b the b matrix
  \param c the c matrix
  \param bvi the list of basic variable indices (changes throughout)
  \param initialBVI the original list of basic variable indices (does not change)
  \param dvi the list of nonbasic variable indices (changes throughout)
  \param numVars the number of variables in the program
  \param phase 1 or 2, corresponding to what phase of the Simplex Method we are on
  \param minimize true if min, false if max
*/
void revisedSimplex(vector<vector<float> > *a, vector<float> *b, vector<float> *c, vector<int> *bvi, vector<int> initialBVI, vector<int> *dvi, int numVars, int phase, bool minimize)
{
  // Phase 1
  if (phase == 1)
  {
    // make the phase 1 c row
    vector<float> p1c;
    // make original variables 0
    for (int i=0; i<(*a)[0].size(); i++)
    {
      p1c.push_back(0.0);
    }
    // add as many 1s as rows
    for (int i=0; i<(*a).size(); i++)
    {
      p1c.push_back(1.0);
      // add an identity matrix to A
      for (int j=0; j<(*a).size(); j++)
      {
        if (j == i)
          (*a)[i].push_back(1.0);
        else
          (*a)[i].push_back(0.0);
      }
    }

    vector<vector<float> > initialA = (*a);

    // run the Revised Simplex algorithm
    bool stop = false;
    while (!stop)
    {
      // get lambda
      vector<float> lambda = cBasicB(p1c, (*bvi), (*a), initialBVI);

      // get r_d
      vector<float> rd = cMinusD(p1c, lambdaD(lambda, initialA, (*dvi)), (*dvi));

      // get the pivot column
      float min_rd = 0;
      int piv_coli = 0;
      int piv_col = 0;
      for (int i=0; i<rd.size(); i++)
      {
        // check for the smallest r_d < 0
        if (rd[i] < min_rd)
        {
          min_rd = rd[i];
          piv_coli = i;
          piv_col = (*dvi)[i];
        }
      }

      // check if current solution is optimal
      if (min_rd >= 0)
      {
        stop = true;
        continue;
      }

      // get y_pivcol
      vector<float> y = bInverseA((*a), initialA, initialBVI, piv_col);

      // get the pivot row
      vector<float> ratios;
      // for each row
      for (int i=0; i<y.size(); i++)
      {
        if (y[i] == 0)
          ratios.push_back(-1);
        else
        {
          // denegerate loop ahead, don't pivot here
          if (y[i] <= 0 && (*b)[i] == 0)
            ratios.push_back(-1);
          else
            ratios.push_back((*b)[i] / y[i]);
        }
      }
      float min_ratio = ratios[0];
      int piv_row = 0;
      for (int i=1; i<ratios.size(); i++)
      {
        if (ratios[i] >= 0 && (ratios[i] < min_ratio || min_ratio < 0))
        {
          min_ratio = ratios[i];
          piv_row = i;
        }
      }

      // no positive ratios, problem is unbounded
      if (min_ratio <= 0)
      {
        stop = true;
        outputUnbounded(a, b, bvi, initialBVI, minimize);
        continue;
      }

      // update bvi and dvi
      (*dvi)[piv_coli] = (*bvi)[piv_row];
      (*bvi)[piv_row] = piv_col;

      // update B^-1 and b
      reduce(a, b, y, initialBVI, piv_row);

      // return to step 1 of the Revised Simplex algorithm
    }

    // phase 1 done, cut off the added variables
    vector<int> p2dvi;
    for (int i=0; i<(*dvi).size(); i++)
    {
      if ((*dvi)[i] < numVars)
        p2dvi.push_back((*dvi)[i]);
    }
    // now p2dvi doesn't contain any artificial variables

    // go to phase 2
    revisedSimplex(a, b, c, bvi, initialBVI, &p2dvi, numVars, 2, minimize);
  }
  // Phase 2
  else
  {
    vector<vector<float> > initialA = (*a);

    // run the Revised Simplex algorithm
    bool stop = false;
    while (!stop)
    {
      // get lambda
      vector<float> lambda = cBasicB((*c), (*bvi), (*a), initialBVI);

      // get r_d
      vector<float> rd = cMinusD((*c), lambdaD(lambda, initialA, (*dvi)), (*dvi));

      // get the pivot column
      float min_rd = 0;
      int piv_coli = 0;
      int piv_col = 0;
      for (int i=0; i<rd.size(); i++)
      {
        // check for the smallest r_d < 0
        if (rd[i] < min_rd)
        {
          min_rd = rd[i];
          piv_coli = i;
          piv_col = (*dvi)[i];
        }
      }

      // check if current solution is optimal
      if (min_rd >= 0)
      {
        stop = true;
        // check for infeasibility
        bool infeasible = true;
        for (int i=0; i<rd.size(); i++)
        {
          // any rd > 0, feasible
          if (rd[i] > 0)
          {
            infeasible = false;
            break;
          }
        }

        if (infeasible)
        {
          outputInfeasible(a, b, bvi, initialBVI);
          continue;
        }

        // get x*
        vector<float> x;
        // for each x
        for (int i=0; i<numVars; i++)
        {
          // check if x_i is basic
          bool basic = false;
          int col = 0;
          for (int j=(*bvi)[col]; col<(*bvi).size(); j=(*bvi)[col])
          {
            // basic, x_i = b_j
            if (i == j)
            {
              x.push_back((*b)[col]);
              basic = true;
            }
            col++;
          }
          // not basic, x_i = 0
          if (!basic)
          {
            x.push_back(0.0);
            basic = false;
          }
        }

        // get optimal
        float opt = 0;
        if (minimize)
        {
          // multiply x values with cost function, sum up
          for (int i=0; i<x.size(); i++)
            opt += x[i] * (*c)[i];
        }
        // maximize, flip the cost function
        else
        {
          for (int i=0; i<x.size(); i++)
            opt += x[i] * (-1*(*c)[i]);
        }

        // output the solution
        output(a, b, x, bvi, initialBVI, opt);
        continue;
      }

      // get y_pivcol
      vector<float> y = bInverseA((*a), initialA, initialBVI, piv_col);

      // get the pivot row
      vector<float> ratios;
      // for each row
      for (int i=0; i<y.size(); i++)
      {
        if (y[i] == 0)
          ratios.push_back(-1);
        else
        {
          // denegerate loop ahead, don't pivot here
          if (y[i] <= 0 && (*b)[i] == 0)
            ratios.push_back(-1);
          else
            ratios.push_back((*b)[i] / y[i]);
        }
      }
      float min_ratio = ratios[0];
      int piv_row = 0;
      for (int i=1; i<ratios.size(); i++)
      {
        if (ratios[i] >= 0 && (ratios[i] < min_ratio || min_ratio < 0))
        {
          min_ratio = ratios[i];
          piv_row = i;
        }
      }

      // no positive ratios, problem is unbounded
      if (min_ratio <= 0)
      {
        stop = true;
        outputUnbounded(a, b, bvi, initialBVI, minimize);
        continue;
      }

      // update bvi
      (*dvi)[piv_coli] = (*bvi)[piv_row];
      (*bvi)[piv_row] = piv_col;

      // update B^-1 and b
      reduce(a, b, y, initialBVI, piv_row);

      // return to step 1 of the Revised Simplex algorithm
    }
  }
}

//! Reads input and starts the Revised Simplex procedure
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
    // keep track of less than constraints to potentially skip phase 1
    int numLT = 0;
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
          // keep track of less than constraints
          numLT++;
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
    for (int i=0; i<vec_a[0].size(); i++)
      vec_c.push_back(0.0);

    // initialize basic index storage
    vector<int> vec_bvi;
    vector<int> vec_dvi;

    // assume phase 1
    int phase = 1;
    // check if we made an identity matrix already
    if (numLT == numConstraints)
    {
      phase = 2;
      // push slack variables to basic list, rest to nonbasic list
      for (int i=0; i<vec_a[0].size(); i++)
      {
        if (i>=numVars)
          vec_bvi.push_back(i);
        else
          vec_dvi.push_back(i);
      }
    }
    // phase 1 definitely
    else
    {
      // make basic list artificial variable indices, rest to nonbasic list
      for (int i=0; i<vec_a[0].size()+vec_a.size(); i++)
      {
        if (i<vec_a[0].size())
          vec_dvi.push_back(i);
        else
          vec_bvi.push_back(i);
      }
    }

    // run the Revised Simplex algorithm
    revisedSimplex(&vec_a, &vec_b, &vec_c, &vec_bvi, vec_bvi, &vec_dvi, numVars, phase, minimize);

  }
  // error
  else
    cout << "Could not open input file." << endl;

  return 0;
}
