#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>
#include <stdlib.h>
using namespace std;

void output(vector<vector<float> > *a, vector<float> *b, vector<float> *c, float opt, int numVars)
{
  ofstream outfile;
  outfile.open("out.txt");
  if (outfile.is_open())
  {
    outfile << "Final Tableau:" << endl;
    for (int i=0; i<(*a).size(); i++)
    {
      for (int j=0; j<(*a)[0].size(); j++)
      {
        outfile << (*a)[i][j] << " ";
      }
      outfile << "| " << (*b)[i] << endl;
    }
    for (int i=0; i<(*c).size(); i++)
    {
      outfile << (*c)[i] << " ";
    }
    outfile << "| " << opt << endl;
    outfile << endl;

    outfile << "z* = " << opt << endl;
    vector<float> x;
    for (int i=0; i<numVars; i++)
    {
      bool basic = true;
      int basic_i = 0;
      for (int j=0; j<(*a).size(); j++)
      {
        if ((*a)[j][i] != 0 && (*a)[j][i] != 1)
        {
          basic = false;
        }
        if ((*a)[j][i] == 1)
          basic_i = j;
      }

      if (basic)
        x.push_back((*b)[basic_i]);
      else
        x.push_back(0.0);
    }

    outfile << "x* = (";
    for (int i=0; i<x.size(); i++)
    {
      if (i == x.size()-1)
        outfile << x[i] << ")" << endl;
      else
        outfile << x[i] << ", ";
    }

    outfile.close();
  }
  else
    cout << "Could not open output file." << endl;
}

void outputUnbounded(vector<vector<float> > *a, vector<float> *b, vector<float> *c, bool minimize)
{
  ofstream outfile;
  outfile.open("out.txt");
  if (outfile.is_open())
  {
    outfile << "Final Tableau:" << endl;
    for (int i=0; i<(*a).size(); i++)
    {
      for (int j=0; j<(*a)[0].size(); j++)
      {
        outfile << (*a)[i][j] << " ";
      }
      outfile << "| " << (*b)[i] << endl;
    }
    for (int i=0; i<(*c).size(); i++)
    {
      outfile << (*c)[i] << " ";
    }
    outfile << endl;

    if (minimize)
      outfile << "z* = -infinity" << endl;
    else
      outfile << "z* = infinity" << endl;
    outfile << "The program is unbounded." << endl;

    outfile.close();
  }
  else
    cout << "Could not open output file." << endl;
}

void reduce(vector<vector<float> > *a, vector<float> *b, vector<float> *c, float *opt, int row, int col)
{
  float pivdiv = 1 / (*a)[row][col];
  for (int i=0; i<(*a)[0].size(); i++)
  {
    (*a)[row][i] = (*a)[row][i] * pivdiv;
  }
  (*b)[row] = (*b)[row] * pivdiv;
  for (int i=0; i<(*a).size(); i++)
  {
    if (i == row)
      continue;

    float coeff = (*a)[i][col] * (*a)[row][col];
    for (int j=0; j<(*a)[0].size(); j++)
    {
      (*a)[i][j] = (*a)[i][j] - ((*a)[row][j] * coeff);
    }
    (*b)[i] = (*b)[i] - ((*b)[row] * coeff);
  }
  float coeff = (*c)[col] * (*a)[row][col];
  for (int i=0; i<(*c).size(); i++)
  {
    (*c)[i] = (*c)[i] - ((*a)[row][i] * coeff);
  }
  (*opt) = (*opt) - ((*b)[row] * coeff);
}

void reduceC(vector<vector<float> > *a, vector<float> *b, vector<float> *c, float *opt)
{
  for (int i=0; i<(*a).size(); i++)
  {
    float coeff = 1;
    for (int j=0; j<(*a)[0].size(); j++)
    {
      if ((*a)[i][j] == 1)
      {
        bool basic = true;
        for (int k=0; k<(*a).size(); k++)
        {
          if ((*a)[k][j] != 0 && k != i)
            basic = false;
        }
        if (basic)
          coeff = (*c)[j] * (*a)[i][j];
      }
      (*c)[j] = (*c)[j] - ((*a)[i][j] * coeff);
    }
    (*opt) = (*opt) - ((*b)[i] * coeff);
  }
}

void simplex(vector<vector<float> > *a, vector<float> *b, vector<float> *c, float initialOpt, int numVars, int phase, bool minimize)
{
  // Phase 1
  if (phase == 1)
  {
    vector<float> p1c;
    for (int i=0; i<(*a)[0].size(); i++)
    {
      p1c.push_back(0.0);
    }
    for (int i=0; i<(*a).size(); i++)
    {
      p1c.push_back(1.0);
      for (int j=0; j<(*a).size(); j++)
      {
        if (j == i)
          (*a)[i].push_back(1.0);
        else
          (*a)[i].push_back(0.0);
      }
    }

    float opt = 0;
    reduceC(a, b, &p1c, &opt);

    // for (int i=0; i<(*a).size(); i++)
    // {
    //   for (int j=0; j<(*a)[0].size(); j++)
    //   {
    //     cout << (*a)[i][j] << " ";
    //   }
    //   cout << "| " << (*b)[i] << endl;
    // }
    // for (int i=0; i<p1c.size(); i++)
    // {
    //   cout << p1c[i] << " ";
    // }
    // cout << "| " << opt << endl;
    // cout << endl;

    bool stop = false;
    while (!stop)
    {
      float min_r = 0;
      int piv_col = 0;
      for (int i=0; i<p1c.size(); i++)
      {
        if (p1c[i] < min_r)
        {
          min_r = p1c[i];
          piv_col = i;
        }
      }

      if (min_r >= 0)
      {
        stop = true;
        continue;
      }

      float min_ratio = 0;
      int piv_row = 0;
      for (int i=0; i<(*a).size(); i++)
      {
        if ((*a)[i][piv_col] > 0 && ((*b)[i] / (*a)[i][piv_col] < min_ratio || min_ratio == 0))
        {
          min_ratio = (*b)[i] / (*a)[i][piv_col];
          piv_row = i;
        }
      }

      reduce(a, b, &p1c, &opt, piv_row, piv_col);
    }

    for (int i=0; i<(*a).size(); i++)
    {
      (*a)[i].resize(numVars);
    }

    opt = 0;
    reduceC(a, b, c, &opt);
    simplex(a, b, c, opt, numVars, 2, minimize);
  }
  // Phase 2
  else
  {
    float opt = initialOpt;
    bool stop = false;
    while (!stop)
    {
      float min_r = 0;
      int piv_col = 0;
      for (int i=0; i<(*c).size(); i++)
      {
        if ((*c)[i] < min_r)
        {
          min_r = (*c)[i];
          piv_col = i;
        }
      }

      if (min_r >= 0)
      {
        stop = true;
        if (minimize)
          output(a, b, c, -1*opt, numVars);
        else
          output(a, b, c, opt, numVars);
        continue;
      }

      float min_ratio = 0;
      int piv_row = 0;
      for (int i=0; i<(*a).size(); i++)
      {
        if ((*a)[i][piv_col] > 0 && ((*b)[i] / (*a)[i][piv_col] < min_ratio || min_ratio == 0))
        {
          min_ratio = (*b)[i] / (*a)[i][piv_col];
          piv_row = i;
        }
      }

      if (min_ratio <= 0)
      {
        stop = true;
        outputUnbounded(a, b, c, minimize);
        continue;
      }

      reduce(a, b, c, &opt, piv_row, piv_col);
    }
  }
}

int main()
{
  string line;
  ifstream infile;
  infile.open("in.txt");
  if (infile.is_open())
  {
    int numConstraints;
    if (getline(infile, line))
      numConstraints = atoi(line.c_str());

    int numVars;
    if (getline(infile, line))
      numVars = atoi(line.c_str());

    bool minimize;
    if (getline(infile, line))
    {
      if (line.compare(0, 3, "min") == 0)
        minimize = true;
      else
        minimize = false;
    }

    int c[numVars];
    if (getline(infile, line))
    {
      int i = 0;
      char *token = strtok((char *)line.c_str(), " ");
      while (token)
      {
        c[i] = atoi(token);
        if (!minimize)
          c[i] = -1 * c[i];
        i++;
        token = strtok(NULL, " ");
      }
    }

    vector<vector<float> > vec_a(numConstraints);
    vector<float> vec_b;
    int row = 0;
    int numLT = 0;
    while (getline(infile, line))
    {
      bool past = false;
      char *token = strtok((char *)line.c_str(), " ");
      while (token)
      {
        if (strcmp(token, "=") == 0)
          past = true;
        else if (strcmp(token, "<=") == 0)
        {
          for (int i=0; i<numConstraints; i++)
          {
            if (row == i)
              vec_a[row].push_back(1.0);
            else
              vec_a[row].push_back(0.0);
          }
          numLT++;
          past = true;
        }
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
        else
        {
          if (past)
            vec_b.push_back(atof(token));
          else
            vec_a[row].push_back(atof(token));
        }
        token = strtok(NULL, " ");
      }
      row++;
    }

    infile.close();

    vector<float> vec_c;
    for (int i=0; i<vec_a[0].size(); i++)
    {
      if (i<sizeof(c)/sizeof(c[0]))
        vec_c.push_back(c[i]);
      else
        vec_c.push_back(0.0);
    }

    int phase = 1;
    if (numLT == numConstraints)
      phase = 2;

    simplex(&vec_a, &vec_b, &vec_c, 0, numVars, phase, minimize);

  }
  else
    cout << "Could not open input file." << endl;

  return 0;
}
