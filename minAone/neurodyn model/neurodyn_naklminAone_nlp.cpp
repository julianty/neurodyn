// neurodyn_nakl.cpp
// Nonlinear Ipopt program

// Author: Bryan A. Toth
// btoth@physics.ucsd.edu

#include "neurodyn_naklminAone_nlp.hpp"
#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <cstring>
#ifndef HAVE_CSTDIO
#define HAVE_CSTDIO
# include <cstdio>
# include <iostream>
# include "myfunctions.hpp"
# include <fstream>
# include <string>
# include <stdlib.h>
#else
# ifndef HAVE_STDIO_H
#define HAVE_STDIO_H
#  include <stdio.h>
# else
#  error "don't have header file for stdio"
# endif
#endif

using namespace Ipopt;

using namespace std;

double Alphaspline(double, double);
double Alphasplinejac(double, double, int);
double Alphasplinehes(double, double, int, int);
double Betaspline(double, double);
double Betasplinejac(double, double, int);
double Betasplinehes(double, double, int, int);
// constructor
NEURODYN_NAKL_NLP::NEURODYN_NAKL_NLP(int id)
{
  nU=0;
  nP=6;
  nY=4;
  nM=4;
  nI=1;

  K11val = new double[nU];
  K11val2 = new double[nU];
  K11valp1 = new double[nU];
  dK11val = new double[nU];
  dK11val2 = new double[nU];
  dK11valp1 = new double[nU];
  Xdval = new double[nM];
  Xdval2 = new double[nM];
  Xdvalp1 = new double[nM];
  Xval = new double[nY];
  Xval2 = new double[nY];
  Xvalp1 = new double[nY];
  Pval = new double[nP];
  Ival = new double[nI];
  Ival2 = new double[nI];
  Ivalp1 = new double[nI];
  Rf0 = new double[nY];

  string buffer;
  specs = new string[9+nP+nY+nI+nU+nM+1+1];

  int count;
  count = 0;
  
  ifstream fin ("specs.txt");
  if (fin.is_open()){
    while (! fin.eof()){
      getline (fin,buffer);
      if (buffer[0] !='#'){
        specs[count] = buffer;
        count++;
      }
    }
    fin.close();
  }
  else cout << "Unable to open file";

  Time = atoi(specs[0].c_str());
  skip = atoi(specs[1].c_str());
  hstep = atof(specs[2].c_str());

  string filename;
  int ret;
  double tempdata;
  Ntotal = (2*Time+1)*nY +nP;
  solution = new double[Ntotal];
  int toggle = 0;
  file_fmt = atoi(specs[3].c_str());

  // Read in observation files

    if (file_fmt == 1) toggle = 1;
    if (file_fmt == 2) toggle = 2;
    if (file_fmt == 3) toggle = 3;

    if (file_fmt < 2) {
      // Read in files one by one if file_fmt = 0,1

      taskid = id;

      VDATA0 = new double[2*Time+1];
      VDATA0dummy = new double[skip];
      FILE *pFile0;
      filename = specs[4+toggle];
      pFile0 = fopen(filename.c_str(),"r");
      for(Index jt=0;jt<skip;jt++)
      {
        ret = fscanf (pFile0, "%lf", &VDATA0dummy[jt]);
        if (ret == EOF) break;
      }
      for(Index jt=0;jt<2*Time+1;jt++)
      {
        ret = fscanf (pFile0, "%lf", &VDATA0[jt]);
        if (ret == EOF) break;
      }
      fclose (pFile0);

      MDATA0 = new double[2*Time+1];
      MDATA0dummy = new double[skip];
      FILE *pFile1;
      filename = specs[5+toggle];
      pFile1 = fopen(filename.c_str(),"r");
      for(Index jt=0;jt<skip;jt++)
      {
        ret = fscanf (pFile1, "%lf", &MDATA0dummy[jt]);
        if (ret == EOF) break;
      }
      for(Index jt=0;jt<2*Time+1;jt++)
      {
        ret = fscanf (pFile1, "%lf", &MDATA0[jt]);
        if (ret == EOF) break;
      }
      fclose (pFile1);

      HDATA0 = new double[2*Time+1];
      HDATA0dummy = new double[skip];
      FILE *pFile2;
      filename = specs[6+toggle];
      pFile2 = fopen(filename.c_str(),"r");
      for(Index jt=0;jt<skip;jt++)
      {
        ret = fscanf (pFile2, "%lf", &HDATA0dummy[jt]);
        if (ret == EOF) break;
      }
      for(Index jt=0;jt<2*Time+1;jt++)
      {
        ret = fscanf (pFile2, "%lf", &HDATA0[jt]);
        if (ret == EOF) break;
      }
      fclose (pFile2);

      NDATA0 = new double[2*Time+1];
      NDATA0dummy = new double[skip];
      FILE *pFile3;
      filename = specs[7+toggle];
      pFile3 = fopen(filename.c_str(),"r");
      for(Index jt=0;jt<skip;jt++)
      {
        ret = fscanf (pFile3, "%lf", &NDATA0dummy[jt]);
        if (ret == EOF) break;
      }
      for(Index jt=0;jt<2*Time+1;jt++)
      {
        ret = fscanf (pFile3, "%lf", &NDATA0[jt]);
        if (ret == EOF) break;
      }
      fclose (pFile3);

    }
    else if (file_fmt >= 2){
      // Read in files from single file if file_fmt = 2,3

      VDATA0 = new double[2*Time+1];
      VDATA0dummy = new double[skip];
      MDATA0 = new double[2*Time+1];
      MDATA0dummy = new double[skip];
      HDATA0 = new double[2*Time+1];
      HDATA0dummy = new double[skip];
      NDATA0 = new double[2*Time+1];
      NDATA0dummy = new double[skip];

      FILE *pFile0;
      filename = specs[4+toggle];
      char init_idx[21]; 
      modid = atoi(specs[4].c_str());
      if (modid != 0){
        pathid = id / modid;
        taskid = id % modid;
        sprintf(init_idx, "%d.%s", pathid, specs[5].c_str());
      }
      else{
        pathid = 0;
        taskid = id;
        sprintf(init_idx, ".%s", specs[5].c_str());
      }
      filename += init_idx;
      pFile0 = fopen(filename.c_str(),"r");

      for(Index rows=0;rows<skip;rows++){
        for (Index cols=0;cols<nY;cols++){
          if (cols == 0){
            ret = fscanf (pFile0, "%lf", &VDATA0dummy[rows]);
          }
          else if (cols == 1){
            ret = fscanf (pFile0, "%lf", &MDATA0dummy[rows]);
          }
          else if (cols == 2){
            ret = fscanf (pFile0, "%lf", &HDATA0dummy[rows]);
          }
          else if (cols == 3){
            ret = fscanf (pFile0, "%lf", &NDATA0dummy[rows]);
          }
          else {
            ret = fscanf (pFile0, "%lf", &tempdata);
          }          
           if (ret == EOF) break;
        }
      }
      for(Index rows=0;rows<2*Time+1;rows++){
        for (Index cols=0;cols<nY;cols++){
          if (cols == 0){
            ret = fscanf (pFile0, "%lf", &VDATA0[rows]);
          }
          else if (cols == 1){
            ret = fscanf (pFile0, "%lf", &MDATA0[rows]);
          }
          else if (cols == 2){
            ret = fscanf (pFile0, "%lf", &HDATA0[rows]);
          }
          else if (cols == 3){
            ret = fscanf (pFile0, "%lf", &NDATA0[rows]);
          }
          else {
            ret = fscanf (pFile0, "%lf", &tempdata);
          }          
           if (ret == EOF) break;
        }
      }
      fclose (pFile0);
    }

  // Read in stimulus files

    if (file_fmt < 2) {
      // Read in files one by one if file_fmt = 0,1

      Iinj = new double[2*Time+1];
      Iinjdummy = new double[skip];
      FILE *qFile0;
      filename = specs[8+toggle];
      qFile0 = fopen(filename.c_str(),"r");
      for(Index jt=0;jt<skip;jt++)
      {
        ret = fscanf (qFile0, "%lf", &Iinjdummy[jt]);
        if (ret == EOF) break;
      }
      for(Index jt=0;jt<2*Time+1;jt++)
      {
        ret = fscanf (qFile0, "%lf", &Iinj[jt]);
        if (ret == EOF) break;
      }
      fclose (qFile0);

    }

    else if (file_fmt >= 2){
      // Read in data from a single file if file_fmt = 2,3

      Iinj = new double[2*Time+1];
      Iinjdummy = new double[skip];

      FILE *qFile0;
      filename = specs[5+toggle];
      char init_idx[21]; 
      modid = atoi(specs[4].c_str());
      if (modid != 0){
        pathid = id / modid;
        taskid = id % modid;
        sprintf(init_idx, "%d.%s", pathid, specs[5].c_str());
      }
      else{
        pathid = 0;
        taskid = id;
        sprintf(init_idx, ".%s", specs[5].c_str());
      }
      filename += init_idx;
      qFile0 = fopen(filename.c_str(),"r");

      for(Index rows=0;rows<skip;rows++){
        for (Index cols=0;cols<nI;cols++){
          if (cols == 0){
            ret = fscanf (qFile0, "%lf", &Iinjdummy[rows]);
          }
        }
      }
      for(Index rows=0;rows<2*Time+1;rows++){
        for (Index cols=0;cols<nI;cols++){
          if (cols == 0){
            ret = fscanf (qFile0, "%lf", &Iinj[rows]);
          }
        }
      }
      fclose (qFile0);
    }

  //Read in boundary conditions

  if (file_fmt == 0) toggle = nM + nI;
  if (file_fmt == 1) toggle = nM + nI + 1;
  if (file_fmt == 2) toggle = (nM > 0) + (nI > 0) + 2;
  if (file_fmt == 3) toggle = (nM > 0) + (nI > 0) + 3;

  int rows = nY+nU+nP+1;
  bounds = new double*[rows];
  for (Index i=0;i<rows;i++) bounds[i] = new double[4];
  int counter;
  for(Index k=0;k<rows;k++){
    counter=0;
    char* tmp = new char[specs[4+toggle+k].size()+1];
    strcpy( tmp, specs[4+toggle+k].c_str() );
    char *ptr = strtok(tmp,",");
    bounds[k][3] = 0.0;
    while(ptr != 0){
      if(counter<3){
        bounds[k][counter] = atof(ptr);
      }
      if(counter==3) {
        bounds[k][counter] = atof(ptr);
      }
      ptr = strtok(0,",");
      counter++;
    }
  }
  for (Index i=0;i<nY;i++){
    Rf0[i]=bounds[i][2];
    bounds[i][3]=bounds[i][2];
  }

  beta=0;
  alpha = bounds[nY+nU+nP][0];
  delta_beta=(int) bounds[nY+nU+nP][1];
  max_beta=(int) bounds[nY+nU+nP][2];

  //Read in output_fmt (zero is default)
  output_fmt = atoi(specs[4+toggle+rows].c_str());
}

// destructor
NEURODYN_NAKL_NLP::~NEURODYN_NAKL_NLP()
{
  delete [] K11val;
  delete [] K11val2;
  delete [] K11valp1;
  delete [] dK11val;
  delete [] dK11val2;
  delete [] dK11valp1;
  delete [] Xdval;
  delete [] Xdval2;
  delete [] Xdvalp1;
  delete [] Xval;
  delete [] Xval2;
  delete [] Xvalp1;
  delete [] Pval;
  delete [] Ival;
  delete [] Ival2;
  delete [] Ivalp1;
  delete [] specs;
  delete [] VDATA0;
  delete [] VDATA0dummy;
  delete [] MDATA0;
  delete [] MDATA0dummy;
  delete [] HDATA0;
  delete [] HDATA0dummy;
  delete [] NDATA0;
  delete [] NDATA0dummy;
  delete [] Iinj;
  delete [] Iinjdummy;
  int rows = nY+nU+nP;
  for (Index i=0;i<rows;i++) delete [] bounds[i];
  delete [] bounds;
}

bool NEURODYN_NAKL_NLP::changeRf(){
  if((beta+delta_beta)>(max_beta-1))
    return false;
  else
    beta = beta + delta_beta;
  for (Index i=0;i<nY;i++) {
    bounds[i][3]=pow(alpha,beta)*Rf0[i];
  }
  printf("\n\n\n\n<------       Beta=%d       ------->\n\n\n\n",beta);
  return true;
}

// returns the size of the problem
bool NEURODYN_NAKL_NLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
        Index& nnz_h_lag, IndexStyleEnum& index_style)

{
  // Number of variables
  n = 8*Time+10;

  // Number of equality constraints
  m = 0;

  // Number of Jacobian nonzero entries
  nnz_jac_g = 0;

  // Number of Hessian nonzero entries
  nnz_h_lag = 34*(Time+1)+82*Time+21;

  // use the C style indexing (0-based)
  index_style = TNLP::C_STYLE;

  return true;
}


// returns the variable bounds
bool NEURODYN_NAKL_NLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                Index m, Number* g_l, Number* g_u)
{
  // Here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
  assert(n == 8*Time+10);
  assert(m == 0);

  for(Index jt=0;jt<Time+1;jt++) {
    for(Index var=0;var<nY;var++) {
      // Bounds for x
      x_l[(Time+1)*var+jt]=bounds[var][0];
      x_u[(Time+1)*var+jt]=bounds[var][1];
      // Bounds for midpoints
      if(jt<Time){
        x_l[(Time+1)*(nY+nU)+Time*var+jt]=bounds[var][0];
        x_u[(Time+1)*(nY+nU)+Time*var+jt]=bounds[var][1];
      }
    }
    for(Index cup=0;cup<nU;cup++) {
      // Bounds for k
      x_l[(Time+1)*(nY+cup)+jt]=bounds[nY+cup][0];
      x_u[(Time+1)*(nY+cup)+jt]=bounds[nY+cup][1];
      // Bounds for midpoints
      if(jt<Time) {
        x_l[(Time+1)*(nY+nU)+Time*(nY+cup)+jt]=bounds[nY+cup][0];
        x_u[(Time+1)*(nY+nU)+Time*(nY+cup)+jt]=bounds[nY+cup][1];
      }
    }
  } 

  for(Index par=0;par<nP;par++) {
    // Bounds for parameters
    x_l[2*Time*(nY+nU)+nY+nU+par]=bounds[nY+nU+par][0];
    x_u[2*Time*(nY+nU)+nY+nU+par]=bounds[nY+nU+par][1];
  }

  return true;
}
// returns the initial point for the problem
bool NEURODYN_NAKL_NLP::get_starting_point(Index n, bool init_x, Number* x,
                                   bool init_z, Number* z_L, Number* z_U,
                                   Index m, bool init_lambda,
                                   Number* lambda)
{
  assert(init_x == true);
  assert(init_z == false);
  assert(init_lambda == false);

  for (Index i=0; i<n; i++) {
    x[i] = 0.0;
  }

  if(beta==0){
    double **skipinit = new double* [skip];
    double **init = new double* [2*Time + 1];
    double *param_init = new double [nP];

    for(Index i=0;i<skip;i++) skipinit[i] = new double[nY];
    for(Index i=0;i<2*Time+1;i++) init[i] = new double[nY];

    string filename;
    int ret;

    // Read in initial data files for file_fmt = 1,3

    if (file_fmt == 1){
      filename = specs[4];
      FILE *initFILE;
      initFILE = fopen(filename.c_str(),"r");
      for(Index jt=0;jt<skip;jt++){
        for (Index jy=0;jy<nY;jy++){
          ret = fscanf (initFILE,"%lf",&skipinit[jt][jy]);
          if (ret == EOF) break;
        }
      }
      for(Index jt=0;jt<(2*Time+1);jt++){
        for (Index jy=0;jy<nY;jy++){
          ret = fscanf (initFILE,"%lf",&init[jt][jy]);
          if (ret == EOF) break;
        }
      }
      for(Index i=0;i<nP;i++){
        ret = fscanf (initFILE,"%lf",&param_init[i]);
        if (ret == EOF) break;
      }
      fclose (initFILE);
    }
    else if (file_fmt == 3){
      filename = specs[6];
      FILE *initFILE;
      char init_idx[21];
      if (modid != 0){
        sprintf(init_idx, "%d.%s", taskid, specs[5].c_str());
      }
      else{
        sprintf(init_idx, ".%s", specs[5].c_str());
      }
      filename += init_idx;
      initFILE = fopen(filename.c_str(),"r");

      for(Index jt=0;jt<skip;jt++){
        for (Index jy=0;jy<nY;jy++){
          ret = fscanf (initFILE,"%lf",&skipinit[jt][jy]);
          if (ret == EOF) break;
        }
      }
      for(Index jt=0;jt<(2*Time+1);jt++){
        for (Index jy=0;jy<nY;jy++){
          ret = fscanf (initFILE,"%lf",&init[jt][jy]);
          if (ret == EOF) break;
        }
      }
      for(Index i=0;i<nP;i++){
        ret = fscanf (initFILE,"%lf",&param_init[i]);
        if (ret == EOF) break;
      }
      fclose (initFILE);
    }

    // Save init_data to variables for file_fmt = 0,1,2,3

    if ((file_fmt == 1) || (file_fmt == 3)){
      for (Index jt=0;jt<Time+1;jt++){
        for (Index var=0;var<nY;var++){

        // Initial conditions for x and midpoints
          x[(Time+1)*var+jt] = init[2*jt][var];
          if (jt<Time) x[(Time+1)*(nY+nU)+Time*var+jt] = init[2*jt+1][var];
        }

        // Initial conditions for controls and midpoints
        for (Index cup=0;cup<nU;cup++){
          x[(Time+1)*(cup+nY)+jt]=bounds[cup+nY][2];
          if (jt<Time) x[(Time+1)*(nY+nU)+Time*(cup+nY)+jt]=bounds[cup+nY][2];
        }
      }

      // Initial conditions for parameters
      for (Index par=0;par<nP;par++){
        x[(2*Time+1)*(nY+nU)+par] = param_init[par];
      }
    }
    else if ((file_fmt == 0) || (file_fmt == 2)){
      srand(taskid);
      for (Index jt=0;jt<Time+1;jt++){
        for(Index var=0;var<nY;var++){
          x[(Time+1)*var+jt] = rand()*1.0/RAND_MAX*(bounds[var][1]-bounds[var][0])+bounds[var][0];
          if (jt<Time) x[(Time+1)*(nY+nU)+Time*var+jt] = rand()*1.0/RAND_MAX*(bounds[var][1]-bounds[var][0])+bounds[var][0];
        }
        for (Index cup=0;cup<nU;cup++){
          x[(Time+1)*(cup+nY)+jt]=bounds[cup+nY][2];
          if (jt<Time) x[(Time+1)*(nY+nU)+Time*(cup+nY)+jt]=bounds[cup+nY][2];
        }
        for (Index par=0;par<nP;par++){
          x[(2*Time+1)*(nY+nU)+par] = rand()*1.0/RAND_MAX*(bounds[nY+nU+par][1]-bounds[nY+nU+par][0])+bounds[nY+nU+par][0];
        }
      }
    }
    for(Index i=0;i<2*Time+1;i++) delete [] init[i];
    delete [] init;
  }
  else{
    for(Index i=0;i<Ntotal-nP;i++){
      x[i] = solution[i];
    }
    for (Index jt=0;jt<Time+1;jt++){
      for (Index cup=0;cup<nU;cup++){
        x[(Time+1)*(cup+nY)+jt]=bounds[cup+nY][2];
        if (jt<Time) x[(Time+1)*(nY+nU)+Time*(cup+nY)+jt]=bounds[cup+nY][2];
      }
    }
    for(Index i=0;i<nP;i++){
      x[(2*Time+1)*(nU+nY)+i] = solution[Ntotal-nP+i];
    }
  }
  return true;
}

// returns the value of the objective function
bool NEURODYN_NAKL_NLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  assert(n == 8*Time+10);
  obj_value = 0;

  for(Index jt=0;jt<Time;jt++) {

     for(Index i=0;i<nY;i++) {
        Xval[i] = x[i*(Time+1) + jt];
        Xvalp1[i] = x[i*(Time+1) + jt + 1];
        Xval2[i] = x[(Time+1)*(nY+nU) + i*(Time) + jt];
     } //end for loop


     for(Index cup=0;cup<nU;cup++) {
        K11val[cup] = x[nY*(Time+1) + cup*(Time+1) + jt];
        K11valp1[cup] = x[nY*(Time+1) + cup*(Time+1) + jt + 1];
        K11val2[cup] = x[(Time+1)*(nY+nU) + Time*nY + cup*(Time) + jt];
     } //end for loop

     Xdval[0] = VDATA0[2*jt];
     Xdval2[0] = VDATA0[2*jt+1];
     Xdvalp1[0] = VDATA0[2*jt+2];
     Xdval[1] = MDATA0[2*jt];
     Xdval2[1] = MDATA0[2*jt+1];
     Xdvalp1[1] = MDATA0[2*jt+2];
     Xdval[2] = HDATA0[2*jt];
     Xdval2[2] = HDATA0[2*jt+1];
     Xdvalp1[2] = HDATA0[2*jt+2];
     Xdval[3] = NDATA0[2*jt];
     Xdval2[3] = NDATA0[2*jt+1];
     Xdvalp1[3] = NDATA0[2*jt+2];
     Ival[0] = Iinj[2*jt];
     Ival2[0] = Iinj[2*jt+1];
     Ivalp1[0] = Iinj[2*jt+2];

     for(Index i=0;i<nP;i++) {
        Pval[i] = x[(2*Time+1)*(nY+nU)+i];
     } //end for loop


    obj_value += pow(Xdval[2] - Xval[2], 2) + pow(Xdval2[2] - Xval2[2], 2) + pow(Xdval[1] - Xval[1], 2) + pow(Xdval2[1] - Xval2[1], 2) + pow(Xdval[3] - Xval[3], 2) + pow(Xdval2[3] - Xval2[3], 2) + pow(Xdval[0] - Xval[0], 2) + pow(Xdval2[0] - Xval2[0], 2); 

    obj_value += bounds[0][3]*(pow(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0], 2) + pow(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0], 2)); 

    obj_value += bounds[1][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
// Betaspline
pow(0.166666666666667*hstep*(-Xval[1]*Betaspline(Xval[0], 1) - 4*Xval2[1]*Betaspline(Xval2[0], 1) - Xvalp1[1]*Betaspline(Xvalp1[0], 1) + (-Xval[1] + 1)*Alphaspline(Xval[0], 1) + (-4*Xval2[1] + 4)*Alphaspline(Xval2[0], 1) + (-Xvalp1[1] + 1)*Alphaspline(Xvalp1[0], 1)) + Xval[1] - Xvalp1[1], 2) + pow(0.125*hstep*(-Xval[1]*Betaspline(Xval[0], 1) + Xvalp1[1]*Betaspline(Xvalp1[0], 1) + (-Xval[1] + 1)*Alphaspline(Xval[0], 1) - (-Xvalp1[1] + 1)*Alphaspline(Xvalp1[0], 1)) + 0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1], 2)); 

    obj_value += bounds[2][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
// Betaspline
pow(Xval[2] - Xvalp1[2] + 0.166666666666667*hstep*(-Xval[2]*Betaspline(Xval[0], 2) - 4*Xval2[2]*Betaspline(Xval2[0], 2) - Xvalp1[2]*Betaspline(Xvalp1[0], 2) + (-Xval[2] + 1)*Alphaspline(Xval[0], 2) + (-4*Xval2[2] + 4)*Alphaspline(Xval2[0], 2) + (-Xvalp1[2] + 1)*Alphaspline(Xvalp1[0], 2)), 2) + pow(0.5*Xval[2] - Xval2[2] + 0.5*Xvalp1[2] + 0.125*hstep*(-Xval[2]*Betaspline(Xval[0], 2) + Xvalp1[2]*Betaspline(Xvalp1[0], 2) + (-Xval[2] + 1)*Alphaspline(Xval[0], 2) - (-Xvalp1[2] + 1)*Alphaspline(Xvalp1[0], 2)), 2)); 

    obj_value += bounds[3][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
// Betaspline
pow(0.166666666666667*hstep*(-Xval[3]*Betaspline(Xval[0], 3) - 4*Xval2[3]*Betaspline(Xval2[0], 3) - Xvalp1[3]*Betaspline(Xvalp1[0], 3) + (-Xval[3] + 1)*Alphaspline(Xval[0], 3) + (-4*Xval2[3] + 4)*Alphaspline(Xval2[0], 3) + (-Xvalp1[3] + 1)*Alphaspline(Xvalp1[0], 3)) + Xval[3] - Xvalp1[3], 2) + pow(0.125*hstep*(-Xval[3]*Betaspline(Xval[0], 3) + Xvalp1[3]*Betaspline(Xvalp1[0], 3) + (-Xval[3] + 1)*Alphaspline(Xval[0], 3) - (-Xvalp1[3] + 1)*Alphaspline(Xvalp1[0], 3)) + 0.5*Xval[3] - Xval2[3] + 0.5*Xvalp1[3], 2)); 

  } //end for loop

// Add last element
     for(Index i=0;i<nY;i++) {
        Xval[i] = x[Time + i*(Time+1)];
        Xvalp1[i] = 0;
        Xval2[i] = 0;
     } //end for loop


     for(Index cup=0;cup<nU;cup++) {
        K11val[cup] = x[nY*(Time+1) + Time + cup*(Time+1)];
        K11valp1[cup] = 0;
        K11val2[cup] = 0;
     } //end for loop

     Xdval[0] = VDATA0[2*Time];
     Xdval2[0] = 0;
     Xdvalp1[0] = 0;
     Xdval[1] = MDATA0[2*Time];
     Xdval2[1] = 0;
     Xdvalp1[1] = 0;
     Xdval[2] = HDATA0[2*Time];
     Xdval2[2] = 0;
     Xdvalp1[2] = 0;
     Xdval[3] = NDATA0[2*Time];
     Xdval2[3] = 0;
     Xdvalp1[3] = 0;
     Ival[0] = Iinj[2*Time];
     Ival2[0] = 0;
     Ivalp1[0] = 0;

     for(Index i=0;i<nP;i++) {
        Pval[i] = x[(2*Time+1)*(nY+nU)+i];
     } //end for loop


  obj_value += pow(Xdval[2] - Xval[2], 2) + pow(Xdval2[2] - Xval2[2], 2) + pow(Xdval[1] - Xval[1], 2) + pow(Xdval2[1] - Xval2[1], 2) + pow(Xdval[3] - Xval[3], 2) + pow(Xdval2[3] - Xval2[3], 2) + pow(Xdval[0] - Xval[0], 2) + pow(Xdval2[0] - Xval2[0], 2);

  obj_value = obj_value/(2*Time+1);

  return true;
}


// return the gradient of the objective function grad_{x} f(x)
bool NEURODYN_NAKL_NLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  assert(n == 8*Time+10);

  for(Index i=0;i<n;i++) {
     grad_f[i] = 0;
  }

  for(Index jt=0;jt<Time;jt++) {

     for(Index i=0;i<nY;i++) {
        Xval[i] = x[jt + i*(Time+1)];
        Xvalp1[i] = x[jt + i*(Time+1) + 1];
        Xval2[i] = x[(Time+1)*(nY+nU) + i*(Time)+ jt];
     } //end for loop

     for(Index cup=0;cup<nU;cup++) {
        K11val[cup] = x[nY*(Time+1) + cup*(Time+1) + jt];
        K11valp1[cup] = x[nY*(Time+1) + cup*(Time+1) + jt + 1];
        K11val2[cup] = x[(Time+1)*(nY+nU) + (nY)*Time + cup*(Time) +  jt];
     } //end for loop

     Xdval[0] = VDATA0[2*jt];
     Xdval2[0] = VDATA0[2*jt+1];
     Xdvalp1[0] = VDATA0[2*jt+2];
     Xdval[1] = MDATA0[2*jt];
     Xdval2[1] = MDATA0[2*jt+1];
     Xdvalp1[1] = MDATA0[2*jt+2];
     Xdval[2] = HDATA0[2*jt];
     Xdval2[2] = HDATA0[2*jt+1];
     Xdvalp1[2] = HDATA0[2*jt+2];
     Xdval[3] = NDATA0[2*jt];
     Xdval2[3] = NDATA0[2*jt+1];
     Xdvalp1[3] = NDATA0[2*jt+2];
     Ival[0] = Iinj[2*jt];
     Ival2[0] = Iinj[2*jt+1];
     Ivalp1[0] = Iinj[2*jt+2];

     for(Index i=0;i<nP;i++) {
        Pval[i] = x[(2*Time+1)*(nY+nU)+i];
     } //end for loop

    grad_f[jt+0*(Time+1)] += (-2*Xdval[0] + 2*Xval[0])/(2*Time+1);
    grad_f[jt+1*(Time+1)] += (-2*Xdval[1] + 2*Xval[1])/(2*Time+1);
    grad_f[jt+2*(Time+1)] += (-2*Xdval[2] + 2*Xval[2])/(2*Time+1);
    grad_f[jt+3*(Time+1)] += (-2*Xdval[3] + 2*Xval[3])/(2*Time+1);
    grad_f[(Time+1)*(nU+nY) + 0*Time + jt] += (-2*Xdval2[0] + 2*Xval2[0])/(2*Time+1);
    grad_f[(Time+1)*(nU+nY) + 1*Time + jt] += (-2*Xdval2[1] + 2*Xval2[1])/(2*Time+1);
    grad_f[(Time+1)*(nU+nY) + 2*Time + jt] += (-2*Xdval2[2] + 2*Xval2[2])/(2*Time+1);
    grad_f[(Time+1)*(nU+nY) + 3*Time + jt] += (-2*Xdval2[3] + 2*Xval2[3])/(2*Time+1);

    grad_f[jt+0*(Time+1)] += bounds[0][3]*((0.25*hstep*(-1.0L/4.0L*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xval[2]*pow(Xval[1], 3)) + 1.0)*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]) + (0.333333333333333*hstep*(-1.0L/4.0L*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xval[2]*pow(Xval[1], 3)) + 2)*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]))/(2*Time+1);
    grad_f[jt+1*(Time+1)] += bounds[0][3]*(0.25*Pval[0]*Xval[2]*hstep*pow(Xval[1], 2)*(Pval[1] - Xval[0])*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]) + 0.1875*Pval[0]*Xval[2]*hstep*pow(Xval[1], 2)*(Pval[1] - Xval[0])*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
    grad_f[jt+2*(Time+1)] += bounds[0][3]*(0.0833333333333333*Pval[0]*hstep*pow(Xval[1], 3)*(Pval[1] - Xval[0])*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]) + 0.0625*Pval[0]*hstep*pow(Xval[1], 3)*(Pval[1] - Xval[0])*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
    grad_f[jt+3*(Time+1)] += bounds[0][3]*(0.333333333333333*Pval[2]*hstep*pow(Xval[3], 3)*(Pval[3] - Xval[0])*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]) + 0.25*Pval[2]*hstep*pow(Xval[3], 3)*(Pval[3] - Xval[0])*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
    grad_f[jt+1+0*(Time+1)] += bounds[0][3]*((0.333333333333333*hstep*(-1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) - 2)*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]) + (0.25*hstep*((1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4) + (1.0L/4.0L)*Pval[4] + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) + 1.0)*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
    grad_f[jt+1+1*(Time+1)] += bounds[0][3]*(0.25*Pval[0]*Xvalp1[2]*hstep*pow(Xvalp1[1], 2)*(Pval[1] - Xvalp1[0])*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]) - 0.1875*Pval[0]*Xvalp1[2]*hstep*pow(Xvalp1[1], 2)*(Pval[1] - Xvalp1[0])*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
    grad_f[jt+1+2*(Time+1)] += bounds[0][3]*(0.0833333333333333*Pval[0]*hstep*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]) - 0.0625*Pval[0]*hstep*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
    grad_f[jt+1+3*(Time+1)] += bounds[0][3]*(0.333333333333333*Pval[2]*hstep*pow(Xvalp1[3], 3)*(Pval[3] - Xvalp1[0])*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]) - 0.25*Pval[2]*hstep*pow(Xvalp1[3], 3)*(Pval[3] - Xvalp1[0])*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
    grad_f[(Time+1)*(nU+nY) + 0*Time + jt] += bounds[0][3]*(0.333333333333333*hstep*(-Pval[2]*pow(Xval2[3], 4) - Pval[4] - Pval[0]*Xval2[2]*pow(Xval2[1], 3))*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]) - 0.25*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) - 1.0*Xval[0] + 2*Xval2[0] - 1.0*Xvalp1[0])/(2*Time+1);
    grad_f[(Time+1)*(nU+nY) + 1*Time + jt] += bounds[0][3]*(1.0*Pval[0]*Xval2[2]*hstep*pow(Xval2[1], 2)*(Pval[1] - Xval2[0])*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]))/(2*Time+1);
    grad_f[(Time+1)*(nU+nY) + 2*Time + jt] += bounds[0][3]*(0.333333333333333*Pval[0]*hstep*pow(Xval2[1], 3)*(Pval[1] - Xval2[0])*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]))/(2*Time+1);
    grad_f[(Time+1)*(nU+nY) + 3*Time + jt] += bounds[0][3]*(1.33333333333333*Pval[2]*hstep*pow(Xval2[3], 3)*(Pval[3] - Xval2[0])*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]))/(2*Time+1);
     grad_f[(2*Time+1)*(nY+nU)+0] += bounds[0][3]*(0.25*hstep*((1.0L/4.0L)*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0]))*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]) + 0.333333333333333*hstep*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0])*((1.0L/4.0L)*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])))/(2*Time+1);
     grad_f[(2*Time+1)*(nY+nU)+1] += bounds[0][3]*(0.25*hstep*((1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3))*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]) + 0.333333333333333*hstep*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0])*((1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3) + Pval[0]*Xval2[2]*pow(Xval2[1], 3) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)))/(2*Time+1);
     grad_f[(2*Time+1)*(nY+nU)+2] += bounds[0][3]*(0.25*hstep*((1.0L/4.0L)*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]))*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]) + 0.333333333333333*hstep*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0])*((1.0L/4.0L)*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0])))/(2*Time+1);
     grad_f[(2*Time+1)*(nY+nU)+3] += bounds[0][3]*(0.25*hstep*((1.0L/4.0L)*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4))*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]) + 0.333333333333333*hstep*((1.0L/4.0L)*Pval[2]*pow(Xval[3], 4) + Pval[2]*pow(Xval2[3], 4) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4))*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]))/(2*Time+1);
     grad_f[(2*Time+1)*(nY+nU)+4] += bounds[0][3]*(0.25*hstep*(-1.0L/4.0L*Xval[0] + (1.0L/4.0L)*Xvalp1[0])*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]) + 0.333333333333333*hstep*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0])*((3.0L/2.0L)*Pval[5] - 1.0L/4.0L*Xval[0] - Xval2[0] - 1.0L/4.0L*Xvalp1[0]))/(2*Time+1);
     grad_f[(2*Time+1)*(nY+nU)+5] += bounds[0][3]*(0.5*Pval[4]*hstep*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]))/(2*Time+1);
    grad_f[jt+0*(Time+1)] += bounds[1][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
// Betaspline
// Derivative
// Derivative
0.333333333333333*hstep*(-Xval[1]*Betasplinejac(Xval[0],1,1) + (-Xval[1] + 1)*Alphasplinejac(Xval[0],1,1))*(0.166666666666667*hstep*(-Xval[1]*Betaspline(Xval[0], 1) - 4*Xval2[1]*Betaspline(Xval2[0], 1) - Xvalp1[1]*Betaspline(Xvalp1[0], 1) + (-Xval[1] + 1)*Alphaspline(Xval[0], 1) + (-4*Xval2[1] + 4)*Alphaspline(Xval2[0], 1) + (-Xvalp1[1] + 1)*Alphaspline(Xvalp1[0], 1)) + Xval[1] - Xvalp1[1]) + 0.25*hstep*(-Xval[1]*Betasplinejac(Xval[0],1,1) + (-Xval[1] + 1)*Alphasplinejac(Xval[0],1,1))*(0.125*hstep*(-Xval[1]*Betaspline(Xval[0], 1) + Xvalp1[1]*Betaspline(Xvalp1[0], 1) + (-Xval[1] + 1)*Alphaspline(Xval[0], 1) - (-Xvalp1[1] + 1)*Alphaspline(Xvalp1[0], 1)) + 0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1]))/(2*Time+1);
    grad_f[jt+1*(Time+1)] += bounds[1][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
// Betaspline
(0.25*hstep*(-Alphaspline(Xval[0], 1) - Betaspline(Xval[0], 1)) + 1.0)*(0.125*hstep*(-Xval[1]*Betaspline(Xval[0], 1) + Xvalp1[1]*Betaspline(Xvalp1[0], 1) + (-Xval[1] + 1)*Alphaspline(Xval[0], 1) - (-Xvalp1[1] + 1)*Alphaspline(Xvalp1[0], 1)) + 0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1]) + (0.333333333333333*hstep*(-Alphaspline(Xval[0], 1) - Betaspline(Xval[0], 1)) + 2)*(0.166666666666667*hstep*(-Xval[1]*Betaspline(Xval[0], 1) - 4*Xval2[1]*Betaspline(Xval2[0], 1) - Xvalp1[1]*Betaspline(Xvalp1[0], 1) + (-Xval[1] + 1)*Alphaspline(Xval[0], 1) + (-4*Xval2[1] + 4)*Alphaspline(Xval2[0], 1) + (-Xvalp1[1] + 1)*Alphaspline(Xvalp1[0], 1)) + Xval[1] - Xvalp1[1]))/(2*Time+1);
    grad_f[jt+1+0*(Time+1)] += bounds[1][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
// Betaspline
// Derivative
// Derivative
0.333333333333333*hstep*(-Xvalp1[1]*Betasplinejac(Xvalp1[0],1,1) + (-Xvalp1[1] + 1)*Alphasplinejac(Xvalp1[0],1,1))*(0.166666666666667*hstep*(-Xval[1]*Betaspline(Xval[0], 1) - 4*Xval2[1]*Betaspline(Xval2[0], 1) - Xvalp1[1]*Betaspline(Xvalp1[0], 1) + (-Xval[1] + 1)*Alphaspline(Xval[0], 1) + (-4*Xval2[1] + 4)*Alphaspline(Xval2[0], 1) + (-Xvalp1[1] + 1)*Alphaspline(Xvalp1[0], 1)) + Xval[1] - Xvalp1[1]) + 0.25*hstep*(Xvalp1[1]*Betasplinejac(Xvalp1[0],1,1) - (-Xvalp1[1] + 1)*Alphasplinejac(Xvalp1[0],1,1))*(0.125*hstep*(-Xval[1]*Betaspline(Xval[0], 1) + Xvalp1[1]*Betaspline(Xvalp1[0], 1) + (-Xval[1] + 1)*Alphaspline(Xval[0], 1) - (-Xvalp1[1] + 1)*Alphaspline(Xvalp1[0], 1)) + 0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1]))/(2*Time+1);
    grad_f[jt+1+1*(Time+1)] += bounds[1][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
// Betaspline
(0.333333333333333*hstep*(-Alphaspline(Xvalp1[0], 1) - Betaspline(Xvalp1[0], 1)) - 2)*(0.166666666666667*hstep*(-Xval[1]*Betaspline(Xval[0], 1) - 4*Xval2[1]*Betaspline(Xval2[0], 1) - Xvalp1[1]*Betaspline(Xvalp1[0], 1) + (-Xval[1] + 1)*Alphaspline(Xval[0], 1) + (-4*Xval2[1] + 4)*Alphaspline(Xval2[0], 1) + (-Xvalp1[1] + 1)*Alphaspline(Xvalp1[0], 1)) + Xval[1] - Xvalp1[1]) + (0.25*hstep*(Alphaspline(Xvalp1[0], 1) + Betaspline(Xvalp1[0], 1)) + 1.0)*(0.125*hstep*(-Xval[1]*Betaspline(Xval[0], 1) + Xvalp1[1]*Betaspline(Xvalp1[0], 1) + (-Xval[1] + 1)*Alphaspline(Xval[0], 1) - (-Xvalp1[1] + 1)*Alphaspline(Xvalp1[0], 1)) + 0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1]))/(2*Time+1);
    grad_f[(Time+1)*(nU+nY) + 0*Time + jt] += bounds[1][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
// Betaspline
// Derivative
// Derivative
0.333333333333333*hstep*(-4*Xval2[1]*Betasplinejac(Xval2[0],1,1) + (-4*Xval2[1] + 4)*Alphasplinejac(Xval2[0],1,1))*(0.166666666666667*hstep*(-Xval[1]*Betaspline(Xval[0], 1) - 4*Xval2[1]*Betaspline(Xval2[0], 1) - Xvalp1[1]*Betaspline(Xvalp1[0], 1) + (-Xval[1] + 1)*Alphaspline(Xval[0], 1) + (-4*Xval2[1] + 4)*Alphaspline(Xval2[0], 1) + (-Xvalp1[1] + 1)*Alphaspline(Xvalp1[0], 1)) + Xval[1] - Xvalp1[1]))/(2*Time+1);
    grad_f[(Time+1)*(nU+nY) + 1*Time + jt] += bounds[1][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
// Betaspline
0.333333333333333*hstep*(-4*Alphaspline(Xval2[0], 1) - 4*Betaspline(Xval2[0], 1))*(0.166666666666667*hstep*(-Xval[1]*Betaspline(Xval[0], 1) - 4*Xval2[1]*Betaspline(Xval2[0], 1) - Xvalp1[1]*Betaspline(Xvalp1[0], 1) + (-Xval[1] + 1)*Alphaspline(Xval[0], 1) + (-4*Xval2[1] + 4)*Alphaspline(Xval2[0], 1) + (-Xvalp1[1] + 1)*Alphaspline(Xvalp1[0], 1)) + Xval[1] - Xvalp1[1]) - 0.25*hstep*(-Xval[1]*Betaspline(Xval[0], 1) + Xvalp1[1]*Betaspline(Xvalp1[0], 1) + (-Xval[1] + 1)*Alphaspline(Xval[0], 1) - (-Xvalp1[1] + 1)*Alphaspline(Xvalp1[0], 1)) - 1.0*Xval[1] + 2*Xval2[1] - 1.0*Xvalp1[1])/(2*Time+1);
    grad_f[jt+0*(Time+1)] += bounds[2][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
// Betaspline
// Derivative
// Derivative
0.333333333333333*hstep*(-Xval[2]*Betasplinejac(Xval[0],2,1) + (-Xval[2] + 1)*Alphasplinejac(Xval[0],2,1))*(Xval[2] - Xvalp1[2] + 0.166666666666667*hstep*(-Xval[2]*Betaspline(Xval[0], 2) - 4*Xval2[2]*Betaspline(Xval2[0], 2) - Xvalp1[2]*Betaspline(Xvalp1[0], 2) + (-Xval[2] + 1)*Alphaspline(Xval[0], 2) + (-4*Xval2[2] + 4)*Alphaspline(Xval2[0], 2) + (-Xvalp1[2] + 1)*Alphaspline(Xvalp1[0], 2))) + 0.25*hstep*(-Xval[2]*Betasplinejac(Xval[0],2,1) + (-Xval[2] + 1)*Alphasplinejac(Xval[0],2,1))*(0.5*Xval[2] - Xval2[2] + 0.5*Xvalp1[2] + 0.125*hstep*(-Xval[2]*Betaspline(Xval[0], 2) + Xvalp1[2]*Betaspline(Xvalp1[0], 2) + (-Xval[2] + 1)*Alphaspline(Xval[0], 2) - (-Xvalp1[2] + 1)*Alphaspline(Xvalp1[0], 2))))/(2*Time+1);
    grad_f[jt+2*(Time+1)] += bounds[2][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
// Betaspline
(0.25*hstep*(-Alphaspline(Xval[0], 2) - Betaspline(Xval[0], 2)) + 1.0)*(0.5*Xval[2] - Xval2[2] + 0.5*Xvalp1[2] + 0.125*hstep*(-Xval[2]*Betaspline(Xval[0], 2) + Xvalp1[2]*Betaspline(Xvalp1[0], 2) + (-Xval[2] + 1)*Alphaspline(Xval[0], 2) - (-Xvalp1[2] + 1)*Alphaspline(Xvalp1[0], 2))) + (0.333333333333333*hstep*(-Alphaspline(Xval[0], 2) - Betaspline(Xval[0], 2)) + 2)*(Xval[2] - Xvalp1[2] + 0.166666666666667*hstep*(-Xval[2]*Betaspline(Xval[0], 2) - 4*Xval2[2]*Betaspline(Xval2[0], 2) - Xvalp1[2]*Betaspline(Xvalp1[0], 2) + (-Xval[2] + 1)*Alphaspline(Xval[0], 2) + (-4*Xval2[2] + 4)*Alphaspline(Xval2[0], 2) + (-Xvalp1[2] + 1)*Alphaspline(Xvalp1[0], 2))))/(2*Time+1);
    grad_f[jt+1+0*(Time+1)] += bounds[2][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
// Betaspline
// Derivative
// Derivative
0.333333333333333*hstep*(-Xvalp1[2]*Betasplinejac(Xvalp1[0],2,1) + (-Xvalp1[2] + 1)*Alphasplinejac(Xvalp1[0],2,1))*(Xval[2] - Xvalp1[2] + 0.166666666666667*hstep*(-Xval[2]*Betaspline(Xval[0], 2) - 4*Xval2[2]*Betaspline(Xval2[0], 2) - Xvalp1[2]*Betaspline(Xvalp1[0], 2) + (-Xval[2] + 1)*Alphaspline(Xval[0], 2) + (-4*Xval2[2] + 4)*Alphaspline(Xval2[0], 2) + (-Xvalp1[2] + 1)*Alphaspline(Xvalp1[0], 2))) + 0.25*hstep*(Xvalp1[2]*Betasplinejac(Xvalp1[0],2,1) - (-Xvalp1[2] + 1)*Alphasplinejac(Xvalp1[0],2,1))*(0.5*Xval[2] - Xval2[2] + 0.5*Xvalp1[2] + 0.125*hstep*(-Xval[2]*Betaspline(Xval[0], 2) + Xvalp1[2]*Betaspline(Xvalp1[0], 2) + (-Xval[2] + 1)*Alphaspline(Xval[0], 2) - (-Xvalp1[2] + 1)*Alphaspline(Xvalp1[0], 2))))/(2*Time+1);
    grad_f[jt+1+2*(Time+1)] += bounds[2][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
// Betaspline
(0.333333333333333*hstep*(-Alphaspline(Xvalp1[0], 2) - Betaspline(Xvalp1[0], 2)) - 2)*(Xval[2] - Xvalp1[2] + 0.166666666666667*hstep*(-Xval[2]*Betaspline(Xval[0], 2) - 4*Xval2[2]*Betaspline(Xval2[0], 2) - Xvalp1[2]*Betaspline(Xvalp1[0], 2) + (-Xval[2] + 1)*Alphaspline(Xval[0], 2) + (-4*Xval2[2] + 4)*Alphaspline(Xval2[0], 2) + (-Xvalp1[2] + 1)*Alphaspline(Xvalp1[0], 2))) + (0.25*hstep*(Alphaspline(Xvalp1[0], 2) + Betaspline(Xvalp1[0], 2)) + 1.0)*(0.5*Xval[2] - Xval2[2] + 0.5*Xvalp1[2] + 0.125*hstep*(-Xval[2]*Betaspline(Xval[0], 2) + Xvalp1[2]*Betaspline(Xvalp1[0], 2) + (-Xval[2] + 1)*Alphaspline(Xval[0], 2) - (-Xvalp1[2] + 1)*Alphaspline(Xvalp1[0], 2))))/(2*Time+1);
    grad_f[(Time+1)*(nU+nY) + 0*Time + jt] += bounds[2][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
// Betaspline
// Derivative
// Derivative
0.333333333333333*hstep*(-4*Xval2[2]*Betasplinejac(Xval2[0],2,1) + (-4*Xval2[2] + 4)*Alphasplinejac(Xval2[0],2,1))*(Xval[2] - Xvalp1[2] + 0.166666666666667*hstep*(-Xval[2]*Betaspline(Xval[0], 2) - 4*Xval2[2]*Betaspline(Xval2[0], 2) - Xvalp1[2]*Betaspline(Xvalp1[0], 2) + (-Xval[2] + 1)*Alphaspline(Xval[0], 2) + (-4*Xval2[2] + 4)*Alphaspline(Xval2[0], 2) + (-Xvalp1[2] + 1)*Alphaspline(Xvalp1[0], 2))))/(2*Time+1);
    grad_f[(Time+1)*(nU+nY) + 2*Time + jt] += bounds[2][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
// Betaspline
-1.0*Xval[2] + 2*Xval2[2] - 1.0*Xvalp1[2] + 0.333333333333333*hstep*(-4*Alphaspline(Xval2[0], 2) - 4*Betaspline(Xval2[0], 2))*(Xval[2] - Xvalp1[2] + 0.166666666666667*hstep*(-Xval[2]*Betaspline(Xval[0], 2) - 4*Xval2[2]*Betaspline(Xval2[0], 2) - Xvalp1[2]*Betaspline(Xvalp1[0], 2) + (-Xval[2] + 1)*Alphaspline(Xval[0], 2) + (-4*Xval2[2] + 4)*Alphaspline(Xval2[0], 2) + (-Xvalp1[2] + 1)*Alphaspline(Xvalp1[0], 2))) - 0.25*hstep*(-Xval[2]*Betaspline(Xval[0], 2) + Xvalp1[2]*Betaspline(Xvalp1[0], 2) + (-Xval[2] + 1)*Alphaspline(Xval[0], 2) - (-Xvalp1[2] + 1)*Alphaspline(Xvalp1[0], 2)))/(2*Time+1);
    grad_f[jt+0*(Time+1)] += bounds[3][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
// Betaspline
// Derivative
// Derivative
0.333333333333333*hstep*(-Xval[3]*Betasplinejac(Xval[0],3,1) + (-Xval[3] + 1)*Alphasplinejac(Xval[0],3,1))*(0.166666666666667*hstep*(-Xval[3]*Betaspline(Xval[0], 3) - 4*Xval2[3]*Betaspline(Xval2[0], 3) - Xvalp1[3]*Betaspline(Xvalp1[0], 3) + (-Xval[3] + 1)*Alphaspline(Xval[0], 3) + (-4*Xval2[3] + 4)*Alphaspline(Xval2[0], 3) + (-Xvalp1[3] + 1)*Alphaspline(Xvalp1[0], 3)) + Xval[3] - Xvalp1[3]) + 0.25*hstep*(-Xval[3]*Betasplinejac(Xval[0],3,1) + (-Xval[3] + 1)*Alphasplinejac(Xval[0],3,1))*(0.125*hstep*(-Xval[3]*Betaspline(Xval[0], 3) + Xvalp1[3]*Betaspline(Xvalp1[0], 3) + (-Xval[3] + 1)*Alphaspline(Xval[0], 3) - (-Xvalp1[3] + 1)*Alphaspline(Xvalp1[0], 3)) + 0.5*Xval[3] - Xval2[3] + 0.5*Xvalp1[3]))/(2*Time+1);
    grad_f[jt+3*(Time+1)] += bounds[3][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
// Betaspline
(0.25*hstep*(-Alphaspline(Xval[0], 3) - Betaspline(Xval[0], 3)) + 1.0)*(0.125*hstep*(-Xval[3]*Betaspline(Xval[0], 3) + Xvalp1[3]*Betaspline(Xvalp1[0], 3) + (-Xval[3] + 1)*Alphaspline(Xval[0], 3) - (-Xvalp1[3] + 1)*Alphaspline(Xvalp1[0], 3)) + 0.5*Xval[3] - Xval2[3] + 0.5*Xvalp1[3]) + (0.333333333333333*hstep*(-Alphaspline(Xval[0], 3) - Betaspline(Xval[0], 3)) + 2)*(0.166666666666667*hstep*(-Xval[3]*Betaspline(Xval[0], 3) - 4*Xval2[3]*Betaspline(Xval2[0], 3) - Xvalp1[3]*Betaspline(Xvalp1[0], 3) + (-Xval[3] + 1)*Alphaspline(Xval[0], 3) + (-4*Xval2[3] + 4)*Alphaspline(Xval2[0], 3) + (-Xvalp1[3] + 1)*Alphaspline(Xvalp1[0], 3)) + Xval[3] - Xvalp1[3]))/(2*Time+1);
    grad_f[jt+1+0*(Time+1)] += bounds[3][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
// Betaspline
// Derivative
// Derivative
0.333333333333333*hstep*(-Xvalp1[3]*Betasplinejac(Xvalp1[0],3,1) + (-Xvalp1[3] + 1)*Alphasplinejac(Xvalp1[0],3,1))*(0.166666666666667*hstep*(-Xval[3]*Betaspline(Xval[0], 3) - 4*Xval2[3]*Betaspline(Xval2[0], 3) - Xvalp1[3]*Betaspline(Xvalp1[0], 3) + (-Xval[3] + 1)*Alphaspline(Xval[0], 3) + (-4*Xval2[3] + 4)*Alphaspline(Xval2[0], 3) + (-Xvalp1[3] + 1)*Alphaspline(Xvalp1[0], 3)) + Xval[3] - Xvalp1[3]) + 0.25*hstep*(Xvalp1[3]*Betasplinejac(Xvalp1[0],3,1) - (-Xvalp1[3] + 1)*Alphasplinejac(Xvalp1[0],3,1))*(0.125*hstep*(-Xval[3]*Betaspline(Xval[0], 3) + Xvalp1[3]*Betaspline(Xvalp1[0], 3) + (-Xval[3] + 1)*Alphaspline(Xval[0], 3) - (-Xvalp1[3] + 1)*Alphaspline(Xvalp1[0], 3)) + 0.5*Xval[3] - Xval2[3] + 0.5*Xvalp1[3]))/(2*Time+1);
    grad_f[jt+1+3*(Time+1)] += bounds[3][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
// Betaspline
(0.333333333333333*hstep*(-Alphaspline(Xvalp1[0], 3) - Betaspline(Xvalp1[0], 3)) - 2)*(0.166666666666667*hstep*(-Xval[3]*Betaspline(Xval[0], 3) - 4*Xval2[3]*Betaspline(Xval2[0], 3) - Xvalp1[3]*Betaspline(Xvalp1[0], 3) + (-Xval[3] + 1)*Alphaspline(Xval[0], 3) + (-4*Xval2[3] + 4)*Alphaspline(Xval2[0], 3) + (-Xvalp1[3] + 1)*Alphaspline(Xvalp1[0], 3)) + Xval[3] - Xvalp1[3]) + (0.25*hstep*(Alphaspline(Xvalp1[0], 3) + Betaspline(Xvalp1[0], 3)) + 1.0)*(0.125*hstep*(-Xval[3]*Betaspline(Xval[0], 3) + Xvalp1[3]*Betaspline(Xvalp1[0], 3) + (-Xval[3] + 1)*Alphaspline(Xval[0], 3) - (-Xvalp1[3] + 1)*Alphaspline(Xvalp1[0], 3)) + 0.5*Xval[3] - Xval2[3] + 0.5*Xvalp1[3]))/(2*Time+1);
    grad_f[(Time+1)*(nU+nY) + 0*Time + jt] += bounds[3][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
// Betaspline
// Derivative
// Derivative
0.333333333333333*hstep*(-4*Xval2[3]*Betasplinejac(Xval2[0],3,1) + (-4*Xval2[3] + 4)*Alphasplinejac(Xval2[0],3,1))*(0.166666666666667*hstep*(-Xval[3]*Betaspline(Xval[0], 3) - 4*Xval2[3]*Betaspline(Xval2[0], 3) - Xvalp1[3]*Betaspline(Xvalp1[0], 3) + (-Xval[3] + 1)*Alphaspline(Xval[0], 3) + (-4*Xval2[3] + 4)*Alphaspline(Xval2[0], 3) + (-Xvalp1[3] + 1)*Alphaspline(Xvalp1[0], 3)) + Xval[3] - Xvalp1[3]))/(2*Time+1);
    grad_f[(Time+1)*(nU+nY) + 3*Time + jt] += bounds[3][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
// Betaspline
0.333333333333333*hstep*(-4*Alphaspline(Xval2[0], 3) - 4*Betaspline(Xval2[0], 3))*(0.166666666666667*hstep*(-Xval[3]*Betaspline(Xval[0], 3) - 4*Xval2[3]*Betaspline(Xval2[0], 3) - Xvalp1[3]*Betaspline(Xvalp1[0], 3) + (-Xval[3] + 1)*Alphaspline(Xval[0], 3) + (-4*Xval2[3] + 4)*Alphaspline(Xval2[0], 3) + (-Xvalp1[3] + 1)*Alphaspline(Xvalp1[0], 3)) + Xval[3] - Xvalp1[3]) - 0.25*hstep*(-Xval[3]*Betaspline(Xval[0], 3) + Xvalp1[3]*Betaspline(Xvalp1[0], 3) + (-Xval[3] + 1)*Alphaspline(Xval[0], 3) - (-Xvalp1[3] + 1)*Alphaspline(Xvalp1[0], 3)) - 1.0*Xval[3] + 2*Xval2[3] - 1.0*Xvalp1[3])/(2*Time+1);

  } //end for loop

// Add last element
     for(Index i=0;i<nY;i++) {
        Xval[i] = x[i*(Time+1) + Time];
        Xvalp1[i] = 0;
        Xval2[i] = 0;
     } //end for loop

     for(Index cup=0;cup<nU;cup++) {
        K11val[cup] = x[nY*(Time+1) + cup*(Time+1) + Time];
        K11valp1[cup] = 0;
        K11val2[cup] = 0;
     } //end for loop

     Xdval[0] = VDATA0[2*Time];
     Xdval2[0] = 0;
     Xdvalp1[0] = 0;
     Xdval[1] = MDATA0[2*Time];
     Xdval2[1] = 0;
     Xdvalp1[1] = 0;
     Xdval[2] = HDATA0[2*Time];
     Xdval2[2] = 0;
     Xdvalp1[2] = 0;
     Xdval[3] = NDATA0[2*Time];
     Xdval2[3] = 0;
     Xdvalp1[3] = 0;
     Ival[0] = Iinj[2*Time];
     Ival2[0] = 0;
     Ivalp1[0] = 0;

     for(Index i=0;i<nP;i++) {
        Pval[i] = x[(2*Time+1)*(nY+nU)+i];
     } //end for loop


    grad_f[Time+0*(Time+1)] += (-2*Xdval[0] + 2*Xval[0])/(2*Time+1);
    grad_f[Time+1*(Time+1)] += (-2*Xdval[1] + 2*Xval[1])/(2*Time+1);
    grad_f[Time+2*(Time+1)] += (-2*Xdval[2] + 2*Xval[2])/(2*Time+1);
    grad_f[Time+3*(Time+1)] += (-2*Xdval[3] + 2*Xval[3])/(2*Time+1);

  return true;
}


// return the value of the constraints: g(x)
bool NEURODYN_NAKL_NLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  assert(n == 8*Time+10);
  assert(m == 0);

  return true;
}


// return the structure or values of the jacobian
bool NEURODYN_NAKL_NLP::eval_jac_g(Index n, const Number* x, bool new_x,
                            Index m, Index nele_jac, Index* iRow, Index* jCol,
                            Number* values)
{

return true;
}


// return the structure or values of the hessian
bool NEURODYN_NAKL_NLP::eval_h(Index n, const Number* x, bool new_x,
                       Number obj_factor, Index m, const Number* lambda,
                       bool new_lambda, Index nele_hess, Index* iRow,
                       Index* jCol, Number* values)
{

if (values == NULL) {
   // return the structure.  This is a symmetric matrix, fill in the lower left
   // triangle only.

   // Each non-one Hessian element has its own explicit loop
   // since each element needs a different number of matrix elements

     for(Index jt=0;jt<Time+1;jt++) {
     iRow[0*(Time+1)+0*(Time)+0+jt] = (Time+1)*0+jt;
     jCol[0*(Time+1)+0*(Time)+0+jt] = (Time+1)*0+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[1*(Time+1)+0*(Time)+0+jt] = (Time+1)*0+jt+1;
     jCol[1*(Time+1)+0*(Time)+0+jt] = (Time+1)*0+jt;
   }
     for(Index jt=0;jt<Time+1;jt++) {
     iRow[1*(Time+1)+1*(Time)+0+jt] = (Time+1)*1+jt;
     jCol[1*(Time+1)+1*(Time)+0+jt] = (Time+1)*0+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[2*(Time+1)+1*(Time)+0+jt] = (Time+1)*1+jt;
     jCol[2*(Time+1)+1*(Time)+0+jt] = (Time+1)*0+jt+1;
   }
     for(Index jt=0;jt<Time+1;jt++) {
     iRow[2*(Time+1)+2*(Time)+0+jt] = (Time+1)*1+jt;
     jCol[2*(Time+1)+2*(Time)+0+jt] = (Time+1)*1+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[3*(Time+1)+2*(Time)+0+jt] = (Time+1)*1+jt+1;
     jCol[3*(Time+1)+2*(Time)+0+jt] = (Time+1)*0+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[3*(Time+1)+3*(Time)+0+jt] = (Time+1)*1+jt+1;
     jCol[3*(Time+1)+3*(Time)+0+jt] = (Time+1)*1+jt;
   }
     for(Index jt=0;jt<Time+1;jt++) {
     iRow[3*(Time+1)+4*(Time)+0+jt] = (Time+1)*2+jt;
     jCol[3*(Time+1)+4*(Time)+0+jt] = (Time+1)*0+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[4*(Time+1)+4*(Time)+0+jt] = (Time+1)*2+jt;
     jCol[4*(Time+1)+4*(Time)+0+jt] = (Time+1)*0+jt+1;
   }
     for(Index jt=0;jt<Time+1;jt++) {
     iRow[4*(Time+1)+5*(Time)+0+jt] = (Time+1)*2+jt;
     jCol[4*(Time+1)+5*(Time)+0+jt] = (Time+1)*1+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[5*(Time+1)+5*(Time)+0+jt] = (Time+1)*2+jt;
     jCol[5*(Time+1)+5*(Time)+0+jt] = (Time+1)*1+jt+1;
   }
     for(Index jt=0;jt<Time+1;jt++) {
     iRow[5*(Time+1)+6*(Time)+0+jt] = (Time+1)*2+jt;
     jCol[5*(Time+1)+6*(Time)+0+jt] = (Time+1)*2+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[6*(Time+1)+6*(Time)+0+jt] = (Time+1)*2+jt+1;
     jCol[6*(Time+1)+6*(Time)+0+jt] = (Time+1)*0+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[6*(Time+1)+7*(Time)+0+jt] = (Time+1)*2+jt+1;
     jCol[6*(Time+1)+7*(Time)+0+jt] = (Time+1)*1+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[6*(Time+1)+8*(Time)+0+jt] = (Time+1)*2+jt+1;
     jCol[6*(Time+1)+8*(Time)+0+jt] = (Time+1)*2+jt;
   }
     for(Index jt=0;jt<Time+1;jt++) {
     iRow[6*(Time+1)+9*(Time)+0+jt] = (Time+1)*3+jt;
     jCol[6*(Time+1)+9*(Time)+0+jt] = (Time+1)*0+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[7*(Time+1)+9*(Time)+0+jt] = (Time+1)*3+jt;
     jCol[7*(Time+1)+9*(Time)+0+jt] = (Time+1)*0+jt+1;
   }
     for(Index jt=0;jt<Time+1;jt++) {
     iRow[7*(Time+1)+10*(Time)+0+jt] = (Time+1)*3+jt;
     jCol[7*(Time+1)+10*(Time)+0+jt] = (Time+1)*1+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[8*(Time+1)+10*(Time)+0+jt] = (Time+1)*3+jt;
     jCol[8*(Time+1)+10*(Time)+0+jt] = (Time+1)*1+jt+1;
   }
     for(Index jt=0;jt<Time+1;jt++) {
     iRow[8*(Time+1)+11*(Time)+0+jt] = (Time+1)*3+jt;
     jCol[8*(Time+1)+11*(Time)+0+jt] = (Time+1)*2+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[9*(Time+1)+11*(Time)+0+jt] = (Time+1)*3+jt;
     jCol[9*(Time+1)+11*(Time)+0+jt] = (Time+1)*2+jt+1;
   }
     for(Index jt=0;jt<Time+1;jt++) {
     iRow[9*(Time+1)+12*(Time)+0+jt] = (Time+1)*3+jt;
     jCol[9*(Time+1)+12*(Time)+0+jt] = (Time+1)*3+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+12*(Time)+0+jt] = (Time+1)*3+jt+1;
     jCol[10*(Time+1)+12*(Time)+0+jt] = (Time+1)*0+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+13*(Time)+0+jt] = (Time+1)*3+jt+1;
     jCol[10*(Time+1)+13*(Time)+0+jt] = (Time+1)*1+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+14*(Time)+0+jt] = (Time+1)*3+jt+1;
     jCol[10*(Time+1)+14*(Time)+0+jt] = (Time+1)*2+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+15*(Time)+0+jt] = (Time+1)*3+jt+1;
     jCol[10*(Time+1)+15*(Time)+0+jt] = (Time+1)*3+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+16*(Time)+0+jt] = (Time+1)*4+Time*0+jt;
     jCol[10*(Time+1)+16*(Time)+0+jt] = (Time+1)*0+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+17*(Time)+0+jt] = (Time+1)*4+Time*0+jt;
     jCol[10*(Time+1)+17*(Time)+0+jt] = (Time+1)*0+jt+1;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+18*(Time)+0+jt] = (Time+1)*4+Time*0+jt;
     jCol[10*(Time+1)+18*(Time)+0+jt] = (Time+1)*1+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+19*(Time)+0+jt] = (Time+1)*4+Time*0+jt;
     jCol[10*(Time+1)+19*(Time)+0+jt] = (Time+1)*1+jt+1;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+20*(Time)+0+jt] = (Time+1)*4+Time*0+jt;
     jCol[10*(Time+1)+20*(Time)+0+jt] = (Time+1)*2+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+21*(Time)+0+jt] = (Time+1)*4+Time*0+jt;
     jCol[10*(Time+1)+21*(Time)+0+jt] = (Time+1)*2+jt+1;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+22*(Time)+0+jt] = (Time+1)*4+Time*0+jt;
     jCol[10*(Time+1)+22*(Time)+0+jt] = (Time+1)*3+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+23*(Time)+0+jt] = (Time+1)*4+Time*0+jt;
     jCol[10*(Time+1)+23*(Time)+0+jt] = (Time+1)*3+jt+1;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+24*(Time)+0+jt] = (Time+1)*4+Time*0+jt;
     jCol[10*(Time+1)+24*(Time)+0+jt] = (Time+1)*4+Time*0+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+25*(Time)+0+jt] = (Time+1)*4+Time*1+jt;
     jCol[10*(Time+1)+25*(Time)+0+jt] = (Time+1)*0+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+26*(Time)+0+jt] = (Time+1)*4+Time*1+jt;
     jCol[10*(Time+1)+26*(Time)+0+jt] = (Time+1)*0+jt+1;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+27*(Time)+0+jt] = (Time+1)*4+Time*1+jt;
     jCol[10*(Time+1)+27*(Time)+0+jt] = (Time+1)*1+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+28*(Time)+0+jt] = (Time+1)*4+Time*1+jt;
     jCol[10*(Time+1)+28*(Time)+0+jt] = (Time+1)*1+jt+1;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+29*(Time)+0+jt] = (Time+1)*4+Time*1+jt;
     jCol[10*(Time+1)+29*(Time)+0+jt] = (Time+1)*2+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+30*(Time)+0+jt] = (Time+1)*4+Time*1+jt;
     jCol[10*(Time+1)+30*(Time)+0+jt] = (Time+1)*2+jt+1;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+31*(Time)+0+jt] = (Time+1)*4+Time*1+jt;
     jCol[10*(Time+1)+31*(Time)+0+jt] = (Time+1)*3+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+32*(Time)+0+jt] = (Time+1)*4+Time*1+jt;
     jCol[10*(Time+1)+32*(Time)+0+jt] = (Time+1)*3+jt+1;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+33*(Time)+0+jt] = (Time+1)*4+Time*1+jt;
     jCol[10*(Time+1)+33*(Time)+0+jt] = (Time+1)*4+Time*0+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+34*(Time)+0+jt] = (Time+1)*4+Time*1+jt;
     jCol[10*(Time+1)+34*(Time)+0+jt] = (Time+1)*4+Time*1+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+35*(Time)+0+jt] = (Time+1)*4+Time*2+jt;
     jCol[10*(Time+1)+35*(Time)+0+jt] = (Time+1)*0+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+36*(Time)+0+jt] = (Time+1)*4+Time*2+jt;
     jCol[10*(Time+1)+36*(Time)+0+jt] = (Time+1)*0+jt+1;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+37*(Time)+0+jt] = (Time+1)*4+Time*2+jt;
     jCol[10*(Time+1)+37*(Time)+0+jt] = (Time+1)*1+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+38*(Time)+0+jt] = (Time+1)*4+Time*2+jt;
     jCol[10*(Time+1)+38*(Time)+0+jt] = (Time+1)*1+jt+1;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+39*(Time)+0+jt] = (Time+1)*4+Time*2+jt;
     jCol[10*(Time+1)+39*(Time)+0+jt] = (Time+1)*2+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+40*(Time)+0+jt] = (Time+1)*4+Time*2+jt;
     jCol[10*(Time+1)+40*(Time)+0+jt] = (Time+1)*2+jt+1;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+41*(Time)+0+jt] = (Time+1)*4+Time*2+jt;
     jCol[10*(Time+1)+41*(Time)+0+jt] = (Time+1)*3+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+42*(Time)+0+jt] = (Time+1)*4+Time*2+jt;
     jCol[10*(Time+1)+42*(Time)+0+jt] = (Time+1)*3+jt+1;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+43*(Time)+0+jt] = (Time+1)*4+Time*2+jt;
     jCol[10*(Time+1)+43*(Time)+0+jt] = (Time+1)*4+Time*0+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+44*(Time)+0+jt] = (Time+1)*4+Time*2+jt;
     jCol[10*(Time+1)+44*(Time)+0+jt] = (Time+1)*4+Time*1+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+45*(Time)+0+jt] = (Time+1)*4+Time*2+jt;
     jCol[10*(Time+1)+45*(Time)+0+jt] = (Time+1)*4+Time*2+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+46*(Time)+0+jt] = (Time+1)*4+Time*3+jt;
     jCol[10*(Time+1)+46*(Time)+0+jt] = (Time+1)*0+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+47*(Time)+0+jt] = (Time+1)*4+Time*3+jt;
     jCol[10*(Time+1)+47*(Time)+0+jt] = (Time+1)*0+jt+1;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+48*(Time)+0+jt] = (Time+1)*4+Time*3+jt;
     jCol[10*(Time+1)+48*(Time)+0+jt] = (Time+1)*1+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+49*(Time)+0+jt] = (Time+1)*4+Time*3+jt;
     jCol[10*(Time+1)+49*(Time)+0+jt] = (Time+1)*1+jt+1;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+50*(Time)+0+jt] = (Time+1)*4+Time*3+jt;
     jCol[10*(Time+1)+50*(Time)+0+jt] = (Time+1)*2+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+51*(Time)+0+jt] = (Time+1)*4+Time*3+jt;
     jCol[10*(Time+1)+51*(Time)+0+jt] = (Time+1)*2+jt+1;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+52*(Time)+0+jt] = (Time+1)*4+Time*3+jt;
     jCol[10*(Time+1)+52*(Time)+0+jt] = (Time+1)*3+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+53*(Time)+0+jt] = (Time+1)*4+Time*3+jt;
     jCol[10*(Time+1)+53*(Time)+0+jt] = (Time+1)*3+jt+1;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+54*(Time)+0+jt] = (Time+1)*4+Time*3+jt;
     jCol[10*(Time+1)+54*(Time)+0+jt] = (Time+1)*4+Time*0+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+55*(Time)+0+jt] = (Time+1)*4+Time*3+jt;
     jCol[10*(Time+1)+55*(Time)+0+jt] = (Time+1)*4+Time*1+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+56*(Time)+0+jt] = (Time+1)*4+Time*3+jt;
     jCol[10*(Time+1)+56*(Time)+0+jt] = (Time+1)*4+Time*2+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[10*(Time+1)+57*(Time)+0+jt] = (Time+1)*4+Time*3+jt;
     jCol[10*(Time+1)+57*(Time)+0+jt] = (Time+1)*4+Time*3+jt;
   }
     for(Index jt=0;jt<Time+1;jt++) {
     iRow[10*(Time+1)+58*(Time)+0+jt] = 2*Time*4+4;
     jCol[10*(Time+1)+58*(Time)+0+jt] = (Time+1)*0+jt;
   }
     for(Index jt=0;jt<Time+1;jt++) {
     iRow[11*(Time+1)+58*(Time)+0+jt] = 2*Time*4+4;
     jCol[11*(Time+1)+58*(Time)+0+jt] = (Time+1)*1+jt;
   }
     for(Index jt=0;jt<Time+1;jt++) {
     iRow[12*(Time+1)+58*(Time)+0+jt] = 2*Time*4+4;
     jCol[12*(Time+1)+58*(Time)+0+jt] = (Time+1)*2+jt;
   }
     for(Index jt=0;jt<Time+1;jt++) {
     iRow[13*(Time+1)+58*(Time)+0+jt] = 2*Time*4+4;
     jCol[13*(Time+1)+58*(Time)+0+jt] = (Time+1)*3+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[14*(Time+1)+58*(Time)+0+jt] = 2*Time*4+4;
     jCol[14*(Time+1)+58*(Time)+0+jt] = (Time+1)*4+Time*0+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[14*(Time+1)+59*(Time)+0+jt] = 2*Time*4+4;
     jCol[14*(Time+1)+59*(Time)+0+jt] = (Time+1)*4+Time*1+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[14*(Time+1)+60*(Time)+0+jt] = 2*Time*4+4;
     jCol[14*(Time+1)+60*(Time)+0+jt] = (Time+1)*4+Time*2+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[14*(Time+1)+61*(Time)+0+jt] = 2*Time*4+4;
     jCol[14*(Time+1)+61*(Time)+0+jt] = (Time+1)*4+Time*3+jt;
   }
     for(Index jt=0;jt<1;jt++) {
     iRow[14*(Time+1)+62*(Time)+0+jt] = 2*Time*4+4;
     jCol[14*(Time+1)+62*(Time)+0+jt] = 2*Time*4+4;
   }
     for(Index jt=0;jt<Time+1;jt++) {
     iRow[14*(Time+1)+62*(Time)+1+jt] = 2*Time*4+5;
     jCol[14*(Time+1)+62*(Time)+1+jt] = (Time+1)*0+jt;
   }
     for(Index jt=0;jt<Time+1;jt++) {
     iRow[15*(Time+1)+62*(Time)+1+jt] = 2*Time*4+5;
     jCol[15*(Time+1)+62*(Time)+1+jt] = (Time+1)*1+jt;
   }
     for(Index jt=0;jt<Time+1;jt++) {
     iRow[16*(Time+1)+62*(Time)+1+jt] = 2*Time*4+5;
     jCol[16*(Time+1)+62*(Time)+1+jt] = (Time+1)*2+jt;
   }
     for(Index jt=0;jt<Time+1;jt++) {
     iRow[17*(Time+1)+62*(Time)+1+jt] = 2*Time*4+5;
     jCol[17*(Time+1)+62*(Time)+1+jt] = (Time+1)*3+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[18*(Time+1)+62*(Time)+1+jt] = 2*Time*4+5;
     jCol[18*(Time+1)+62*(Time)+1+jt] = (Time+1)*4+Time*0+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[18*(Time+1)+63*(Time)+1+jt] = 2*Time*4+5;
     jCol[18*(Time+1)+63*(Time)+1+jt] = (Time+1)*4+Time*1+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[18*(Time+1)+64*(Time)+1+jt] = 2*Time*4+5;
     jCol[18*(Time+1)+64*(Time)+1+jt] = (Time+1)*4+Time*2+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[18*(Time+1)+65*(Time)+1+jt] = 2*Time*4+5;
     jCol[18*(Time+1)+65*(Time)+1+jt] = (Time+1)*4+Time*3+jt;
   }
     for(Index jt=0;jt<1;jt++) {
     iRow[18*(Time+1)+66*(Time)+1+jt] = 2*Time*4+5;
     jCol[18*(Time+1)+66*(Time)+1+jt] = 2*Time*4+4;
   }
     for(Index jt=0;jt<1;jt++) {
     iRow[18*(Time+1)+66*(Time)+2+jt] = 2*Time*4+5;
     jCol[18*(Time+1)+66*(Time)+2+jt] = 2*Time*4+5;
   }
     for(Index jt=0;jt<Time+1;jt++) {
     iRow[18*(Time+1)+66*(Time)+3+jt] = 2*Time*4+6;
     jCol[18*(Time+1)+66*(Time)+3+jt] = (Time+1)*0+jt;
   }
     for(Index jt=0;jt<Time+1;jt++) {
     iRow[19*(Time+1)+66*(Time)+3+jt] = 2*Time*4+6;
     jCol[19*(Time+1)+66*(Time)+3+jt] = (Time+1)*1+jt;
   }
     for(Index jt=0;jt<Time+1;jt++) {
     iRow[20*(Time+1)+66*(Time)+3+jt] = 2*Time*4+6;
     jCol[20*(Time+1)+66*(Time)+3+jt] = (Time+1)*2+jt;
   }
     for(Index jt=0;jt<Time+1;jt++) {
     iRow[21*(Time+1)+66*(Time)+3+jt] = 2*Time*4+6;
     jCol[21*(Time+1)+66*(Time)+3+jt] = (Time+1)*3+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+66*(Time)+3+jt] = 2*Time*4+6;
     jCol[22*(Time+1)+66*(Time)+3+jt] = (Time+1)*4+Time*0+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+67*(Time)+3+jt] = 2*Time*4+6;
     jCol[22*(Time+1)+67*(Time)+3+jt] = (Time+1)*4+Time*1+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+68*(Time)+3+jt] = 2*Time*4+6;
     jCol[22*(Time+1)+68*(Time)+3+jt] = (Time+1)*4+Time*2+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[22*(Time+1)+69*(Time)+3+jt] = 2*Time*4+6;
     jCol[22*(Time+1)+69*(Time)+3+jt] = (Time+1)*4+Time*3+jt;
   }
     for(Index jt=0;jt<1;jt++) {
     iRow[22*(Time+1)+70*(Time)+3+jt] = 2*Time*4+6;
     jCol[22*(Time+1)+70*(Time)+3+jt] = 2*Time*4+4;
   }
     for(Index jt=0;jt<1;jt++) {
     iRow[22*(Time+1)+70*(Time)+4+jt] = 2*Time*4+6;
     jCol[22*(Time+1)+70*(Time)+4+jt] = 2*Time*4+5;
   }
     for(Index jt=0;jt<1;jt++) {
     iRow[22*(Time+1)+70*(Time)+5+jt] = 2*Time*4+6;
     jCol[22*(Time+1)+70*(Time)+5+jt] = 2*Time*4+6;
   }
     for(Index jt=0;jt<Time+1;jt++) {
     iRow[22*(Time+1)+70*(Time)+6+jt] = 2*Time*4+7;
     jCol[22*(Time+1)+70*(Time)+6+jt] = (Time+1)*0+jt;
   }
     for(Index jt=0;jt<Time+1;jt++) {
     iRow[23*(Time+1)+70*(Time)+6+jt] = 2*Time*4+7;
     jCol[23*(Time+1)+70*(Time)+6+jt] = (Time+1)*1+jt;
   }
     for(Index jt=0;jt<Time+1;jt++) {
     iRow[24*(Time+1)+70*(Time)+6+jt] = 2*Time*4+7;
     jCol[24*(Time+1)+70*(Time)+6+jt] = (Time+1)*2+jt;
   }
     for(Index jt=0;jt<Time+1;jt++) {
     iRow[25*(Time+1)+70*(Time)+6+jt] = 2*Time*4+7;
     jCol[25*(Time+1)+70*(Time)+6+jt] = (Time+1)*3+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[26*(Time+1)+70*(Time)+6+jt] = 2*Time*4+7;
     jCol[26*(Time+1)+70*(Time)+6+jt] = (Time+1)*4+Time*0+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[26*(Time+1)+71*(Time)+6+jt] = 2*Time*4+7;
     jCol[26*(Time+1)+71*(Time)+6+jt] = (Time+1)*4+Time*1+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[26*(Time+1)+72*(Time)+6+jt] = 2*Time*4+7;
     jCol[26*(Time+1)+72*(Time)+6+jt] = (Time+1)*4+Time*2+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[26*(Time+1)+73*(Time)+6+jt] = 2*Time*4+7;
     jCol[26*(Time+1)+73*(Time)+6+jt] = (Time+1)*4+Time*3+jt;
   }
     for(Index jt=0;jt<1;jt++) {
     iRow[26*(Time+1)+74*(Time)+6+jt] = 2*Time*4+7;
     jCol[26*(Time+1)+74*(Time)+6+jt] = 2*Time*4+4;
   }
     for(Index jt=0;jt<1;jt++) {
     iRow[26*(Time+1)+74*(Time)+7+jt] = 2*Time*4+7;
     jCol[26*(Time+1)+74*(Time)+7+jt] = 2*Time*4+5;
   }
     for(Index jt=0;jt<1;jt++) {
     iRow[26*(Time+1)+74*(Time)+8+jt] = 2*Time*4+7;
     jCol[26*(Time+1)+74*(Time)+8+jt] = 2*Time*4+6;
   }
     for(Index jt=0;jt<1;jt++) {
     iRow[26*(Time+1)+74*(Time)+9+jt] = 2*Time*4+7;
     jCol[26*(Time+1)+74*(Time)+9+jt] = 2*Time*4+7;
   }
     for(Index jt=0;jt<Time+1;jt++) {
     iRow[26*(Time+1)+74*(Time)+10+jt] = 2*Time*4+8;
     jCol[26*(Time+1)+74*(Time)+10+jt] = (Time+1)*0+jt;
   }
     for(Index jt=0;jt<Time+1;jt++) {
     iRow[27*(Time+1)+74*(Time)+10+jt] = 2*Time*4+8;
     jCol[27*(Time+1)+74*(Time)+10+jt] = (Time+1)*1+jt;
   }
     for(Index jt=0;jt<Time+1;jt++) {
     iRow[28*(Time+1)+74*(Time)+10+jt] = 2*Time*4+8;
     jCol[28*(Time+1)+74*(Time)+10+jt] = (Time+1)*2+jt;
   }
     for(Index jt=0;jt<Time+1;jt++) {
     iRow[29*(Time+1)+74*(Time)+10+jt] = 2*Time*4+8;
     jCol[29*(Time+1)+74*(Time)+10+jt] = (Time+1)*3+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[30*(Time+1)+74*(Time)+10+jt] = 2*Time*4+8;
     jCol[30*(Time+1)+74*(Time)+10+jt] = (Time+1)*4+Time*0+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[30*(Time+1)+75*(Time)+10+jt] = 2*Time*4+8;
     jCol[30*(Time+1)+75*(Time)+10+jt] = (Time+1)*4+Time*1+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[30*(Time+1)+76*(Time)+10+jt] = 2*Time*4+8;
     jCol[30*(Time+1)+76*(Time)+10+jt] = (Time+1)*4+Time*2+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[30*(Time+1)+77*(Time)+10+jt] = 2*Time*4+8;
     jCol[30*(Time+1)+77*(Time)+10+jt] = (Time+1)*4+Time*3+jt;
   }
     for(Index jt=0;jt<1;jt++) {
     iRow[30*(Time+1)+78*(Time)+10+jt] = 2*Time*4+8;
     jCol[30*(Time+1)+78*(Time)+10+jt] = 2*Time*4+4;
   }
     for(Index jt=0;jt<1;jt++) {
     iRow[30*(Time+1)+78*(Time)+11+jt] = 2*Time*4+8;
     jCol[30*(Time+1)+78*(Time)+11+jt] = 2*Time*4+5;
   }
     for(Index jt=0;jt<1;jt++) {
     iRow[30*(Time+1)+78*(Time)+12+jt] = 2*Time*4+8;
     jCol[30*(Time+1)+78*(Time)+12+jt] = 2*Time*4+6;
   }
     for(Index jt=0;jt<1;jt++) {
     iRow[30*(Time+1)+78*(Time)+13+jt] = 2*Time*4+8;
     jCol[30*(Time+1)+78*(Time)+13+jt] = 2*Time*4+7;
   }
     for(Index jt=0;jt<1;jt++) {
     iRow[30*(Time+1)+78*(Time)+14+jt] = 2*Time*4+8;
     jCol[30*(Time+1)+78*(Time)+14+jt] = 2*Time*4+8;
   }
     for(Index jt=0;jt<Time+1;jt++) {
     iRow[30*(Time+1)+78*(Time)+15+jt] = 2*Time*4+9;
     jCol[30*(Time+1)+78*(Time)+15+jt] = (Time+1)*0+jt;
   }
     for(Index jt=0;jt<Time+1;jt++) {
     iRow[31*(Time+1)+78*(Time)+15+jt] = 2*Time*4+9;
     jCol[31*(Time+1)+78*(Time)+15+jt] = (Time+1)*1+jt;
   }
     for(Index jt=0;jt<Time+1;jt++) {
     iRow[32*(Time+1)+78*(Time)+15+jt] = 2*Time*4+9;
     jCol[32*(Time+1)+78*(Time)+15+jt] = (Time+1)*2+jt;
   }
     for(Index jt=0;jt<Time+1;jt++) {
     iRow[33*(Time+1)+78*(Time)+15+jt] = 2*Time*4+9;
     jCol[33*(Time+1)+78*(Time)+15+jt] = (Time+1)*3+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[34*(Time+1)+78*(Time)+15+jt] = 2*Time*4+9;
     jCol[34*(Time+1)+78*(Time)+15+jt] = (Time+1)*4+Time*0+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[34*(Time+1)+79*(Time)+15+jt] = 2*Time*4+9;
     jCol[34*(Time+1)+79*(Time)+15+jt] = (Time+1)*4+Time*1+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[34*(Time+1)+80*(Time)+15+jt] = 2*Time*4+9;
     jCol[34*(Time+1)+80*(Time)+15+jt] = (Time+1)*4+Time*2+jt;
   }
     for(Index jt=0;jt<Time;jt++) {
     iRow[34*(Time+1)+81*(Time)+15+jt] = 2*Time*4+9;
     jCol[34*(Time+1)+81*(Time)+15+jt] = (Time+1)*4+Time*3+jt;
   }
     for(Index jt=0;jt<1;jt++) {
     iRow[34*(Time+1)+82*(Time)+15+jt] = 2*Time*4+9;
     jCol[34*(Time+1)+82*(Time)+15+jt] = 2*Time*4+4;
   }
     for(Index jt=0;jt<1;jt++) {
     iRow[34*(Time+1)+82*(Time)+16+jt] = 2*Time*4+9;
     jCol[34*(Time+1)+82*(Time)+16+jt] = 2*Time*4+5;
   }
     for(Index jt=0;jt<1;jt++) {
     iRow[34*(Time+1)+82*(Time)+17+jt] = 2*Time*4+9;
     jCol[34*(Time+1)+82*(Time)+17+jt] = 2*Time*4+6;
   }
     for(Index jt=0;jt<1;jt++) {
     iRow[34*(Time+1)+82*(Time)+18+jt] = 2*Time*4+9;
     jCol[34*(Time+1)+82*(Time)+18+jt] = 2*Time*4+7;
   }
     for(Index jt=0;jt<1;jt++) {
     iRow[34*(Time+1)+82*(Time)+19+jt] = 2*Time*4+9;
     jCol[34*(Time+1)+82*(Time)+19+jt] = 2*Time*4+8;
   }
     for(Index jt=0;jt<1;jt++) {
     iRow[34*(Time+1)+82*(Time)+20+jt] = 2*Time*4+9;
     jCol[34*(Time+1)+82*(Time)+20+jt] = 2*Time*4+9;
   }
}
else {
  // return the values.  This is a symmetric matrix, fill the lower left
  // triangle only
  // initialize the values array
  // Point to the initial starting spot for the Hessian elements

  for(Index jt=0;jt<34*(Time+1)+82*Time+21;jt++) values[jt] = 0.; // Initialize matrix

   // fill the objective portion

  for(Index jt=0;jt<Time;jt++) {

     for(Index i=0;i<nY;i++) {
        Xval[i] = x[jt + i*(Time+1)];
        Xvalp1[i] = x[jt + i*(Time+1) + 1];
        Xval2[i] = x[(Time+1)*(nY+nU) + jt + i*(Time)];
     } //end for loop

     for(Index cup=0;cup<nU;cup++) {
        K11val[cup] = x[jt + nY*(Time+1) + cup*(Time+1)];
        K11valp1[cup] = x[jt + nY*(Time+1) + cup*(Time+1) + 1];
        K11val2[cup] = x[(Time+1)*(nY+nU) + (nY+cup)*Time + jt];
     } //end for loop

     Xdval[0] = VDATA0[2*jt];
     Xdval2[0] = VDATA0[2*jt+1];
     Xdvalp1[0] = VDATA0[2*jt+2];
     Xdval[1] = MDATA0[2*jt];
     Xdval2[1] = MDATA0[2*jt+1];
     Xdvalp1[1] = MDATA0[2*jt+2];
     Xdval[2] = HDATA0[2*jt];
     Xdval2[2] = HDATA0[2*jt+1];
     Xdvalp1[2] = HDATA0[2*jt+2];
     Xdval[3] = NDATA0[2*jt];
     Xdval2[3] = NDATA0[2*jt+1];
     Xdvalp1[3] = NDATA0[2*jt+2];
     Ival[0] = Iinj[2*jt];
     Ival2[0] = Iinj[2*jt+1];
     Ivalp1[0] = Iinj[2*jt+2];

     for(Index i=0;i<nP;i++) {
        Pval[i] = x[(2*Time+1)*(nY+nU)+i];
     } //end for loop

    values[14*(Time+1)+62*(Time)+0] += obj_factor*bounds[0][3]*(0.03125*pow(hstep, 2)*pow((1.0L/4.0L)*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0]), 2) + 0.0555555555555556*pow(hstep, 2)*pow((1.0L/4.0L)*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0]), 2))/(2*Time+1);
    values[18*(Time+1)+66*(Time)+1] += obj_factor*bounds[0][3]*(0.03125*pow(hstep, 2)*((1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3))*((1.0L/4.0L)*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.0555555555555556*pow(hstep, 2)*((1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3) + Pval[0]*Xval2[2]*pow(Xval2[1], 3) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3))*((1.0L/4.0L)*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.25*hstep*((1.0L/4.0L)*Xval[2]*pow(Xval[1], 3) - 1.0L/4.0L*Xvalp1[2]*pow(Xvalp1[1], 3))*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]) + 0.333333333333333*hstep*((1.0L/4.0L)*Xval[2]*pow(Xval[1], 3) + Xval2[2]*pow(Xval2[1], 3) + (1.0L/4.0L)*Xvalp1[2]*pow(Xvalp1[1], 3))*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]))/(2*Time+1);
    values[18*(Time+1)+66*(Time)+2] += obj_factor*bounds[0][3]*(0.03125*pow(hstep, 2)*pow((1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3), 2) + 0.0555555555555556*pow(hstep, 2)*pow((1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3) + Pval[0]*Xval2[2]*pow(Xval2[1], 3) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3), 2))/(2*Time+1);
    values[22*(Time+1)+70*(Time)+3] += obj_factor*bounds[0][3]*(0.03125*pow(hstep, 2)*((1.0L/4.0L)*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]))*((1.0L/4.0L)*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.0555555555555556*pow(hstep, 2)*((1.0L/4.0L)*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]))*((1.0L/4.0L)*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])))/(2*Time+1);
    values[22*(Time+1)+70*(Time)+4] += obj_factor*bounds[0][3]*(0.03125*pow(hstep, 2)*((1.0L/4.0L)*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]))*((1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) + 0.0555555555555556*pow(hstep, 2)*((1.0L/4.0L)*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]))*((1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3) + Pval[0]*Xval2[2]*pow(Xval2[1], 3) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)))/(2*Time+1);
    values[22*(Time+1)+70*(Time)+5] += obj_factor*bounds[0][3]*(0.03125*pow(hstep, 2)*pow((1.0L/4.0L)*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]), 2) + 0.0555555555555556*pow(hstep, 2)*pow((1.0L/4.0L)*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]), 2))/(2*Time+1);
    values[26*(Time+1)+74*(Time)+6] += obj_factor*bounds[0][3]*(0.03125*pow(hstep, 2)*((1.0L/4.0L)*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4))*((1.0L/4.0L)*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.0555555555555556*pow(hstep, 2)*((1.0L/4.0L)*Pval[2]*pow(Xval[3], 4) + Pval[2]*pow(Xval2[3], 4) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4))*((1.0L/4.0L)*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])))/(2*Time+1);
    values[26*(Time+1)+74*(Time)+7] += obj_factor*bounds[0][3]*(0.03125*pow(hstep, 2)*((1.0L/4.0L)*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4))*((1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) + 0.0555555555555556*pow(hstep, 2)*((1.0L/4.0L)*Pval[2]*pow(Xval[3], 4) + Pval[2]*pow(Xval2[3], 4) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4))*((1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3) + Pval[0]*Xval2[2]*pow(Xval2[1], 3) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)))/(2*Time+1);
    values[26*(Time+1)+74*(Time)+8] += obj_factor*bounds[0][3]*(0.03125*pow(hstep, 2)*((1.0L/4.0L)*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4))*((1.0L/4.0L)*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0])) + 0.0555555555555556*pow(hstep, 2)*((1.0L/4.0L)*Pval[2]*pow(Xval[3], 4) + Pval[2]*pow(Xval2[3], 4) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4))*((1.0L/4.0L)*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0])) + 0.25*hstep*((1.0L/4.0L)*pow(Xval[3], 4) - 1.0L/4.0L*pow(Xvalp1[3], 4))*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]) + 0.333333333333333*hstep*((1.0L/4.0L)*pow(Xval[3], 4) + pow(Xval2[3], 4) + (1.0L/4.0L)*pow(Xvalp1[3], 4))*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]))/(2*Time+1);
    values[26*(Time+1)+74*(Time)+9] += obj_factor*bounds[0][3]*(0.03125*pow(hstep, 2)*pow((1.0L/4.0L)*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4), 2) + 0.0555555555555556*pow(hstep, 2)*pow((1.0L/4.0L)*Pval[2]*pow(Xval[3], 4) + Pval[2]*pow(Xval2[3], 4) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4), 2))/(2*Time+1);
    values[30*(Time+1)+78*(Time)+10] += obj_factor*bounds[0][3]*(0.03125*pow(hstep, 2)*(-1.0L/4.0L*Xval[0] + (1.0L/4.0L)*Xvalp1[0])*((1.0L/4.0L)*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.0555555555555556*pow(hstep, 2)*((1.0L/4.0L)*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0]))*((3.0L/2.0L)*Pval[5] - 1.0L/4.0L*Xval[0] - Xval2[0] - 1.0L/4.0L*Xvalp1[0]))/(2*Time+1);
    values[30*(Time+1)+78*(Time)+11] += obj_factor*bounds[0][3]*(0.03125*pow(hstep, 2)*(-1.0L/4.0L*Xval[0] + (1.0L/4.0L)*Xvalp1[0])*((1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) + 0.0555555555555556*pow(hstep, 2)*((1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3) + Pval[0]*Xval2[2]*pow(Xval2[1], 3) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3))*((3.0L/2.0L)*Pval[5] - 1.0L/4.0L*Xval[0] - Xval2[0] - 1.0L/4.0L*Xvalp1[0]))/(2*Time+1);
    values[30*(Time+1)+78*(Time)+12] += obj_factor*bounds[0][3]*(0.03125*pow(hstep, 2)*(-1.0L/4.0L*Xval[0] + (1.0L/4.0L)*Xvalp1[0])*((1.0L/4.0L)*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0])) + 0.0555555555555556*pow(hstep, 2)*((1.0L/4.0L)*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]))*((3.0L/2.0L)*Pval[5] - 1.0L/4.0L*Xval[0] - Xval2[0] - 1.0L/4.0L*Xvalp1[0]))/(2*Time+1);
    values[30*(Time+1)+78*(Time)+13] += obj_factor*bounds[0][3]*(0.03125*pow(hstep, 2)*(-1.0L/4.0L*Xval[0] + (1.0L/4.0L)*Xvalp1[0])*((1.0L/4.0L)*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)) + 0.0555555555555556*pow(hstep, 2)*((1.0L/4.0L)*Pval[2]*pow(Xval[3], 4) + Pval[2]*pow(Xval2[3], 4) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4))*((3.0L/2.0L)*Pval[5] - 1.0L/4.0L*Xval[0] - Xval2[0] - 1.0L/4.0L*Xvalp1[0]))/(2*Time+1);
    values[30*(Time+1)+78*(Time)+14] += obj_factor*bounds[0][3]*(0.03125*pow(hstep, 2)*pow(-1.0L/4.0L*Xval[0] + (1.0L/4.0L)*Xvalp1[0], 2) + 0.0555555555555556*pow(hstep, 2)*pow((3.0L/2.0L)*Pval[5] - 1.0L/4.0L*Xval[0] - Xval2[0] - 1.0L/4.0L*Xvalp1[0], 2))/(2*Time+1);
    values[34*(Time+1)+82*(Time)+15] += obj_factor*bounds[0][3]*(0.0833333333333333*Pval[4]*pow(hstep, 2)*((1.0L/4.0L)*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])))/(2*Time+1);
    values[34*(Time+1)+82*(Time)+16] += obj_factor*bounds[0][3]*(0.0833333333333333*Pval[4]*pow(hstep, 2)*((1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3) + Pval[0]*Xval2[2]*pow(Xval2[1], 3) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)))/(2*Time+1);
    values[34*(Time+1)+82*(Time)+17] += obj_factor*bounds[0][3]*(0.0833333333333333*Pval[4]*pow(hstep, 2)*((1.0L/4.0L)*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0])))/(2*Time+1);
    values[34*(Time+1)+82*(Time)+18] += obj_factor*bounds[0][3]*(0.0833333333333333*Pval[4]*pow(hstep, 2)*((1.0L/4.0L)*Pval[2]*pow(Xval[3], 4) + Pval[2]*pow(Xval2[3], 4) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)))/(2*Time+1);
    values[34*(Time+1)+82*(Time)+19] += obj_factor*bounds[0][3]*(0.0833333333333333*Pval[4]*pow(hstep, 2)*((3.0L/2.0L)*Pval[5] - 1.0L/4.0L*Xval[0] - Xval2[0] - 1.0L/4.0L*Xvalp1[0]) + 0.5*hstep*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]))/(2*Time+1);
    values[34*(Time+1)+82*(Time)+20] += obj_factor*bounds[0][3]*(0.125*pow(Pval[4], 2)*pow(hstep, 2))/(2*Time+1);
   values[0*(Time+1)+0*(Time)+0+jt] += obj_factor*1*(2)/(2*Time+1);
   values[2*(Time+1)+2*(Time)+0+jt] += obj_factor*1*(2)/(2*Time+1);
   values[5*(Time+1)+6*(Time)+0+jt] += obj_factor*1*(2)/(2*Time+1);
   values[9*(Time+1)+12*(Time)+0+jt] += obj_factor*1*(2)/(2*Time+1);
   values[10*(Time+1)+24*(Time)+0+jt] += obj_factor*1*(2)/(2*Time+1);
   values[10*(Time+1)+34*(Time)+0+jt] += obj_factor*1*(2)/(2*Time+1);
   values[10*(Time+1)+45*(Time)+0+jt] += obj_factor*1*(2)/(2*Time+1);
   values[10*(Time+1)+57*(Time)+0+jt] += obj_factor*1*(2)/(2*Time+1);
   values[0*(Time+1)+0*(Time)+0+jt] += obj_factor*bounds[0][3]*((0.125*hstep*(-1.0L/4.0L*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xval[2]*pow(Xval[1], 3)) + 0.5)*(0.25*hstep*(-1.0L/4.0L*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xval[2]*pow(Xval[1], 3)) + 1.0) + (0.166666666666667*hstep*(-1.0L/4.0L*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xval[2]*pow(Xval[1], 3)) + 1)*(0.333333333333333*hstep*(-1.0L/4.0L*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xval[2]*pow(Xval[1], 3)) + 2))/(2*Time+1);
   values[1*(Time+1)+0*(Time)+0+jt] += obj_factor*bounds[0][3]*((0.125*hstep*(-1.0L/4.0L*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xval[2]*pow(Xval[1], 3)) + 0.5)*(0.25*hstep*((1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4) + (1.0L/4.0L)*Pval[4] + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) + 1.0) + (0.166666666666667*hstep*(-1.0L/4.0L*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xval[2]*pow(Xval[1], 3)) + 1)*(0.333333333333333*hstep*(-1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) - 2))/(2*Time+1);
   values[0*(Time+1)+0*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*((0.166666666666667*hstep*(-1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) - 1)*(0.333333333333333*hstep*(-1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) - 2) + (0.125*hstep*((1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4) + (1.0L/4.0L)*Pval[4] + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) + 0.5)*(0.25*hstep*((1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4) + (1.0L/4.0L)*Pval[4] + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) + 1.0))/(2*Time+1);
   values[1*(Time+1)+1*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.1875*Pval[0]*Xval[2]*hstep*pow(Xval[1], 2)*(Pval[1] - Xval[0])*(0.125*hstep*(-1.0L/4.0L*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xval[2]*pow(Xval[1], 3)) + 0.5) + 0.25*Pval[0]*Xval[2]*hstep*pow(Xval[1], 2)*(Pval[1] - Xval[0])*(0.166666666666667*hstep*(-1.0L/4.0L*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xval[2]*pow(Xval[1], 3)) + 1) - 0.25*Pval[0]*Xval[2]*hstep*pow(Xval[1], 2)*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]) - 0.1875*Pval[0]*Xval[2]*hstep*pow(Xval[1], 2)*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
   values[2*(Time+1)+1*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.25*Pval[0]*Xval[2]*hstep*pow(Xval[1], 2)*(Pval[1] - Xval[0])*(0.166666666666667*hstep*(-1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) - 1) + 0.1875*Pval[0]*Xval[2]*hstep*pow(Xval[1], 2)*(Pval[1] - Xval[0])*(0.125*hstep*((1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4) + (1.0L/4.0L)*Pval[4] + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) + 0.5))/(2*Time+1);
   values[2*(Time+1)+2*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.048828125*pow(Pval[0], 2)*pow(Xval[2], 2)*pow(hstep, 2)*pow(Xval[1], 4)*pow(Pval[1] - Xval[0], 2) + 0.5*Pval[0]*Xval[2]*hstep*Xval[1]*(Pval[1] - Xval[0])*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]) + 0.375*Pval[0]*Xval[2]*hstep*Xval[1]*(Pval[1] - Xval[0])*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
   values[3*(Time+1)+2*(Time)+0+jt] += obj_factor*bounds[0][3]*(-0.1875*Pval[0]*Xvalp1[2]*hstep*pow(Xvalp1[1], 2)*(Pval[1] - Xvalp1[0])*(0.125*hstep*(-1.0L/4.0L*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xval[2]*pow(Xval[1], 3)) + 0.5) + 0.25*Pval[0]*Xvalp1[2]*hstep*pow(Xvalp1[1], 2)*(Pval[1] - Xvalp1[0])*(0.166666666666667*hstep*(-1.0L/4.0L*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xval[2]*pow(Xval[1], 3)) + 1))/(2*Time+1);
   values[1*(Time+1)+1*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(0.25*Pval[0]*Xvalp1[2]*hstep*pow(Xvalp1[1], 2)*(Pval[1] - Xvalp1[0])*(0.166666666666667*hstep*(-1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) - 1) - 0.1875*Pval[0]*Xvalp1[2]*hstep*pow(Xvalp1[1], 2)*(Pval[1] - Xvalp1[0])*(0.125*hstep*((1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4) + (1.0L/4.0L)*Pval[4] + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) + 0.5) - 0.25*Pval[0]*Xvalp1[2]*hstep*pow(Xvalp1[1], 2)*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]) + 0.1875*Pval[0]*Xvalp1[2]*hstep*pow(Xvalp1[1], 2)*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
   values[3*(Time+1)+3*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.013671875*pow(Pval[0], 2)*Xval[2]*Xvalp1[2]*pow(hstep, 2)*pow(Xval[1], 2)*pow(Xvalp1[1], 2)*(Pval[1] - Xval[0])*(Pval[1] - Xvalp1[0]))/(2*Time+1);
   values[2*(Time+1)+2*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(0.048828125*pow(Pval[0], 2)*pow(Xvalp1[2], 2)*pow(hstep, 2)*pow(Xvalp1[1], 4)*pow(Pval[1] - Xvalp1[0], 2) + 0.5*Pval[0]*Xvalp1[2]*hstep*Xvalp1[1]*(Pval[1] - Xvalp1[0])*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]) - 0.375*Pval[0]*Xvalp1[2]*hstep*Xvalp1[1]*(Pval[1] - Xvalp1[0])*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
   values[3*(Time+1)+4*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0625*Pval[0]*hstep*pow(Xval[1], 3)*(Pval[1] - Xval[0])*(0.125*hstep*(-1.0L/4.0L*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xval[2]*pow(Xval[1], 3)) + 0.5) + 0.0833333333333333*Pval[0]*hstep*pow(Xval[1], 3)*(Pval[1] - Xval[0])*(0.166666666666667*hstep*(-1.0L/4.0L*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xval[2]*pow(Xval[1], 3)) + 1) - 0.0833333333333333*Pval[0]*hstep*pow(Xval[1], 3)*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]) - 0.0625*Pval[0]*hstep*pow(Xval[1], 3)*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
   values[4*(Time+1)+4*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0833333333333333*Pval[0]*hstep*pow(Xval[1], 3)*(Pval[1] - Xval[0])*(0.166666666666667*hstep*(-1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) - 1) + 0.0625*Pval[0]*hstep*pow(Xval[1], 3)*(Pval[1] - Xval[0])*(0.125*hstep*((1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4) + (1.0L/4.0L)*Pval[4] + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) + 0.5))/(2*Time+1);
   values[4*(Time+1)+5*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0162760416666667*pow(Pval[0], 2)*Xval[2]*pow(hstep, 2)*pow(Xval[1], 5)*pow(Pval[1] - Xval[0], 2) + 0.25*Pval[0]*hstep*pow(Xval[1], 2)*(Pval[1] - Xval[0])*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]) + 0.1875*Pval[0]*hstep*pow(Xval[1], 2)*(Pval[1] - Xval[0])*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
   values[5*(Time+1)+5*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.00455729166666667*pow(Pval[0], 2)*Xvalp1[2]*pow(hstep, 2)*pow(Xval[1], 3)*pow(Xvalp1[1], 2)*(Pval[1] - Xval[0])*(Pval[1] - Xvalp1[0]))/(2*Time+1);
   values[5*(Time+1)+6*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.00542534722222222*pow(Pval[0], 2)*pow(hstep, 2)*pow(Xval[1], 6)*pow(Pval[1] - Xval[0], 2))/(2*Time+1);
   values[6*(Time+1)+6*(Time)+0+jt] += obj_factor*bounds[0][3]*(-0.0625*Pval[0]*hstep*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])*(0.125*hstep*(-1.0L/4.0L*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xval[2]*pow(Xval[1], 3)) + 0.5) + 0.0833333333333333*Pval[0]*hstep*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])*(0.166666666666667*hstep*(-1.0L/4.0L*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xval[2]*pow(Xval[1], 3)) + 1))/(2*Time+1);
   values[3*(Time+1)+4*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(0.0833333333333333*Pval[0]*hstep*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])*(0.166666666666667*hstep*(-1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) - 1) - 0.0625*Pval[0]*hstep*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])*(0.125*hstep*((1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4) + (1.0L/4.0L)*Pval[4] + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) + 0.5) - 0.0833333333333333*Pval[0]*hstep*pow(Xvalp1[1], 3)*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]) + 0.0625*Pval[0]*hstep*pow(Xvalp1[1], 3)*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
   values[6*(Time+1)+7*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.00455729166666667*pow(Pval[0], 2)*Xval[2]*pow(hstep, 2)*pow(Xval[1], 2)*pow(Xvalp1[1], 3)*(Pval[1] - Xval[0])*(Pval[1] - Xvalp1[0]))/(2*Time+1);
   values[4*(Time+1)+5*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(0.0162760416666667*pow(Pval[0], 2)*Xvalp1[2]*pow(hstep, 2)*pow(Xvalp1[1], 5)*pow(Pval[1] - Xvalp1[0], 2) + 0.25*Pval[0]*hstep*pow(Xvalp1[1], 2)*(Pval[1] - Xvalp1[0])*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]) - 0.1875*Pval[0]*hstep*pow(Xvalp1[1], 2)*(Pval[1] - Xvalp1[0])*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
   values[6*(Time+1)+8*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.00151909722222222*pow(Pval[0], 2)*pow(hstep, 2)*pow(Xval[1], 3)*pow(Xvalp1[1], 3)*(Pval[1] - Xval[0])*(Pval[1] - Xvalp1[0]))/(2*Time+1);
   values[5*(Time+1)+6*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(0.00542534722222222*pow(Pval[0], 2)*pow(hstep, 2)*pow(Xvalp1[1], 6)*pow(Pval[1] - Xvalp1[0], 2))/(2*Time+1);
   values[6*(Time+1)+9*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.25*Pval[2]*hstep*pow(Xval[3], 3)*(Pval[3] - Xval[0])*(0.125*hstep*(-1.0L/4.0L*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xval[2]*pow(Xval[1], 3)) + 0.5) + 0.333333333333333*Pval[2]*hstep*pow(Xval[3], 3)*(Pval[3] - Xval[0])*(0.166666666666667*hstep*(-1.0L/4.0L*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xval[2]*pow(Xval[1], 3)) + 1) - 0.333333333333333*Pval[2]*hstep*pow(Xval[3], 3)*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]) - 0.25*Pval[2]*hstep*pow(Xval[3], 3)*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
   values[7*(Time+1)+9*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.333333333333333*Pval[2]*hstep*pow(Xval[3], 3)*(Pval[3] - Xval[0])*(0.166666666666667*hstep*(-1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) - 1) + 0.25*Pval[2]*hstep*pow(Xval[3], 3)*(Pval[3] - Xval[0])*(0.125*hstep*((1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4) + (1.0L/4.0L)*Pval[4] + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) + 0.5))/(2*Time+1);
   values[7*(Time+1)+10*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0651041666666667*Pval[2]*Pval[0]*Xval[2]*pow(hstep, 2)*pow(Xval[1], 2)*pow(Xval[3], 3)*(Pval[3] - Xval[0])*(Pval[1] - Xval[0]))/(2*Time+1);
   values[8*(Time+1)+10*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0182291666666667*Pval[2]*Pval[0]*Xvalp1[2]*pow(hstep, 2)*pow(Xvalp1[1], 2)*pow(Xval[3], 3)*(Pval[3] - Xval[0])*(Pval[1] - Xvalp1[0]))/(2*Time+1);
   values[8*(Time+1)+11*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0217013888888889*Pval[2]*Pval[0]*pow(hstep, 2)*pow(Xval[1], 3)*pow(Xval[3], 3)*(Pval[3] - Xval[0])*(Pval[1] - Xval[0]))/(2*Time+1);
   values[9*(Time+1)+11*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.00607638888888889*Pval[2]*Pval[0]*pow(hstep, 2)*pow(Xvalp1[1], 3)*pow(Xval[3], 3)*(Pval[3] - Xval[0])*(Pval[1] - Xvalp1[0]))/(2*Time+1);
   values[9*(Time+1)+12*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0868055555555556*pow(Pval[2], 2)*pow(hstep, 2)*pow(Xval[3], 6)*pow(Pval[3] - Xval[0], 2) + 1.0*Pval[2]*hstep*pow(Xval[3], 2)*(Pval[3] - Xval[0])*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]) + 0.75*Pval[2]*hstep*pow(Xval[3], 2)*(Pval[3] - Xval[0])*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
   values[10*(Time+1)+12*(Time)+0+jt] += obj_factor*bounds[0][3]*(-0.25*Pval[2]*hstep*pow(Xvalp1[3], 3)*(Pval[3] - Xvalp1[0])*(0.125*hstep*(-1.0L/4.0L*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xval[2]*pow(Xval[1], 3)) + 0.5) + 0.333333333333333*Pval[2]*hstep*pow(Xvalp1[3], 3)*(Pval[3] - Xvalp1[0])*(0.166666666666667*hstep*(-1.0L/4.0L*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xval[2]*pow(Xval[1], 3)) + 1))/(2*Time+1);
   values[6*(Time+1)+9*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(0.333333333333333*Pval[2]*hstep*pow(Xvalp1[3], 3)*(Pval[3] - Xvalp1[0])*(0.166666666666667*hstep*(-1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) - 1) - 0.25*Pval[2]*hstep*pow(Xvalp1[3], 3)*(Pval[3] - Xvalp1[0])*(0.125*hstep*((1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4) + (1.0L/4.0L)*Pval[4] + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) + 0.5) - 0.333333333333333*Pval[2]*hstep*pow(Xvalp1[3], 3)*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]) + 0.25*Pval[2]*hstep*pow(Xvalp1[3], 3)*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
   values[10*(Time+1)+13*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0182291666666667*Pval[2]*Pval[0]*Xval[2]*pow(hstep, 2)*pow(Xval[1], 2)*pow(Xvalp1[3], 3)*(Pval[3] - Xvalp1[0])*(Pval[1] - Xval[0]))/(2*Time+1);
   values[7*(Time+1)+10*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(0.0651041666666667*Pval[2]*Pval[0]*Xvalp1[2]*pow(hstep, 2)*pow(Xvalp1[1], 2)*pow(Xvalp1[3], 3)*(Pval[3] - Xvalp1[0])*(Pval[1] - Xvalp1[0]))/(2*Time+1);
   values[10*(Time+1)+14*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.00607638888888889*Pval[2]*Pval[0]*pow(hstep, 2)*pow(Xval[1], 3)*pow(Xvalp1[3], 3)*(Pval[3] - Xvalp1[0])*(Pval[1] - Xval[0]))/(2*Time+1);
   values[8*(Time+1)+11*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(0.0217013888888889*Pval[2]*Pval[0]*pow(hstep, 2)*pow(Xvalp1[1], 3)*pow(Xvalp1[3], 3)*(Pval[3] - Xvalp1[0])*(Pval[1] - Xvalp1[0]))/(2*Time+1);
   values[10*(Time+1)+15*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0243055555555556*pow(Pval[2], 2)*pow(hstep, 2)*pow(Xval[3], 3)*pow(Xvalp1[3], 3)*(Pval[3] - Xval[0])*(Pval[3] - Xvalp1[0]))/(2*Time+1);
   values[9*(Time+1)+12*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(0.0868055555555556*pow(Pval[2], 2)*pow(hstep, 2)*pow(Xvalp1[3], 6)*pow(Pval[3] - Xvalp1[0], 2) + 1.0*Pval[2]*hstep*pow(Xvalp1[3], 2)*(Pval[3] - Xvalp1[0])*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]) - 0.75*Pval[2]*hstep*pow(Xvalp1[3], 2)*(Pval[3] - Xvalp1[0])*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
   values[10*(Time+1)+16*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.333333333333333*hstep*(0.166666666666667*hstep*(-1.0L/4.0L*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xval[2]*pow(Xval[1], 3)) + 1)*(-Pval[2]*pow(Xval2[3], 4) - Pval[4] - Pval[0]*Xval2[2]*pow(Xval2[1], 3)) - 0.25*hstep*(-1.0L/4.0L*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xval[2]*pow(Xval[1], 3)) - 1.0)/(2*Time+1);
   values[10*(Time+1)+17*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.333333333333333*hstep*(0.166666666666667*hstep*(-1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) - 1)*(-Pval[2]*pow(Xval2[3], 4) - Pval[4] - Pval[0]*Xval2[2]*pow(Xval2[1], 3)) - 0.25*hstep*((1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4) + (1.0L/4.0L)*Pval[4] + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) - 1.0)/(2*Time+1);
   values[10*(Time+1)+18*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0416666666666667*Pval[0]*Xval[2]*pow(hstep, 2)*pow(Xval[1], 2)*(Pval[1] - Xval[0])*(-Pval[2]*pow(Xval2[3], 4) - Pval[4] - Pval[0]*Xval2[2]*pow(Xval2[1], 3)) - 0.1875*Pval[0]*Xval[2]*hstep*pow(Xval[1], 2)*(Pval[1] - Xval[0]))/(2*Time+1);
   values[10*(Time+1)+19*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0416666666666667*Pval[0]*Xvalp1[2]*pow(hstep, 2)*pow(Xvalp1[1], 2)*(Pval[1] - Xvalp1[0])*(-Pval[2]*pow(Xval2[3], 4) - Pval[4] - Pval[0]*Xval2[2]*pow(Xval2[1], 3)) + 0.1875*Pval[0]*Xvalp1[2]*hstep*pow(Xvalp1[1], 2)*(Pval[1] - Xvalp1[0]))/(2*Time+1);
   values[10*(Time+1)+20*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0138888888888889*Pval[0]*pow(hstep, 2)*pow(Xval[1], 3)*(Pval[1] - Xval[0])*(-Pval[2]*pow(Xval2[3], 4) - Pval[4] - Pval[0]*Xval2[2]*pow(Xval2[1], 3)) - 0.0625*Pval[0]*hstep*pow(Xval[1], 3)*(Pval[1] - Xval[0]))/(2*Time+1);
   values[10*(Time+1)+21*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0138888888888889*Pval[0]*pow(hstep, 2)*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])*(-Pval[2]*pow(Xval2[3], 4) - Pval[4] - Pval[0]*Xval2[2]*pow(Xval2[1], 3)) + 0.0625*Pval[0]*hstep*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0]))/(2*Time+1);
   values[10*(Time+1)+22*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0555555555555556*Pval[2]*pow(hstep, 2)*pow(Xval[3], 3)*(Pval[3] - Xval[0])*(-Pval[2]*pow(Xval2[3], 4) - Pval[4] - Pval[0]*Xval2[2]*pow(Xval2[1], 3)) - 0.25*Pval[2]*hstep*pow(Xval[3], 3)*(Pval[3] - Xval[0]))/(2*Time+1);
   values[10*(Time+1)+23*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0555555555555556*Pval[2]*pow(hstep, 2)*pow(Xvalp1[3], 3)*(Pval[3] - Xvalp1[0])*(-Pval[2]*pow(Xval2[3], 4) - Pval[4] - Pval[0]*Xval2[2]*pow(Xval2[1], 3)) + 0.25*Pval[2]*hstep*pow(Xvalp1[3], 3)*(Pval[3] - Xvalp1[0]))/(2*Time+1);
   values[10*(Time+1)+24*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0555555555555556*pow(hstep, 2)*pow(-Pval[2]*pow(Xval2[3], 4) - Pval[4] - Pval[0]*Xval2[2]*pow(Xval2[1], 3), 2) + 2)/(2*Time+1);
   values[10*(Time+1)+25*(Time)+0+jt] += obj_factor*bounds[0][3]*(1.0*Pval[0]*Xval2[2]*hstep*pow(Xval2[1], 2)*(Pval[1] - Xval2[0])*(0.166666666666667*hstep*(-1.0L/4.0L*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xval[2]*pow(Xval[1], 3)) + 1))/(2*Time+1);
   values[10*(Time+1)+26*(Time)+0+jt] += obj_factor*bounds[0][3]*(1.0*Pval[0]*Xval2[2]*hstep*pow(Xval2[1], 2)*(Pval[1] - Xval2[0])*(0.166666666666667*hstep*(-1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) - 1))/(2*Time+1);
   values[10*(Time+1)+27*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.125*pow(Pval[0], 2)*Xval[2]*Xval2[2]*pow(hstep, 2)*pow(Xval[1], 2)*pow(Xval2[1], 2)*(Pval[1] - Xval[0])*(Pval[1] - Xval2[0]))/(2*Time+1);
   values[10*(Time+1)+28*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.125*pow(Pval[0], 2)*Xval2[2]*Xvalp1[2]*pow(hstep, 2)*pow(Xval2[1], 2)*pow(Xvalp1[1], 2)*(Pval[1] - Xval2[0])*(Pval[1] - Xvalp1[0]))/(2*Time+1);
   values[10*(Time+1)+29*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0416666666666667*pow(Pval[0], 2)*Xval2[2]*pow(hstep, 2)*pow(Xval[1], 3)*pow(Xval2[1], 2)*(Pval[1] - Xval[0])*(Pval[1] - Xval2[0]))/(2*Time+1);
   values[10*(Time+1)+30*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0416666666666667*pow(Pval[0], 2)*Xval2[2]*pow(hstep, 2)*pow(Xval2[1], 2)*pow(Xvalp1[1], 3)*(Pval[1] - Xval2[0])*(Pval[1] - Xvalp1[0]))/(2*Time+1);
   values[10*(Time+1)+31*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.166666666666667*Pval[2]*Pval[0]*Xval2[2]*pow(hstep, 2)*pow(Xval2[1], 2)*pow(Xval[3], 3)*(Pval[3] - Xval[0])*(Pval[1] - Xval2[0]))/(2*Time+1);
   values[10*(Time+1)+32*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.166666666666667*Pval[2]*Pval[0]*Xval2[2]*pow(hstep, 2)*pow(Xval2[1], 2)*pow(Xvalp1[3], 3)*(Pval[3] - Xvalp1[0])*(Pval[1] - Xval2[0]))/(2*Time+1);
   values[10*(Time+1)+33*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.166666666666667*Pval[0]*Xval2[2]*pow(hstep, 2)*pow(Xval2[1], 2)*(Pval[1] - Xval2[0])*(-Pval[2]*pow(Xval2[3], 4) - Pval[4] - Pval[0]*Xval2[2]*pow(Xval2[1], 3)) - 1.0*Pval[0]*Xval2[2]*hstep*pow(Xval2[1], 2)*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]))/(2*Time+1);
   values[10*(Time+1)+34*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.5*pow(Pval[0], 2)*pow(Xval2[2], 2)*pow(hstep, 2)*pow(Xval2[1], 4)*pow(Pval[1] - Xval2[0], 2) + 2.0*Pval[0]*Xval2[2]*hstep*Xval2[1]*(Pval[1] - Xval2[0])*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]))/(2*Time+1);
   values[10*(Time+1)+35*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.333333333333333*Pval[0]*hstep*pow(Xval2[1], 3)*(Pval[1] - Xval2[0])*(0.166666666666667*hstep*(-1.0L/4.0L*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xval[2]*pow(Xval[1], 3)) + 1))/(2*Time+1);
   values[10*(Time+1)+36*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.333333333333333*Pval[0]*hstep*pow(Xval2[1], 3)*(Pval[1] - Xval2[0])*(0.166666666666667*hstep*(-1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) - 1))/(2*Time+1);
   values[10*(Time+1)+37*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0416666666666667*pow(Pval[0], 2)*Xval[2]*pow(hstep, 2)*pow(Xval[1], 2)*pow(Xval2[1], 3)*(Pval[1] - Xval[0])*(Pval[1] - Xval2[0]))/(2*Time+1);
   values[10*(Time+1)+38*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0416666666666667*pow(Pval[0], 2)*Xvalp1[2]*pow(hstep, 2)*pow(Xval2[1], 3)*pow(Xvalp1[1], 2)*(Pval[1] - Xval2[0])*(Pval[1] - Xvalp1[0]))/(2*Time+1);
   values[10*(Time+1)+39*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0138888888888889*pow(Pval[0], 2)*pow(hstep, 2)*pow(Xval[1], 3)*pow(Xval2[1], 3)*(Pval[1] - Xval[0])*(Pval[1] - Xval2[0]))/(2*Time+1);
   values[10*(Time+1)+40*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0138888888888889*pow(Pval[0], 2)*pow(hstep, 2)*pow(Xval2[1], 3)*pow(Xvalp1[1], 3)*(Pval[1] - Xval2[0])*(Pval[1] - Xvalp1[0]))/(2*Time+1);
   values[10*(Time+1)+41*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0555555555555556*Pval[2]*Pval[0]*pow(hstep, 2)*pow(Xval2[1], 3)*pow(Xval[3], 3)*(Pval[3] - Xval[0])*(Pval[1] - Xval2[0]))/(2*Time+1);
   values[10*(Time+1)+42*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0555555555555556*Pval[2]*Pval[0]*pow(hstep, 2)*pow(Xval2[1], 3)*pow(Xvalp1[3], 3)*(Pval[3] - Xvalp1[0])*(Pval[1] - Xval2[0]))/(2*Time+1);
   values[10*(Time+1)+43*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0555555555555556*Pval[0]*pow(hstep, 2)*pow(Xval2[1], 3)*(Pval[1] - Xval2[0])*(-Pval[2]*pow(Xval2[3], 4) - Pval[4] - Pval[0]*Xval2[2]*pow(Xval2[1], 3)) - 0.333333333333333*Pval[0]*hstep*pow(Xval2[1], 3)*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]))/(2*Time+1);
   values[10*(Time+1)+44*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.166666666666667*pow(Pval[0], 2)*Xval2[2]*pow(hstep, 2)*pow(Xval2[1], 5)*pow(Pval[1] - Xval2[0], 2) + 1.0*Pval[0]*hstep*pow(Xval2[1], 2)*(Pval[1] - Xval2[0])*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]))/(2*Time+1);
   values[10*(Time+1)+45*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0555555555555556*pow(Pval[0], 2)*pow(hstep, 2)*pow(Xval2[1], 6)*pow(Pval[1] - Xval2[0], 2))/(2*Time+1);
   values[10*(Time+1)+46*(Time)+0+jt] += obj_factor*bounds[0][3]*(1.33333333333333*Pval[2]*hstep*pow(Xval2[3], 3)*(Pval[3] - Xval2[0])*(0.166666666666667*hstep*(-1.0L/4.0L*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xval[2]*pow(Xval[1], 3)) + 1))/(2*Time+1);
   values[10*(Time+1)+47*(Time)+0+jt] += obj_factor*bounds[0][3]*(1.33333333333333*Pval[2]*hstep*pow(Xval2[3], 3)*(Pval[3] - Xval2[0])*(0.166666666666667*hstep*(-1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) - 1))/(2*Time+1);
   values[10*(Time+1)+48*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.166666666666667*Pval[2]*Pval[0]*Xval[2]*pow(hstep, 2)*pow(Xval[1], 2)*pow(Xval2[3], 3)*(Pval[3] - Xval2[0])*(Pval[1] - Xval[0]))/(2*Time+1);
   values[10*(Time+1)+49*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.166666666666667*Pval[2]*Pval[0]*Xvalp1[2]*pow(hstep, 2)*pow(Xvalp1[1], 2)*pow(Xval2[3], 3)*(Pval[3] - Xval2[0])*(Pval[1] - Xvalp1[0]))/(2*Time+1);
   values[10*(Time+1)+50*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0555555555555556*Pval[2]*Pval[0]*pow(hstep, 2)*pow(Xval[1], 3)*pow(Xval2[3], 3)*(Pval[3] - Xval2[0])*(Pval[1] - Xval[0]))/(2*Time+1);
   values[10*(Time+1)+51*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0555555555555556*Pval[2]*Pval[0]*pow(hstep, 2)*pow(Xvalp1[1], 3)*pow(Xval2[3], 3)*(Pval[3] - Xval2[0])*(Pval[1] - Xvalp1[0]))/(2*Time+1);
   values[10*(Time+1)+52*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.222222222222222*pow(Pval[2], 2)*pow(hstep, 2)*pow(Xval[3], 3)*pow(Xval2[3], 3)*(Pval[3] - Xval[0])*(Pval[3] - Xval2[0]))/(2*Time+1);
   values[10*(Time+1)+53*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.222222222222222*pow(Pval[2], 2)*pow(hstep, 2)*pow(Xval2[3], 3)*pow(Xvalp1[3], 3)*(Pval[3] - Xval2[0])*(Pval[3] - Xvalp1[0]))/(2*Time+1);
   values[10*(Time+1)+54*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.222222222222222*Pval[2]*pow(hstep, 2)*pow(Xval2[3], 3)*(Pval[3] - Xval2[0])*(-Pval[2]*pow(Xval2[3], 4) - Pval[4] - Pval[0]*Xval2[2]*pow(Xval2[1], 3)) - 1.33333333333333*Pval[2]*hstep*pow(Xval2[3], 3)*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]))/(2*Time+1);
   values[10*(Time+1)+55*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.666666666666667*Pval[2]*Pval[0]*Xval2[2]*pow(hstep, 2)*pow(Xval2[1], 2)*pow(Xval2[3], 3)*(Pval[3] - Xval2[0])*(Pval[1] - Xval2[0]))/(2*Time+1);
   values[10*(Time+1)+56*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.222222222222222*Pval[2]*Pval[0]*pow(hstep, 2)*pow(Xval2[1], 3)*pow(Xval2[3], 3)*(Pval[3] - Xval2[0])*(Pval[1] - Xval2[0]))/(2*Time+1);
   values[10*(Time+1)+57*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.888888888888889*pow(Pval[2], 2)*pow(hstep, 2)*pow(Xval2[3], 6)*pow(Pval[3] - Xval2[0], 2) + 4.0*Pval[2]*hstep*pow(Xval2[3], 2)*(Pval[3] - Xval2[0])*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]))/(2*Time+1);
   values[10*(Time+1)+58*(Time)+0+jt] += obj_factor*bounds[0][3]*(-0.0833333333333333*Xval[2]*hstep*pow(Xval[1], 3)*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]) - 0.0625*Xval[2]*hstep*pow(Xval[1], 3)*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]) + 0.25*hstep*(0.125*hstep*(-1.0L/4.0L*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xval[2]*pow(Xval[1], 3)) + 0.5)*((1.0L/4.0L)*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.333333333333333*hstep*(0.166666666666667*hstep*(-1.0L/4.0L*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xval[2]*pow(Xval[1], 3)) + 1)*((1.0L/4.0L)*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])))/(2*Time+1);
   values[10*(Time+1)+58*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(-0.0833333333333333*Xvalp1[2]*hstep*pow(Xvalp1[1], 3)*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]) + 0.0625*Xvalp1[2]*hstep*pow(Xvalp1[1], 3)*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]) + 0.333333333333333*hstep*(0.166666666666667*hstep*(-1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) - 1)*((1.0L/4.0L)*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.25*hstep*(0.125*hstep*((1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4) + (1.0L/4.0L)*Pval[4] + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) + 0.5)*((1.0L/4.0L)*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])))/(2*Time+1);
   values[11*(Time+1)+58*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0234375*Pval[0]*Xval[2]*pow(hstep, 2)*pow(Xval[1], 2)*(Pval[1] - Xval[0])*((1.0L/4.0L)*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.0416666666666667*Pval[0]*Xval[2]*pow(hstep, 2)*pow(Xval[1], 2)*(Pval[1] - Xval[0])*((1.0L/4.0L)*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.25*Xval[2]*hstep*pow(Xval[1], 2)*(Pval[1] - Xval[0])*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]) + 0.1875*Xval[2]*hstep*pow(Xval[1], 2)*(Pval[1] - Xval[0])*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
   values[11*(Time+1)+58*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(-0.0234375*Pval[0]*Xvalp1[2]*pow(hstep, 2)*pow(Xvalp1[1], 2)*(Pval[1] - Xvalp1[0])*((1.0L/4.0L)*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.0416666666666667*Pval[0]*Xvalp1[2]*pow(hstep, 2)*pow(Xvalp1[1], 2)*(Pval[1] - Xvalp1[0])*((1.0L/4.0L)*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.25*Xvalp1[2]*hstep*pow(Xvalp1[1], 2)*(Pval[1] - Xvalp1[0])*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]) - 0.1875*Xvalp1[2]*hstep*pow(Xvalp1[1], 2)*(Pval[1] - Xvalp1[0])*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
   values[12*(Time+1)+58*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0078125*Pval[0]*pow(hstep, 2)*pow(Xval[1], 3)*(Pval[1] - Xval[0])*((1.0L/4.0L)*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.0138888888888889*Pval[0]*pow(hstep, 2)*pow(Xval[1], 3)*(Pval[1] - Xval[0])*((1.0L/4.0L)*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.0833333333333333*hstep*pow(Xval[1], 3)*(Pval[1] - Xval[0])*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]) + 0.0625*hstep*pow(Xval[1], 3)*(Pval[1] - Xval[0])*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
   values[12*(Time+1)+58*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(-0.0078125*Pval[0]*pow(hstep, 2)*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])*((1.0L/4.0L)*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.0138888888888889*Pval[0]*pow(hstep, 2)*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])*((1.0L/4.0L)*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.0833333333333333*hstep*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]) - 0.0625*hstep*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
   values[13*(Time+1)+58*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.03125*Pval[2]*pow(hstep, 2)*pow(Xval[3], 3)*(Pval[3] - Xval[0])*((1.0L/4.0L)*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.0555555555555556*Pval[2]*pow(hstep, 2)*pow(Xval[3], 3)*(Pval[3] - Xval[0])*((1.0L/4.0L)*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])))/(2*Time+1);
   values[13*(Time+1)+58*(Time)+0+ 1+jt] += obj_factor*bounds[0][3]*(-0.03125*Pval[2]*pow(hstep, 2)*pow(Xvalp1[3], 3)*(Pval[3] - Xvalp1[0])*((1.0L/4.0L)*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.0555555555555556*Pval[2]*pow(hstep, 2)*pow(Xvalp1[3], 3)*(Pval[3] - Xvalp1[0])*((1.0L/4.0L)*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])))/(2*Time+1);
   values[14*(Time+1)+58*(Time)+0+jt] += obj_factor*bounds[0][3]*(-0.333333333333333*Xval2[2]*hstep*pow(Xval2[1], 3)*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]) + 0.0555555555555556*pow(hstep, 2)*(-Pval[2]*pow(Xval2[3], 4) - Pval[4] - Pval[0]*Xval2[2]*pow(Xval2[1], 3))*((1.0L/4.0L)*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) - 0.25*hstep*((1.0L/4.0L)*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])))/(2*Time+1);
   values[14*(Time+1)+59*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.166666666666667*Pval[0]*Xval2[2]*pow(hstep, 2)*pow(Xval2[1], 2)*(Pval[1] - Xval2[0])*((1.0L/4.0L)*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 1.0*Xval2[2]*hstep*pow(Xval2[1], 2)*(Pval[1] - Xval2[0])*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]))/(2*Time+1);
   values[14*(Time+1)+60*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.0555555555555556*Pval[0]*pow(hstep, 2)*pow(Xval2[1], 3)*(Pval[1] - Xval2[0])*((1.0L/4.0L)*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.333333333333333*hstep*pow(Xval2[1], 3)*(Pval[1] - Xval2[0])*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]))/(2*Time+1);
   values[14*(Time+1)+61*(Time)+0+jt] += obj_factor*bounds[0][3]*(0.222222222222222*Pval[2]*pow(hstep, 2)*pow(Xval2[3], 3)*(Pval[3] - Xval2[0])*((1.0L/4.0L)*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])))/(2*Time+1);
   values[14*(Time+1)+62*(Time)+1+jt] += obj_factor*bounds[0][3]*(0.25*hstep*(0.125*hstep*(-1.0L/4.0L*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xval[2]*pow(Xval[1], 3)) + 0.5)*((1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) + 0.333333333333333*hstep*(0.166666666666667*hstep*(-1.0L/4.0L*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xval[2]*pow(Xval[1], 3)) + 1)*((1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3) + Pval[0]*Xval2[2]*pow(Xval2[1], 3) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)))/(2*Time+1);
   values[14*(Time+1)+62*(Time)+1+ 1+jt] += obj_factor*bounds[0][3]*(0.333333333333333*hstep*(0.166666666666667*hstep*(-1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) - 1)*((1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3) + Pval[0]*Xval2[2]*pow(Xval2[1], 3) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) + 0.25*hstep*(0.125*hstep*((1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4) + (1.0L/4.0L)*Pval[4] + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) + 0.5)*((1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)))/(2*Time+1);
   values[15*(Time+1)+62*(Time)+1+jt] += obj_factor*bounds[0][3]*(0.0234375*Pval[0]*Xval[2]*pow(hstep, 2)*pow(Xval[1], 2)*(Pval[1] - Xval[0])*((1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) + 0.0416666666666667*Pval[0]*Xval[2]*pow(hstep, 2)*pow(Xval[1], 2)*(Pval[1] - Xval[0])*((1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3) + Pval[0]*Xval2[2]*pow(Xval2[1], 3) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) + 0.25*Pval[0]*Xval[2]*hstep*pow(Xval[1], 2)*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]) + 0.1875*Pval[0]*Xval[2]*hstep*pow(Xval[1], 2)*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
   values[15*(Time+1)+62*(Time)+1+ 1+jt] += obj_factor*bounds[0][3]*(-0.0234375*Pval[0]*Xvalp1[2]*pow(hstep, 2)*pow(Xvalp1[1], 2)*(Pval[1] - Xvalp1[0])*((1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) + 0.0416666666666667*Pval[0]*Xvalp1[2]*pow(hstep, 2)*pow(Xvalp1[1], 2)*(Pval[1] - Xvalp1[0])*((1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3) + Pval[0]*Xval2[2]*pow(Xval2[1], 3) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) + 0.25*Pval[0]*Xvalp1[2]*hstep*pow(Xvalp1[1], 2)*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]) - 0.1875*Pval[0]*Xvalp1[2]*hstep*pow(Xvalp1[1], 2)*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
   values[16*(Time+1)+62*(Time)+1+jt] += obj_factor*bounds[0][3]*(0.0078125*Pval[0]*pow(hstep, 2)*pow(Xval[1], 3)*(Pval[1] - Xval[0])*((1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) + 0.0138888888888889*Pval[0]*pow(hstep, 2)*pow(Xval[1], 3)*(Pval[1] - Xval[0])*((1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3) + Pval[0]*Xval2[2]*pow(Xval2[1], 3) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) + 0.0833333333333333*Pval[0]*hstep*pow(Xval[1], 3)*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]) + 0.0625*Pval[0]*hstep*pow(Xval[1], 3)*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
   values[16*(Time+1)+62*(Time)+1+ 1+jt] += obj_factor*bounds[0][3]*(-0.0078125*Pval[0]*pow(hstep, 2)*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])*((1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) + 0.0138888888888889*Pval[0]*pow(hstep, 2)*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])*((1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3) + Pval[0]*Xval2[2]*pow(Xval2[1], 3) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) + 0.0833333333333333*Pval[0]*hstep*pow(Xvalp1[1], 3)*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]) - 0.0625*Pval[0]*hstep*pow(Xvalp1[1], 3)*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
   values[17*(Time+1)+62*(Time)+1+jt] += obj_factor*bounds[0][3]*(0.03125*Pval[2]*pow(hstep, 2)*pow(Xval[3], 3)*(Pval[3] - Xval[0])*((1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) + 0.0555555555555556*Pval[2]*pow(hstep, 2)*pow(Xval[3], 3)*(Pval[3] - Xval[0])*((1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3) + Pval[0]*Xval2[2]*pow(Xval2[1], 3) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)))/(2*Time+1);
   values[17*(Time+1)+62*(Time)+1+ 1+jt] += obj_factor*bounds[0][3]*(-0.03125*Pval[2]*pow(hstep, 2)*pow(Xvalp1[3], 3)*(Pval[3] - Xvalp1[0])*((1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) + 0.0555555555555556*Pval[2]*pow(hstep, 2)*pow(Xvalp1[3], 3)*(Pval[3] - Xvalp1[0])*((1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3) + Pval[0]*Xval2[2]*pow(Xval2[1], 3) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)))/(2*Time+1);
   values[18*(Time+1)+62*(Time)+1+jt] += obj_factor*bounds[0][3]*(0.0555555555555556*pow(hstep, 2)*(-Pval[2]*pow(Xval2[3], 4) - Pval[4] - Pval[0]*Xval2[2]*pow(Xval2[1], 3))*((1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3) + Pval[0]*Xval2[2]*pow(Xval2[1], 3) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) - 0.25*hstep*((1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)))/(2*Time+1);
   values[18*(Time+1)+63*(Time)+1+jt] += obj_factor*bounds[0][3]*(0.166666666666667*Pval[0]*Xval2[2]*pow(hstep, 2)*pow(Xval2[1], 2)*(Pval[1] - Xval2[0])*((1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3) + Pval[0]*Xval2[2]*pow(Xval2[1], 3) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) + 1.0*Pval[0]*Xval2[2]*hstep*pow(Xval2[1], 2)*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]))/(2*Time+1);
   values[18*(Time+1)+64*(Time)+1+jt] += obj_factor*bounds[0][3]*(0.0555555555555556*Pval[0]*pow(hstep, 2)*pow(Xval2[1], 3)*(Pval[1] - Xval2[0])*((1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3) + Pval[0]*Xval2[2]*pow(Xval2[1], 3) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) + 0.333333333333333*Pval[0]*hstep*pow(Xval2[1], 3)*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]))/(2*Time+1);
   values[18*(Time+1)+65*(Time)+1+jt] += obj_factor*bounds[0][3]*(0.222222222222222*Pval[2]*pow(hstep, 2)*pow(Xval2[3], 3)*(Pval[3] - Xval2[0])*((1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3) + Pval[0]*Xval2[2]*pow(Xval2[1], 3) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)))/(2*Time+1);
   values[18*(Time+1)+66*(Time)+3+jt] += obj_factor*bounds[0][3]*(-0.0833333333333333*hstep*pow(Xval[3], 4)*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]) - 0.0625*hstep*pow(Xval[3], 4)*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]) + 0.25*hstep*(0.125*hstep*(-1.0L/4.0L*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xval[2]*pow(Xval[1], 3)) + 0.5)*((1.0L/4.0L)*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0])) + 0.333333333333333*hstep*(0.166666666666667*hstep*(-1.0L/4.0L*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xval[2]*pow(Xval[1], 3)) + 1)*((1.0L/4.0L)*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0])))/(2*Time+1);
   values[18*(Time+1)+66*(Time)+3+ 1+jt] += obj_factor*bounds[0][3]*(-0.0833333333333333*hstep*pow(Xvalp1[3], 4)*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]) + 0.0625*hstep*pow(Xvalp1[3], 4)*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]) + 0.333333333333333*hstep*(0.166666666666667*hstep*(-1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) - 1)*((1.0L/4.0L)*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0])) + 0.25*hstep*(0.125*hstep*((1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4) + (1.0L/4.0L)*Pval[4] + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) + 0.5)*((1.0L/4.0L)*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0])))/(2*Time+1);
   values[19*(Time+1)+66*(Time)+3+jt] += obj_factor*bounds[0][3]*(0.0234375*Pval[0]*Xval[2]*pow(hstep, 2)*pow(Xval[1], 2)*(Pval[1] - Xval[0])*((1.0L/4.0L)*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0])) + 0.0416666666666667*Pval[0]*Xval[2]*pow(hstep, 2)*pow(Xval[1], 2)*(Pval[1] - Xval[0])*((1.0L/4.0L)*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0])))/(2*Time+1);
   values[19*(Time+1)+66*(Time)+3+ 1+jt] += obj_factor*bounds[0][3]*(-0.0234375*Pval[0]*Xvalp1[2]*pow(hstep, 2)*pow(Xvalp1[1], 2)*(Pval[1] - Xvalp1[0])*((1.0L/4.0L)*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0])) + 0.0416666666666667*Pval[0]*Xvalp1[2]*pow(hstep, 2)*pow(Xvalp1[1], 2)*(Pval[1] - Xvalp1[0])*((1.0L/4.0L)*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0])))/(2*Time+1);
   values[20*(Time+1)+66*(Time)+3+jt] += obj_factor*bounds[0][3]*(0.0078125*Pval[0]*pow(hstep, 2)*pow(Xval[1], 3)*(Pval[1] - Xval[0])*((1.0L/4.0L)*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0])) + 0.0138888888888889*Pval[0]*pow(hstep, 2)*pow(Xval[1], 3)*(Pval[1] - Xval[0])*((1.0L/4.0L)*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0])))/(2*Time+1);
   values[20*(Time+1)+66*(Time)+3+ 1+jt] += obj_factor*bounds[0][3]*(-0.0078125*Pval[0]*pow(hstep, 2)*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])*((1.0L/4.0L)*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0])) + 0.0138888888888889*Pval[0]*pow(hstep, 2)*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])*((1.0L/4.0L)*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0])))/(2*Time+1);
   values[21*(Time+1)+66*(Time)+3+jt] += obj_factor*bounds[0][3]*(0.03125*Pval[2]*pow(hstep, 2)*pow(Xval[3], 3)*(Pval[3] - Xval[0])*((1.0L/4.0L)*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0])) + 0.0555555555555556*Pval[2]*pow(hstep, 2)*pow(Xval[3], 3)*(Pval[3] - Xval[0])*((1.0L/4.0L)*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0])) + 0.333333333333333*hstep*pow(Xval[3], 3)*(Pval[3] - Xval[0])*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]) + 0.25*hstep*pow(Xval[3], 3)*(Pval[3] - Xval[0])*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
   values[21*(Time+1)+66*(Time)+3+ 1+jt] += obj_factor*bounds[0][3]*(-0.03125*Pval[2]*pow(hstep, 2)*pow(Xvalp1[3], 3)*(Pval[3] - Xvalp1[0])*((1.0L/4.0L)*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0])) + 0.0555555555555556*Pval[2]*pow(hstep, 2)*pow(Xvalp1[3], 3)*(Pval[3] - Xvalp1[0])*((1.0L/4.0L)*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0])) + 0.333333333333333*hstep*pow(Xvalp1[3], 3)*(Pval[3] - Xvalp1[0])*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]) - 0.25*hstep*pow(Xvalp1[3], 3)*(Pval[3] - Xvalp1[0])*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
   values[22*(Time+1)+66*(Time)+3+jt] += obj_factor*bounds[0][3]*(0.0555555555555556*pow(hstep, 2)*(-Pval[2]*pow(Xval2[3], 4) - Pval[4] - Pval[0]*Xval2[2]*pow(Xval2[1], 3))*((1.0L/4.0L)*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0])) - 0.333333333333333*hstep*pow(Xval2[3], 4)*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]) - 0.25*hstep*((1.0L/4.0L)*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0])))/(2*Time+1);
   values[22*(Time+1)+67*(Time)+3+jt] += obj_factor*bounds[0][3]*(0.166666666666667*Pval[0]*Xval2[2]*pow(hstep, 2)*pow(Xval2[1], 2)*(Pval[1] - Xval2[0])*((1.0L/4.0L)*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0])))/(2*Time+1);
   values[22*(Time+1)+68*(Time)+3+jt] += obj_factor*bounds[0][3]*(0.0555555555555556*Pval[0]*pow(hstep, 2)*pow(Xval2[1], 3)*(Pval[1] - Xval2[0])*((1.0L/4.0L)*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0])))/(2*Time+1);
   values[22*(Time+1)+69*(Time)+3+jt] += obj_factor*bounds[0][3]*(0.222222222222222*Pval[2]*pow(hstep, 2)*pow(Xval2[3], 3)*(Pval[3] - Xval2[0])*((1.0L/4.0L)*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0])) + 1.33333333333333*hstep*pow(Xval2[3], 3)*(Pval[3] - Xval2[0])*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]))/(2*Time+1);
   values[22*(Time+1)+70*(Time)+6+jt] += obj_factor*bounds[0][3]*(0.25*hstep*((1.0L/4.0L)*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4))*(0.125*hstep*(-1.0L/4.0L*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xval[2]*pow(Xval[1], 3)) + 0.5) + 0.333333333333333*hstep*(0.166666666666667*hstep*(-1.0L/4.0L*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xval[2]*pow(Xval[1], 3)) + 1)*((1.0L/4.0L)*Pval[2]*pow(Xval[3], 4) + Pval[2]*pow(Xval2[3], 4) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)))/(2*Time+1);
   values[22*(Time+1)+70*(Time)+6+ 1+jt] += obj_factor*bounds[0][3]*(0.25*hstep*((1.0L/4.0L)*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4))*(0.125*hstep*((1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4) + (1.0L/4.0L)*Pval[4] + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) + 0.5) + 0.333333333333333*hstep*(0.166666666666667*hstep*(-1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) - 1)*((1.0L/4.0L)*Pval[2]*pow(Xval[3], 4) + Pval[2]*pow(Xval2[3], 4) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)))/(2*Time+1);
   values[23*(Time+1)+70*(Time)+6+jt] += obj_factor*bounds[0][3]*(0.0234375*Pval[0]*Xval[2]*pow(hstep, 2)*pow(Xval[1], 2)*(Pval[1] - Xval[0])*((1.0L/4.0L)*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)) + 0.0416666666666667*Pval[0]*Xval[2]*pow(hstep, 2)*pow(Xval[1], 2)*(Pval[1] - Xval[0])*((1.0L/4.0L)*Pval[2]*pow(Xval[3], 4) + Pval[2]*pow(Xval2[3], 4) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)))/(2*Time+1);
   values[23*(Time+1)+70*(Time)+6+ 1+jt] += obj_factor*bounds[0][3]*(-0.0234375*Pval[0]*Xvalp1[2]*pow(hstep, 2)*pow(Xvalp1[1], 2)*(Pval[1] - Xvalp1[0])*((1.0L/4.0L)*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)) + 0.0416666666666667*Pval[0]*Xvalp1[2]*pow(hstep, 2)*pow(Xvalp1[1], 2)*(Pval[1] - Xvalp1[0])*((1.0L/4.0L)*Pval[2]*pow(Xval[3], 4) + Pval[2]*pow(Xval2[3], 4) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)))/(2*Time+1);
   values[24*(Time+1)+70*(Time)+6+jt] += obj_factor*bounds[0][3]*(0.0078125*Pval[0]*pow(hstep, 2)*pow(Xval[1], 3)*(Pval[1] - Xval[0])*((1.0L/4.0L)*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)) + 0.0138888888888889*Pval[0]*pow(hstep, 2)*pow(Xval[1], 3)*(Pval[1] - Xval[0])*((1.0L/4.0L)*Pval[2]*pow(Xval[3], 4) + Pval[2]*pow(Xval2[3], 4) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)))/(2*Time+1);
   values[24*(Time+1)+70*(Time)+6+ 1+jt] += obj_factor*bounds[0][3]*(-0.0078125*Pval[0]*pow(hstep, 2)*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])*((1.0L/4.0L)*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)) + 0.0138888888888889*Pval[0]*pow(hstep, 2)*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])*((1.0L/4.0L)*Pval[2]*pow(Xval[3], 4) + Pval[2]*pow(Xval2[3], 4) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)))/(2*Time+1);
   values[25*(Time+1)+70*(Time)+6+jt] += obj_factor*bounds[0][3]*(0.03125*Pval[2]*pow(hstep, 2)*pow(Xval[3], 3)*(Pval[3] - Xval[0])*((1.0L/4.0L)*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)) + 0.0555555555555556*Pval[2]*pow(hstep, 2)*pow(Xval[3], 3)*(Pval[3] - Xval[0])*((1.0L/4.0L)*Pval[2]*pow(Xval[3], 4) + Pval[2]*pow(Xval2[3], 4) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)) + 0.333333333333333*Pval[2]*hstep*pow(Xval[3], 3)*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]) + 0.25*Pval[2]*hstep*pow(Xval[3], 3)*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
   values[25*(Time+1)+70*(Time)+6+ 1+jt] += obj_factor*bounds[0][3]*(-0.03125*Pval[2]*pow(hstep, 2)*pow(Xvalp1[3], 3)*(Pval[3] - Xvalp1[0])*((1.0L/4.0L)*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)) + 0.0555555555555556*Pval[2]*pow(hstep, 2)*pow(Xvalp1[3], 3)*(Pval[3] - Xvalp1[0])*((1.0L/4.0L)*Pval[2]*pow(Xval[3], 4) + Pval[2]*pow(Xval2[3], 4) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)) + 0.333333333333333*Pval[2]*hstep*pow(Xvalp1[3], 3)*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]) - 0.25*Pval[2]*hstep*pow(Xvalp1[3], 3)*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
   values[26*(Time+1)+70*(Time)+6+jt] += obj_factor*bounds[0][3]*(0.0555555555555556*pow(hstep, 2)*((1.0L/4.0L)*Pval[2]*pow(Xval[3], 4) + Pval[2]*pow(Xval2[3], 4) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4))*(-Pval[2]*pow(Xval2[3], 4) - Pval[4] - Pval[0]*Xval2[2]*pow(Xval2[1], 3)) - 0.25*hstep*((1.0L/4.0L)*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)))/(2*Time+1);
   values[26*(Time+1)+71*(Time)+6+jt] += obj_factor*bounds[0][3]*(0.166666666666667*Pval[0]*Xval2[2]*pow(hstep, 2)*pow(Xval2[1], 2)*(Pval[1] - Xval2[0])*((1.0L/4.0L)*Pval[2]*pow(Xval[3], 4) + Pval[2]*pow(Xval2[3], 4) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)))/(2*Time+1);
   values[26*(Time+1)+72*(Time)+6+jt] += obj_factor*bounds[0][3]*(0.0555555555555556*Pval[0]*pow(hstep, 2)*pow(Xval2[1], 3)*(Pval[1] - Xval2[0])*((1.0L/4.0L)*Pval[2]*pow(Xval[3], 4) + Pval[2]*pow(Xval2[3], 4) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)))/(2*Time+1);
   values[26*(Time+1)+73*(Time)+6+jt] += obj_factor*bounds[0][3]*(0.222222222222222*Pval[2]*pow(hstep, 2)*pow(Xval2[3], 3)*(Pval[3] - Xval2[0])*((1.0L/4.0L)*Pval[2]*pow(Xval[3], 4) + Pval[2]*pow(Xval2[3], 4) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)) + 1.33333333333333*Pval[2]*hstep*pow(Xval2[3], 3)*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]))/(2*Time+1);
   values[26*(Time+1)+74*(Time)+10+jt] += obj_factor*bounds[0][3]*(0.25*hstep*(-1.0L/4.0L*Xval[0] + (1.0L/4.0L)*Xvalp1[0])*(0.125*hstep*(-1.0L/4.0L*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xval[2]*pow(Xval[1], 3)) + 0.5) + 0.333333333333333*hstep*(0.166666666666667*hstep*(-1.0L/4.0L*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xval[2]*pow(Xval[1], 3)) + 1)*((3.0L/2.0L)*Pval[5] - 1.0L/4.0L*Xval[0] - Xval2[0] - 1.0L/4.0L*Xvalp1[0]) - 0.0833333333333333*hstep*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]) - 0.0625*hstep*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
   values[26*(Time+1)+74*(Time)+10+ 1+jt] += obj_factor*bounds[0][3]*(0.25*hstep*(-1.0L/4.0L*Xval[0] + (1.0L/4.0L)*Xvalp1[0])*(0.125*hstep*((1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4) + (1.0L/4.0L)*Pval[4] + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) + 0.5) + 0.333333333333333*hstep*(0.166666666666667*hstep*(-1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) - 1)*((3.0L/2.0L)*Pval[5] - 1.0L/4.0L*Xval[0] - Xval2[0] - 1.0L/4.0L*Xvalp1[0]) - 0.0833333333333333*hstep*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]) + 0.0625*hstep*(0.125*hstep*((1.0L/4.0L)*Ival[0] - 1.0L/4.0L*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) - 1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) - 1.0L/4.0L*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + 0.5*Xval[0] - Xval2[0] + 0.5*Xvalp1[0]))/(2*Time+1);
   values[27*(Time+1)+74*(Time)+10+jt] += obj_factor*bounds[0][3]*(0.0234375*Pval[0]*Xval[2]*pow(hstep, 2)*pow(Xval[1], 2)*(Pval[1] - Xval[0])*(-1.0L/4.0L*Xval[0] + (1.0L/4.0L)*Xvalp1[0]) + 0.0416666666666667*Pval[0]*Xval[2]*pow(hstep, 2)*pow(Xval[1], 2)*(Pval[1] - Xval[0])*((3.0L/2.0L)*Pval[5] - 1.0L/4.0L*Xval[0] - Xval2[0] - 1.0L/4.0L*Xvalp1[0]))/(2*Time+1);
   values[27*(Time+1)+74*(Time)+10+ 1+jt] += obj_factor*bounds[0][3]*(-0.0234375*Pval[0]*Xvalp1[2]*pow(hstep, 2)*pow(Xvalp1[1], 2)*(Pval[1] - Xvalp1[0])*(-1.0L/4.0L*Xval[0] + (1.0L/4.0L)*Xvalp1[0]) + 0.0416666666666667*Pval[0]*Xvalp1[2]*pow(hstep, 2)*pow(Xvalp1[1], 2)*(Pval[1] - Xvalp1[0])*((3.0L/2.0L)*Pval[5] - 1.0L/4.0L*Xval[0] - Xval2[0] - 1.0L/4.0L*Xvalp1[0]))/(2*Time+1);
   values[28*(Time+1)+74*(Time)+10+jt] += obj_factor*bounds[0][3]*(0.0078125*Pval[0]*pow(hstep, 2)*pow(Xval[1], 3)*(Pval[1] - Xval[0])*(-1.0L/4.0L*Xval[0] + (1.0L/4.0L)*Xvalp1[0]) + 0.0138888888888889*Pval[0]*pow(hstep, 2)*pow(Xval[1], 3)*(Pval[1] - Xval[0])*((3.0L/2.0L)*Pval[5] - 1.0L/4.0L*Xval[0] - Xval2[0] - 1.0L/4.0L*Xvalp1[0]))/(2*Time+1);
   values[28*(Time+1)+74*(Time)+10+ 1+jt] += obj_factor*bounds[0][3]*(-0.0078125*Pval[0]*pow(hstep, 2)*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])*(-1.0L/4.0L*Xval[0] + (1.0L/4.0L)*Xvalp1[0]) + 0.0138888888888889*Pval[0]*pow(hstep, 2)*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])*((3.0L/2.0L)*Pval[5] - 1.0L/4.0L*Xval[0] - Xval2[0] - 1.0L/4.0L*Xvalp1[0]))/(2*Time+1);
   values[29*(Time+1)+74*(Time)+10+jt] += obj_factor*bounds[0][3]*(0.03125*Pval[2]*pow(hstep, 2)*pow(Xval[3], 3)*(Pval[3] - Xval[0])*(-1.0L/4.0L*Xval[0] + (1.0L/4.0L)*Xvalp1[0]) + 0.0555555555555556*Pval[2]*pow(hstep, 2)*pow(Xval[3], 3)*(Pval[3] - Xval[0])*((3.0L/2.0L)*Pval[5] - 1.0L/4.0L*Xval[0] - Xval2[0] - 1.0L/4.0L*Xvalp1[0]))/(2*Time+1);
   values[29*(Time+1)+74*(Time)+10+ 1+jt] += obj_factor*bounds[0][3]*(-0.03125*Pval[2]*pow(hstep, 2)*pow(Xvalp1[3], 3)*(Pval[3] - Xvalp1[0])*(-1.0L/4.0L*Xval[0] + (1.0L/4.0L)*Xvalp1[0]) + 0.0555555555555556*Pval[2]*pow(hstep, 2)*pow(Xvalp1[3], 3)*(Pval[3] - Xvalp1[0])*((3.0L/2.0L)*Pval[5] - 1.0L/4.0L*Xval[0] - Xval2[0] - 1.0L/4.0L*Xvalp1[0]))/(2*Time+1);
   values[30*(Time+1)+74*(Time)+10+jt] += obj_factor*bounds[0][3]*(0.0555555555555556*pow(hstep, 2)*(-Pval[2]*pow(Xval2[3], 4) - Pval[4] - Pval[0]*Xval2[2]*pow(Xval2[1], 3))*((3.0L/2.0L)*Pval[5] - 1.0L/4.0L*Xval[0] - Xval2[0] - 1.0L/4.0L*Xvalp1[0]) - 0.25*hstep*(-1.0L/4.0L*Xval[0] + (1.0L/4.0L)*Xvalp1[0]) - 0.333333333333333*hstep*(0.166666666666667*hstep*((1.0L/4.0L)*Ival[0] + Ival2[0] + (1.0L/4.0L)*Ivalp1[0] + (1.0L/4.0L)*Pval[2]*pow(Xval[3], 4)*(Pval[3] - Xval[0]) + Pval[2]*pow(Xval2[3], 4)*(Pval[3] - Xval2[0]) + (1.0L/4.0L)*Pval[2]*pow(Xvalp1[3], 4)*(Pval[3] - Xvalp1[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xval[0]) + Pval[4]*(Pval[5] - Xval2[0]) + (1.0L/4.0L)*Pval[4]*(Pval[5] - Xvalp1[0]) + (1.0L/4.0L)*Pval[0]*Xval[2]*pow(Xval[1], 3)*(Pval[1] - Xval[0]) + Pval[0]*Xval2[2]*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]) + (1.0L/4.0L)*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0])) + Xval[0] - Xvalp1[0]))/(2*Time+1);
   values[30*(Time+1)+75*(Time)+10+jt] += obj_factor*bounds[0][3]*(0.166666666666667*Pval[0]*Xval2[2]*pow(hstep, 2)*pow(Xval2[1], 2)*(Pval[1] - Xval2[0])*((3.0L/2.0L)*Pval[5] - 1.0L/4.0L*Xval[0] - Xval2[0] - 1.0L/4.0L*Xvalp1[0]))/(2*Time+1);
   values[30*(Time+1)+76*(Time)+10+jt] += obj_factor*bounds[0][3]*(0.0555555555555556*Pval[0]*pow(hstep, 2)*pow(Xval2[1], 3)*(Pval[1] - Xval2[0])*((3.0L/2.0L)*Pval[5] - 1.0L/4.0L*Xval[0] - Xval2[0] - 1.0L/4.0L*Xvalp1[0]))/(2*Time+1);
   values[30*(Time+1)+77*(Time)+10+jt] += obj_factor*bounds[0][3]*(0.222222222222222*Pval[2]*pow(hstep, 2)*pow(Xval2[3], 3)*(Pval[3] - Xval2[0])*((3.0L/2.0L)*Pval[5] - 1.0L/4.0L*Xval[0] - Xval2[0] - 1.0L/4.0L*Xvalp1[0]))/(2*Time+1);
   values[30*(Time+1)+78*(Time)+15+jt] += obj_factor*bounds[0][3]*(0.5*Pval[4]*hstep*(0.166666666666667*hstep*(-1.0L/4.0L*Pval[2]*pow(Xval[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xval[2]*pow(Xval[1], 3)) + 1))/(2*Time+1);
   values[30*(Time+1)+78*(Time)+15+ 1+jt] += obj_factor*bounds[0][3]*(0.5*Pval[4]*hstep*(0.166666666666667*hstep*(-1.0L/4.0L*Pval[2]*pow(Xvalp1[3], 4) - 1.0L/4.0L*Pval[4] - 1.0L/4.0L*Pval[0]*Xvalp1[2]*pow(Xvalp1[1], 3)) - 1))/(2*Time+1);
   values[31*(Time+1)+78*(Time)+15+jt] += obj_factor*bounds[0][3]*(0.0625*Pval[4]*Pval[0]*Xval[2]*pow(hstep, 2)*pow(Xval[1], 2)*(Pval[1] - Xval[0]))/(2*Time+1);
   values[31*(Time+1)+78*(Time)+15+ 1+jt] += obj_factor*bounds[0][3]*(0.0625*Pval[4]*Pval[0]*Xvalp1[2]*pow(hstep, 2)*pow(Xvalp1[1], 2)*(Pval[1] - Xvalp1[0]))/(2*Time+1);
   values[32*(Time+1)+78*(Time)+15+jt] += obj_factor*bounds[0][3]*(0.0208333333333333*Pval[4]*Pval[0]*pow(hstep, 2)*pow(Xval[1], 3)*(Pval[1] - Xval[0]))/(2*Time+1);
   values[32*(Time+1)+78*(Time)+15+ 1+jt] += obj_factor*bounds[0][3]*(0.0208333333333333*Pval[4]*Pval[0]*pow(hstep, 2)*pow(Xvalp1[1], 3)*(Pval[1] - Xvalp1[0]))/(2*Time+1);
   values[33*(Time+1)+78*(Time)+15+jt] += obj_factor*bounds[0][3]*(0.0833333333333333*Pval[2]*Pval[4]*pow(hstep, 2)*pow(Xval[3], 3)*(Pval[3] - Xval[0]))/(2*Time+1);
   values[33*(Time+1)+78*(Time)+15+ 1+jt] += obj_factor*bounds[0][3]*(0.0833333333333333*Pval[2]*Pval[4]*pow(hstep, 2)*pow(Xvalp1[3], 3)*(Pval[3] - Xvalp1[0]))/(2*Time+1);
   values[34*(Time+1)+78*(Time)+15+jt] += obj_factor*bounds[0][3]*(0.0833333333333333*Pval[4]*pow(hstep, 2)*(-Pval[2]*pow(Xval2[3], 4) - Pval[4] - Pval[0]*Xval2[2]*pow(Xval2[1], 3)))/(2*Time+1);
   values[34*(Time+1)+79*(Time)+15+jt] += obj_factor*bounds[0][3]*(0.25*Pval[4]*Pval[0]*Xval2[2]*pow(hstep, 2)*pow(Xval2[1], 2)*(Pval[1] - Xval2[0]))/(2*Time+1);
   values[34*(Time+1)+80*(Time)+15+jt] += obj_factor*bounds[0][3]*(0.0833333333333333*Pval[4]*Pval[0]*pow(hstep, 2)*pow(Xval2[1], 3)*(Pval[1] - Xval2[0]))/(2*Time+1);
   values[34*(Time+1)+81*(Time)+15+jt] += obj_factor*bounds[0][3]*(0.333333333333333*Pval[2]*Pval[4]*pow(hstep, 2)*pow(Xval2[3], 3)*(Pval[3] - Xval2[0]))/(2*Time+1);
   values[0*(Time+1)+0*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
// Betaspline
// Derivative
// Derivative
// Derivative
// Derivative
0.0868055555555556*pow(hstep, 2)*pow(-Xval[1]*Betasplinejac(Xval[0],1,1) + (-Xval[1] + 1)*Alphasplinejac(Xval[0],1,1), 2) + 0.333333333333333*hstep*(-Xval[1]*Betasplinehes(Xval[0],1,1,1) + (-Xval[1] + 1)*Alphasplinehes(Xval[0],1,1,1))*(0.166666666666667*hstep*(-Xval[1]*Betaspline(Xval[0], 1) - 4*Xval2[1]*Betaspline(Xval2[0], 1) - Xvalp1[1]*Betaspline(Xvalp1[0], 1) + (-Xval[1] + 1)*Alphaspline(Xval[0], 1) + (-4*Xval2[1] + 4)*Alphaspline(Xval2[0], 1) + (-Xvalp1[1] + 1)*Alphaspline(Xvalp1[0], 1)) + Xval[1] - Xvalp1[1]) + 0.25*hstep*(-Xval[1]*Betasplinehes(Xval[0],1,1,1) + (-Xval[1] + 1)*Alphasplinehes(Xval[0],1,1,1))*(0.125*hstep*(-Xval[1]*Betaspline(Xval[0], 1) + Xvalp1[1]*Betaspline(Xvalp1[0], 1) + (-Xval[1] + 1)*Alphaspline(Xval[0], 1) - (-Xvalp1[1] + 1)*Alphaspline(Xvalp1[0], 1)) + 0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1]))/(2*Time+1);
   values[1*(Time+1)+0*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
0.0555555555555556*pow(hstep, 2)*(-Xval[1]*Betasplinejac(Xval[0],1,1) + (-Xval[1] + 1)*Alphasplinejac(Xval[0],1,1))*(-Xvalp1[1]*Betasplinejac(Xvalp1[0],1,1) + (-Xvalp1[1] + 1)*Alphasplinejac(Xvalp1[0],1,1)) + 0.03125*pow(hstep, 2)*(-Xval[1]*Betasplinejac(Xval[0],1,1) + (-Xval[1] + 1)*Alphasplinejac(Xval[0],1,1))*(Xvalp1[1]*Betasplinejac(Xvalp1[0],1,1) - (-Xvalp1[1] + 1)*Alphasplinejac(Xvalp1[0],1,1)))/(2*Time+1);
   values[0*(Time+1)+0*(Time)+0+ 1+jt] += obj_factor*bounds[1][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
// Betaspline
// Derivative
// Derivative
// Derivative
// Derivative
0.0555555555555556*pow(hstep, 2)*pow(-Xvalp1[1]*Betasplinejac(Xvalp1[0],1,1) + (-Xvalp1[1] + 1)*Alphasplinejac(Xvalp1[0],1,1), 2) + 0.03125*pow(hstep, 2)*pow(Xvalp1[1]*Betasplinejac(Xvalp1[0],1,1) - (-Xvalp1[1] + 1)*Alphasplinejac(Xvalp1[0],1,1), 2) + 0.333333333333333*hstep*(-Xvalp1[1]*Betasplinehes(Xvalp1[0],1,1,1) + (-Xvalp1[1] + 1)*Alphasplinehes(Xvalp1[0],1,1,1))*(0.166666666666667*hstep*(-Xval[1]*Betaspline(Xval[0], 1) - 4*Xval2[1]*Betaspline(Xval2[0], 1) - Xvalp1[1]*Betaspline(Xvalp1[0], 1) + (-Xval[1] + 1)*Alphaspline(Xval[0], 1) + (-4*Xval2[1] + 4)*Alphaspline(Xval2[0], 1) + (-Xvalp1[1] + 1)*Alphaspline(Xvalp1[0], 1)) + Xval[1] - Xvalp1[1]) + 0.25*hstep*(Xvalp1[1]*Betasplinehes(Xvalp1[0],1,1,1) - (-Xvalp1[1] + 1)*Alphasplinehes(Xvalp1[0],1,1,1))*(0.125*hstep*(-Xval[1]*Betaspline(Xval[0], 1) + Xvalp1[1]*Betaspline(Xvalp1[0], 1) + (-Xval[1] + 1)*Alphaspline(Xval[0], 1) - (-Xvalp1[1] + 1)*Alphaspline(Xvalp1[0], 1)) + 0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1]))/(2*Time+1);
   values[1*(Time+1)+1*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
// Betaspline
// Derivative
// Derivative
0.125*hstep*(0.25*hstep*(-Alphaspline(Xval[0], 1) - Betaspline(Xval[0], 1)) + 1.0)*(-Xval[1]*Betasplinejac(Xval[0],1,1) + (-Xval[1] + 1)*Alphasplinejac(Xval[0],1,1)) + 0.166666666666667*hstep*(0.333333333333333*hstep*(-Alphaspline(Xval[0], 1) - Betaspline(Xval[0], 1)) + 2)*(-Xval[1]*Betasplinejac(Xval[0],1,1) + (-Xval[1] + 1)*Alphasplinejac(Xval[0],1,1)) + 0.333333333333333*hstep*(-Alphasplinejac(Xval[0],1,1) - Betasplinejac(Xval[0],1,1))*(0.166666666666667*hstep*(-Xval[1]*Betaspline(Xval[0], 1) - 4*Xval2[1]*Betaspline(Xval2[0], 1) - Xvalp1[1]*Betaspline(Xvalp1[0], 1) + (-Xval[1] + 1)*Alphaspline(Xval[0], 1) + (-4*Xval2[1] + 4)*Alphaspline(Xval2[0], 1) + (-Xvalp1[1] + 1)*Alphaspline(Xvalp1[0], 1)) + Xval[1] - Xvalp1[1]) + 0.25*hstep*(-Alphasplinejac(Xval[0],1,1) - Betasplinejac(Xval[0],1,1))*(0.125*hstep*(-Xval[1]*Betaspline(Xval[0], 1) + Xvalp1[1]*Betaspline(Xvalp1[0], 1) + (-Xval[1] + 1)*Alphaspline(Xval[0], 1) - (-Xvalp1[1] + 1)*Alphaspline(Xvalp1[0], 1)) + 0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1]))/(2*Time+1);
   values[2*(Time+1)+1*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not supported in C:
// Alphaspline
// Betaspline
// Derivative
// Derivative
0.125*hstep*(0.25*hstep*(-Alphaspline(Xval[0], 1) - Betaspline(Xval[0], 1)) + 1.0)*(Xvalp1[1]*Betasplinejac(Xvalp1[0],1,1) - (-Xvalp1[1] + 1)*Alphasplinejac(Xvalp1[0],1,1)) + 0.166666666666667*hstep*(0.333333333333333*hstep*(-Alphaspline(Xval[0], 1) - Betaspline(Xval[0], 1)) + 2)*(-Xvalp1[1]*Betasplinejac(Xvalp1[0],1,1) + (-Xvalp1[1] + 1)*Alphasplinejac(Xvalp1[0],1,1)))/(2*Time+1);
   values[2*(Time+1)+2*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not supported in C:
// Alphaspline
// Betaspline
(0.125*hstep*(-Alphaspline(Xval[0], 1) - Betaspline(Xval[0], 1)) + 0.5)*(0.25*hstep*(-Alphaspline(Xval[0], 1) - Betaspline(Xval[0], 1)) + 1.0) + (0.166666666666667*hstep*(-Alphaspline(Xval[0], 1) - Betaspline(Xval[0], 1)) + 1)*(0.333333333333333*hstep*(-Alphaspline(Xval[0], 1) - Betaspline(Xval[0], 1)) + 2))/(2*Time+1);
   values[3*(Time+1)+2*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not supported in C:
// Alphaspline
// Betaspline
// Derivative
// Derivative
0.166666666666667*hstep*(0.333333333333333*hstep*(-Alphaspline(Xvalp1[0], 1) - Betaspline(Xvalp1[0], 1)) - 2)*(-Xval[1]*Betasplinejac(Xval[0],1,1) + (-Xval[1] + 1)*Alphasplinejac(Xval[0],1,1)) + 0.125*hstep*(0.25*hstep*(Alphaspline(Xvalp1[0], 1) + Betaspline(Xvalp1[0], 1)) + 1.0)*(-Xval[1]*Betasplinejac(Xval[0],1,1) + (-Xval[1] + 1)*Alphasplinejac(Xval[0],1,1)))/(2*Time+1);
   values[1*(Time+1)+1*(Time)+0+ 1+jt] += obj_factor*bounds[1][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
// Betaspline
// Derivative
// Derivative
0.166666666666667*hstep*(0.333333333333333*hstep*(-Alphaspline(Xvalp1[0], 1) - Betaspline(Xvalp1[0], 1)) - 2)*(-Xvalp1[1]*Betasplinejac(Xvalp1[0],1,1) + (-Xvalp1[1] + 1)*Alphasplinejac(Xvalp1[0],1,1)) + 0.125*hstep*(0.25*hstep*(Alphaspline(Xvalp1[0], 1) + Betaspline(Xvalp1[0], 1)) + 1.0)*(Xvalp1[1]*Betasplinejac(Xvalp1[0],1,1) - (-Xvalp1[1] + 1)*Alphasplinejac(Xvalp1[0],1,1)) + 0.333333333333333*hstep*(-Alphasplinejac(Xvalp1[0],1,1) - Betasplinejac(Xvalp1[0],1,1))*(0.166666666666667*hstep*(-Xval[1]*Betaspline(Xval[0], 1) - 4*Xval2[1]*Betaspline(Xval2[0], 1) - Xvalp1[1]*Betaspline(Xvalp1[0], 1) + (-Xval[1] + 1)*Alphaspline(Xval[0], 1) + (-4*Xval2[1] + 4)*Alphaspline(Xval2[0], 1) + (-Xvalp1[1] + 1)*Alphaspline(Xvalp1[0], 1)) + Xval[1] - Xvalp1[1]) + 0.25*hstep*(Alphasplinejac(Xvalp1[0],1,1) + Betasplinejac(Xvalp1[0],1,1))*(0.125*hstep*(-Xval[1]*Betaspline(Xval[0], 1) + Xvalp1[1]*Betaspline(Xvalp1[0], 1) + (-Xval[1] + 1)*Alphaspline(Xval[0], 1) - (-Xvalp1[1] + 1)*Alphaspline(Xvalp1[0], 1)) + 0.5*Xval[1] - Xval2[1] + 0.5*Xvalp1[1]))/(2*Time+1);
   values[3*(Time+1)+3*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
(0.125*hstep*(-Alphaspline(Xval[0], 1) - Betaspline(Xval[0], 1)) + 0.5)*(0.25*hstep*(Alphaspline(Xvalp1[0], 1) + Betaspline(Xvalp1[0], 1)) + 1.0) + (0.166666666666667*hstep*(-Alphaspline(Xval[0], 1) - Betaspline(Xval[0], 1)) + 1)*(0.333333333333333*hstep*(-Alphaspline(Xvalp1[0], 1) - Betaspline(Xvalp1[0], 1)) - 2))/(2*Time+1);
   values[2*(Time+1)+2*(Time)+0+ 1+jt] += obj_factor*bounds[1][3]*(// Not supported in C:
// Alphaspline
// Betaspline
(0.166666666666667*hstep*(-Alphaspline(Xvalp1[0], 1) - Betaspline(Xvalp1[0], 1)) - 1)*(0.333333333333333*hstep*(-Alphaspline(Xvalp1[0], 1) - Betaspline(Xvalp1[0], 1)) - 2) + (0.125*hstep*(Alphaspline(Xvalp1[0], 1) + Betaspline(Xvalp1[0], 1)) + 0.5)*(0.25*hstep*(Alphaspline(Xvalp1[0], 1) + Betaspline(Xvalp1[0], 1)) + 1.0))/(2*Time+1);
   values[10*(Time+1)+16*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
0.0555555555555556*pow(hstep, 2)*(-Xval[1]*Betasplinejac(Xval[0],1,1) + (-Xval[1] + 1)*Alphasplinejac(Xval[0],1,1))*(-4*Xval2[1]*Betasplinejac(Xval2[0],1,1) + (-4*Xval2[1] + 4)*Alphasplinejac(Xval2[0],1,1)))/(2*Time+1);
   values[10*(Time+1)+17*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
0.0555555555555556*pow(hstep, 2)*(-4*Xval2[1]*Betasplinejac(Xval2[0],1,1) + (-4*Xval2[1] + 4)*Alphasplinejac(Xval2[0],1,1))*(-Xvalp1[1]*Betasplinejac(Xvalp1[0],1,1) + (-Xvalp1[1] + 1)*Alphasplinejac(Xvalp1[0],1,1)))/(2*Time+1);
   values[10*(Time+1)+18*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not supported in C:
// Alphaspline
// Betaspline
// Derivative
// Derivative
0.333333333333333*hstep*(0.166666666666667*hstep*(-Alphaspline(Xval[0], 1) - Betaspline(Xval[0], 1)) + 1)*(-4*Xval2[1]*Betasplinejac(Xval2[0],1,1) + (-4*Xval2[1] + 4)*Alphasplinejac(Xval2[0],1,1)))/(2*Time+1);
   values[10*(Time+1)+19*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not supported in C:
// Alphaspline
// Betaspline
// Derivative
// Derivative
0.333333333333333*hstep*(0.166666666666667*hstep*(-Alphaspline(Xvalp1[0], 1) - Betaspline(Xvalp1[0], 1)) - 1)*(-4*Xval2[1]*Betasplinejac(Xval2[0],1,1) + (-4*Xval2[1] + 4)*Alphasplinejac(Xval2[0],1,1)))/(2*Time+1);
   values[10*(Time+1)+24*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
// Betaspline
// Derivative
// Derivative
// Derivative
// Derivative
0.0555555555555556*pow(hstep, 2)*pow(-4*Xval2[1]*Betasplinejac(Xval2[0],1,1) + (-4*Xval2[1] + 4)*Alphasplinejac(Xval2[0],1,1), 2) + 0.333333333333333*hstep*(-4*Xval2[1]*Betasplinehes(Xval2[0],1,1,1) + (-4*Xval2[1] + 4)*Alphasplinehes(Xval2[0],1,1,1))*(0.166666666666667*hstep*(-Xval[1]*Betaspline(Xval[0], 1) - 4*Xval2[1]*Betaspline(Xval2[0], 1) - Xvalp1[1]*Betaspline(Xvalp1[0], 1) + (-Xval[1] + 1)*Alphaspline(Xval[0], 1) + (-4*Xval2[1] + 4)*Alphaspline(Xval2[0], 1) + (-Xvalp1[1] + 1)*Alphaspline(Xvalp1[0], 1)) + Xval[1] - Xvalp1[1]))/(2*Time+1);
   values[10*(Time+1)+25*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not supported in C:
// Alphaspline
// Betaspline
// Derivative
// Derivative
0.0555555555555556*pow(hstep, 2)*(-Xval[1]*Betasplinejac(Xval[0],1,1) + (-Xval[1] + 1)*Alphasplinejac(Xval[0],1,1))*(-4*Alphaspline(Xval2[0], 1) - 4*Betaspline(Xval2[0], 1)) - 0.25*hstep*(-Xval[1]*Betasplinejac(Xval[0],1,1) + (-Xval[1] + 1)*Alphasplinejac(Xval[0],1,1)))/(2*Time+1);
   values[10*(Time+1)+26*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not supported in C:
// Alphaspline
// Betaspline
// Derivative
// Derivative
0.0555555555555556*pow(hstep, 2)*(-Xvalp1[1]*Betasplinejac(Xvalp1[0],1,1) + (-Xvalp1[1] + 1)*Alphasplinejac(Xvalp1[0],1,1))*(-4*Alphaspline(Xval2[0], 1) - 4*Betaspline(Xval2[0], 1)) - 0.25*hstep*(Xvalp1[1]*Betasplinejac(Xvalp1[0],1,1) - (-Xvalp1[1] + 1)*Alphasplinejac(Xvalp1[0],1,1)))/(2*Time+1);
   values[10*(Time+1)+27*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
0.333333333333333*hstep*(0.166666666666667*hstep*(-Alphaspline(Xval[0], 1) - Betaspline(Xval[0], 1)) + 1)*(-4*Alphaspline(Xval2[0], 1) - 4*Betaspline(Xval2[0], 1)) - 0.25*hstep*(-Alphaspline(Xval[0], 1) - Betaspline(Xval[0], 1)) - 1.0)/(2*Time+1);
   values[10*(Time+1)+28*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
0.333333333333333*hstep*(0.166666666666667*hstep*(-Alphaspline(Xvalp1[0], 1) - Betaspline(Xvalp1[0], 1)) - 1)*(-4*Alphaspline(Xval2[0], 1) - 4*Betaspline(Xval2[0], 1)) - 0.25*hstep*(Alphaspline(Xvalp1[0], 1) + Betaspline(Xvalp1[0], 1)) - 1.0)/(2*Time+1);
   values[10*(Time+1)+33*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
// Betaspline
// Derivative
// Derivative
0.0555555555555556*pow(hstep, 2)*(-4*Xval2[1]*Betasplinejac(Xval2[0],1,1) + (-4*Xval2[1] + 4)*Alphasplinejac(Xval2[0],1,1))*(-4*Alphaspline(Xval2[0], 1) - 4*Betaspline(Xval2[0], 1)) + 0.333333333333333*hstep*(-4*Alphasplinejac(Xval2[0],1,1) - 4*Betasplinejac(Xval2[0],1,1))*(0.166666666666667*hstep*(-Xval[1]*Betaspline(Xval[0], 1) - 4*Xval2[1]*Betaspline(Xval2[0], 1) - Xvalp1[1]*Betaspline(Xvalp1[0], 1) + (-Xval[1] + 1)*Alphaspline(Xval[0], 1) + (-4*Xval2[1] + 4)*Alphaspline(Xval2[0], 1) + (-Xvalp1[1] + 1)*Alphaspline(Xvalp1[0], 1)) + Xval[1] - Xvalp1[1]))/(2*Time+1);
   values[10*(Time+1)+34*(Time)+0+jt] += obj_factor*bounds[1][3]*(// Not supported in C:
// Alphaspline
// Betaspline
0.0555555555555556*pow(hstep, 2)*pow(-4*Alphaspline(Xval2[0], 1) - 4*Betaspline(Xval2[0], 1), 2) + 2)/(2*Time+1);
   values[0*(Time+1)+0*(Time)+0+jt] += obj_factor*bounds[2][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
// Betaspline
// Derivative
// Derivative
// Derivative
// Derivative
0.0868055555555556*pow(hstep, 2)*pow(-Xval[2]*Betasplinejac(Xval[0],2,1) + (-Xval[2] + 1)*Alphasplinejac(Xval[0],2,1), 2) + 0.333333333333333*hstep*(-Xval[2]*Betasplinehes(Xval[0],2,1,1) + (-Xval[2] + 1)*Alphasplinehes(Xval[0],2,1,1))*(Xval[2] - Xvalp1[2] + 0.166666666666667*hstep*(-Xval[2]*Betaspline(Xval[0], 2) - 4*Xval2[2]*Betaspline(Xval2[0], 2) - Xvalp1[2]*Betaspline(Xvalp1[0], 2) + (-Xval[2] + 1)*Alphaspline(Xval[0], 2) + (-4*Xval2[2] + 4)*Alphaspline(Xval2[0], 2) + (-Xvalp1[2] + 1)*Alphaspline(Xvalp1[0], 2))) + 0.25*hstep*(-Xval[2]*Betasplinehes(Xval[0],2,1,1) + (-Xval[2] + 1)*Alphasplinehes(Xval[0],2,1,1))*(0.5*Xval[2] - Xval2[2] + 0.5*Xvalp1[2] + 0.125*hstep*(-Xval[2]*Betaspline(Xval[0], 2) + Xvalp1[2]*Betaspline(Xvalp1[0], 2) + (-Xval[2] + 1)*Alphaspline(Xval[0], 2) - (-Xvalp1[2] + 1)*Alphaspline(Xvalp1[0], 2))))/(2*Time+1);
   values[1*(Time+1)+0*(Time)+0+jt] += obj_factor*bounds[2][3]*(// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
0.0555555555555556*pow(hstep, 2)*(-Xval[2]*Betasplinejac(Xval[0],2,1) + (-Xval[2] + 1)*Alphasplinejac(Xval[0],2,1))*(-Xvalp1[2]*Betasplinejac(Xvalp1[0],2,1) + (-Xvalp1[2] + 1)*Alphasplinejac(Xvalp1[0],2,1)) + 0.03125*pow(hstep, 2)*(-Xval[2]*Betasplinejac(Xval[0],2,1) + (-Xval[2] + 1)*Alphasplinejac(Xval[0],2,1))*(Xvalp1[2]*Betasplinejac(Xvalp1[0],2,1) - (-Xvalp1[2] + 1)*Alphasplinejac(Xvalp1[0],2,1)))/(2*Time+1);
   values[0*(Time+1)+0*(Time)+0+ 1+jt] += obj_factor*bounds[2][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
// Betaspline
// Derivative
// Derivative
// Derivative
// Derivative
0.0555555555555556*pow(hstep, 2)*pow(-Xvalp1[2]*Betasplinejac(Xvalp1[0],2,1) + (-Xvalp1[2] + 1)*Alphasplinejac(Xvalp1[0],2,1), 2) + 0.03125*pow(hstep, 2)*pow(Xvalp1[2]*Betasplinejac(Xvalp1[0],2,1) - (-Xvalp1[2] + 1)*Alphasplinejac(Xvalp1[0],2,1), 2) + 0.333333333333333*hstep*(-Xvalp1[2]*Betasplinehes(Xvalp1[0],2,1,1) + (-Xvalp1[2] + 1)*Alphasplinehes(Xvalp1[0],2,1,1))*(Xval[2] - Xvalp1[2] + 0.166666666666667*hstep*(-Xval[2]*Betaspline(Xval[0], 2) - 4*Xval2[2]*Betaspline(Xval2[0], 2) - Xvalp1[2]*Betaspline(Xvalp1[0], 2) + (-Xval[2] + 1)*Alphaspline(Xval[0], 2) + (-4*Xval2[2] + 4)*Alphaspline(Xval2[0], 2) + (-Xvalp1[2] + 1)*Alphaspline(Xvalp1[0], 2))) + 0.25*hstep*(Xvalp1[2]*Betasplinehes(Xvalp1[0],2,1,1) - (-Xvalp1[2] + 1)*Alphasplinehes(Xvalp1[0],2,1,1))*(0.5*Xval[2] - Xval2[2] + 0.5*Xvalp1[2] + 0.125*hstep*(-Xval[2]*Betaspline(Xval[0], 2) + Xvalp1[2]*Betaspline(Xvalp1[0], 2) + (-Xval[2] + 1)*Alphaspline(Xval[0], 2) - (-Xvalp1[2] + 1)*Alphaspline(Xvalp1[0], 2))))/(2*Time+1);
   values[3*(Time+1)+4*(Time)+0+jt] += obj_factor*bounds[2][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
// Betaspline
// Derivative
// Derivative
0.125*hstep*(-Xval[2]*Betasplinejac(Xval[0],2,1) + (-Xval[2] + 1)*Alphasplinejac(Xval[0],2,1))*(0.25*hstep*(-Alphaspline(Xval[0], 2) - Betaspline(Xval[0], 2)) + 1.0) + 0.166666666666667*hstep*(-Xval[2]*Betasplinejac(Xval[0],2,1) + (-Xval[2] + 1)*Alphasplinejac(Xval[0],2,1))*(0.333333333333333*hstep*(-Alphaspline(Xval[0], 2) - Betaspline(Xval[0], 2)) + 2) + 0.333333333333333*hstep*(-Alphasplinejac(Xval[0],2,1) - Betasplinejac(Xval[0],2,1))*(Xval[2] - Xvalp1[2] + 0.166666666666667*hstep*(-Xval[2]*Betaspline(Xval[0], 2) - 4*Xval2[2]*Betaspline(Xval2[0], 2) - Xvalp1[2]*Betaspline(Xvalp1[0], 2) + (-Xval[2] + 1)*Alphaspline(Xval[0], 2) + (-4*Xval2[2] + 4)*Alphaspline(Xval2[0], 2) + (-Xvalp1[2] + 1)*Alphaspline(Xvalp1[0], 2))) + 0.25*hstep*(-Alphasplinejac(Xval[0],2,1) - Betasplinejac(Xval[0],2,1))*(0.5*Xval[2] - Xval2[2] + 0.5*Xvalp1[2] + 0.125*hstep*(-Xval[2]*Betaspline(Xval[0], 2) + Xvalp1[2]*Betaspline(Xvalp1[0], 2) + (-Xval[2] + 1)*Alphaspline(Xval[0], 2) - (-Xvalp1[2] + 1)*Alphaspline(Xvalp1[0], 2))))/(2*Time+1);
   values[4*(Time+1)+4*(Time)+0+jt] += obj_factor*bounds[2][3]*(// Not supported in C:
// Alphaspline
// Betaspline
// Derivative
// Derivative
0.166666666666667*hstep*(-Xvalp1[2]*Betasplinejac(Xvalp1[0],2,1) + (-Xvalp1[2] + 1)*Alphasplinejac(Xvalp1[0],2,1))*(0.333333333333333*hstep*(-Alphaspline(Xval[0], 2) - Betaspline(Xval[0], 2)) + 2) + 0.125*hstep*(Xvalp1[2]*Betasplinejac(Xvalp1[0],2,1) - (-Xvalp1[2] + 1)*Alphasplinejac(Xvalp1[0],2,1))*(0.25*hstep*(-Alphaspline(Xval[0], 2) - Betaspline(Xval[0], 2)) + 1.0))/(2*Time+1);
   values[5*(Time+1)+6*(Time)+0+jt] += obj_factor*bounds[2][3]*(// Not supported in C:
// Alphaspline
// Betaspline
(0.125*hstep*(-Alphaspline(Xval[0], 2) - Betaspline(Xval[0], 2)) + 0.5)*(0.25*hstep*(-Alphaspline(Xval[0], 2) - Betaspline(Xval[0], 2)) + 1.0) + (0.166666666666667*hstep*(-Alphaspline(Xval[0], 2) - Betaspline(Xval[0], 2)) + 1)*(0.333333333333333*hstep*(-Alphaspline(Xval[0], 2) - Betaspline(Xval[0], 2)) + 2))/(2*Time+1);
   values[6*(Time+1)+6*(Time)+0+jt] += obj_factor*bounds[2][3]*(// Not supported in C:
// Alphaspline
// Betaspline
// Derivative
// Derivative
0.166666666666667*hstep*(-Xval[2]*Betasplinejac(Xval[0],2,1) + (-Xval[2] + 1)*Alphasplinejac(Xval[0],2,1))*(0.333333333333333*hstep*(-Alphaspline(Xvalp1[0], 2) - Betaspline(Xvalp1[0], 2)) - 2) + 0.125*hstep*(-Xval[2]*Betasplinejac(Xval[0],2,1) + (-Xval[2] + 1)*Alphasplinejac(Xval[0],2,1))*(0.25*hstep*(Alphaspline(Xvalp1[0], 2) + Betaspline(Xvalp1[0], 2)) + 1.0))/(2*Time+1);
   values[3*(Time+1)+4*(Time)+0+ 1+jt] += obj_factor*bounds[2][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
// Betaspline
// Derivative
// Derivative
0.166666666666667*hstep*(-Xvalp1[2]*Betasplinejac(Xvalp1[0],2,1) + (-Xvalp1[2] + 1)*Alphasplinejac(Xvalp1[0],2,1))*(0.333333333333333*hstep*(-Alphaspline(Xvalp1[0], 2) - Betaspline(Xvalp1[0], 2)) - 2) + 0.125*hstep*(Xvalp1[2]*Betasplinejac(Xvalp1[0],2,1) - (-Xvalp1[2] + 1)*Alphasplinejac(Xvalp1[0],2,1))*(0.25*hstep*(Alphaspline(Xvalp1[0], 2) + Betaspline(Xvalp1[0], 2)) + 1.0) + 0.333333333333333*hstep*(-Alphasplinejac(Xvalp1[0],2,1) - Betasplinejac(Xvalp1[0],2,1))*(Xval[2] - Xvalp1[2] + 0.166666666666667*hstep*(-Xval[2]*Betaspline(Xval[0], 2) - 4*Xval2[2]*Betaspline(Xval2[0], 2) - Xvalp1[2]*Betaspline(Xvalp1[0], 2) + (-Xval[2] + 1)*Alphaspline(Xval[0], 2) + (-4*Xval2[2] + 4)*Alphaspline(Xval2[0], 2) + (-Xvalp1[2] + 1)*Alphaspline(Xvalp1[0], 2))) + 0.25*hstep*(Alphasplinejac(Xvalp1[0],2,1) + Betasplinejac(Xvalp1[0],2,1))*(0.5*Xval[2] - Xval2[2] + 0.5*Xvalp1[2] + 0.125*hstep*(-Xval[2]*Betaspline(Xval[0], 2) + Xvalp1[2]*Betaspline(Xvalp1[0], 2) + (-Xval[2] + 1)*Alphaspline(Xval[0], 2) - (-Xvalp1[2] + 1)*Alphaspline(Xvalp1[0], 2))))/(2*Time+1);
   values[6*(Time+1)+8*(Time)+0+jt] += obj_factor*bounds[2][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
(0.125*hstep*(-Alphaspline(Xval[0], 2) - Betaspline(Xval[0], 2)) + 0.5)*(0.25*hstep*(Alphaspline(Xvalp1[0], 2) + Betaspline(Xvalp1[0], 2)) + 1.0) + (0.166666666666667*hstep*(-Alphaspline(Xval[0], 2) - Betaspline(Xval[0], 2)) + 1)*(0.333333333333333*hstep*(-Alphaspline(Xvalp1[0], 2) - Betaspline(Xvalp1[0], 2)) - 2))/(2*Time+1);
   values[5*(Time+1)+6*(Time)+0+ 1+jt] += obj_factor*bounds[2][3]*(// Not supported in C:
// Alphaspline
// Betaspline
(0.166666666666667*hstep*(-Alphaspline(Xvalp1[0], 2) - Betaspline(Xvalp1[0], 2)) - 1)*(0.333333333333333*hstep*(-Alphaspline(Xvalp1[0], 2) - Betaspline(Xvalp1[0], 2)) - 2) + (0.125*hstep*(Alphaspline(Xvalp1[0], 2) + Betaspline(Xvalp1[0], 2)) + 0.5)*(0.25*hstep*(Alphaspline(Xvalp1[0], 2) + Betaspline(Xvalp1[0], 2)) + 1.0))/(2*Time+1);
   values[10*(Time+1)+16*(Time)+0+jt] += obj_factor*bounds[2][3]*(// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
0.0555555555555556*pow(hstep, 2)*(-Xval[2]*Betasplinejac(Xval[0],2,1) + (-Xval[2] + 1)*Alphasplinejac(Xval[0],2,1))*(-4*Xval2[2]*Betasplinejac(Xval2[0],2,1) + (-4*Xval2[2] + 4)*Alphasplinejac(Xval2[0],2,1)))/(2*Time+1);
   values[10*(Time+1)+17*(Time)+0+jt] += obj_factor*bounds[2][3]*(// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
0.0555555555555556*pow(hstep, 2)*(-4*Xval2[2]*Betasplinejac(Xval2[0],2,1) + (-4*Xval2[2] + 4)*Alphasplinejac(Xval2[0],2,1))*(-Xvalp1[2]*Betasplinejac(Xvalp1[0],2,1) + (-Xvalp1[2] + 1)*Alphasplinejac(Xvalp1[0],2,1)))/(2*Time+1);
   values[10*(Time+1)+20*(Time)+0+jt] += obj_factor*bounds[2][3]*(// Not supported in C:
// Alphaspline
// Betaspline
// Derivative
// Derivative
0.333333333333333*hstep*(-4*Xval2[2]*Betasplinejac(Xval2[0],2,1) + (-4*Xval2[2] + 4)*Alphasplinejac(Xval2[0],2,1))*(0.166666666666667*hstep*(-Alphaspline(Xval[0], 2) - Betaspline(Xval[0], 2)) + 1))/(2*Time+1);
   values[10*(Time+1)+21*(Time)+0+jt] += obj_factor*bounds[2][3]*(// Not supported in C:
// Alphaspline
// Betaspline
// Derivative
// Derivative
0.333333333333333*hstep*(-4*Xval2[2]*Betasplinejac(Xval2[0],2,1) + (-4*Xval2[2] + 4)*Alphasplinejac(Xval2[0],2,1))*(0.166666666666667*hstep*(-Alphaspline(Xvalp1[0], 2) - Betaspline(Xvalp1[0], 2)) - 1))/(2*Time+1);
   values[10*(Time+1)+24*(Time)+0+jt] += obj_factor*bounds[2][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
// Betaspline
// Derivative
// Derivative
// Derivative
// Derivative
0.0555555555555556*pow(hstep, 2)*pow(-4*Xval2[2]*Betasplinejac(Xval2[0],2,1) + (-4*Xval2[2] + 4)*Alphasplinejac(Xval2[0],2,1), 2) + 0.333333333333333*hstep*(-4*Xval2[2]*Betasplinehes(Xval2[0],2,1,1) + (-4*Xval2[2] + 4)*Alphasplinehes(Xval2[0],2,1,1))*(Xval[2] - Xvalp1[2] + 0.166666666666667*hstep*(-Xval[2]*Betaspline(Xval[0], 2) - 4*Xval2[2]*Betaspline(Xval2[0], 2) - Xvalp1[2]*Betaspline(Xvalp1[0], 2) + (-Xval[2] + 1)*Alphaspline(Xval[0], 2) + (-4*Xval2[2] + 4)*Alphaspline(Xval2[0], 2) + (-Xvalp1[2] + 1)*Alphaspline(Xvalp1[0], 2))))/(2*Time+1);
   values[10*(Time+1)+35*(Time)+0+jt] += obj_factor*bounds[2][3]*(// Not supported in C:
// Alphaspline
// Betaspline
// Derivative
// Derivative
0.0555555555555556*pow(hstep, 2)*(-Xval[2]*Betasplinejac(Xval[0],2,1) + (-Xval[2] + 1)*Alphasplinejac(Xval[0],2,1))*(-4*Alphaspline(Xval2[0], 2) - 4*Betaspline(Xval2[0], 2)) - 0.25*hstep*(-Xval[2]*Betasplinejac(Xval[0],2,1) + (-Xval[2] + 1)*Alphasplinejac(Xval[0],2,1)))/(2*Time+1);
   values[10*(Time+1)+36*(Time)+0+jt] += obj_factor*bounds[2][3]*(// Not supported in C:
// Alphaspline
// Betaspline
// Derivative
// Derivative
0.0555555555555556*pow(hstep, 2)*(-Xvalp1[2]*Betasplinejac(Xvalp1[0],2,1) + (-Xvalp1[2] + 1)*Alphasplinejac(Xvalp1[0],2,1))*(-4*Alphaspline(Xval2[0], 2) - 4*Betaspline(Xval2[0], 2)) - 0.25*hstep*(Xvalp1[2]*Betasplinejac(Xvalp1[0],2,1) - (-Xvalp1[2] + 1)*Alphasplinejac(Xvalp1[0],2,1)))/(2*Time+1);
   values[10*(Time+1)+39*(Time)+0+jt] += obj_factor*bounds[2][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
0.333333333333333*hstep*(0.166666666666667*hstep*(-Alphaspline(Xval[0], 2) - Betaspline(Xval[0], 2)) + 1)*(-4*Alphaspline(Xval2[0], 2) - 4*Betaspline(Xval2[0], 2)) - 0.25*hstep*(-Alphaspline(Xval[0], 2) - Betaspline(Xval[0], 2)) - 1.0)/(2*Time+1);
   values[10*(Time+1)+40*(Time)+0+jt] += obj_factor*bounds[2][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
0.333333333333333*hstep*(0.166666666666667*hstep*(-Alphaspline(Xvalp1[0], 2) - Betaspline(Xvalp1[0], 2)) - 1)*(-4*Alphaspline(Xval2[0], 2) - 4*Betaspline(Xval2[0], 2)) - 0.25*hstep*(Alphaspline(Xvalp1[0], 2) + Betaspline(Xvalp1[0], 2)) - 1.0)/(2*Time+1);
   values[10*(Time+1)+43*(Time)+0+jt] += obj_factor*bounds[2][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
// Betaspline
// Derivative
// Derivative
0.0555555555555556*pow(hstep, 2)*(-4*Xval2[2]*Betasplinejac(Xval2[0],2,1) + (-4*Xval2[2] + 4)*Alphasplinejac(Xval2[0],2,1))*(-4*Alphaspline(Xval2[0], 2) - 4*Betaspline(Xval2[0], 2)) + 0.333333333333333*hstep*(-4*Alphasplinejac(Xval2[0],2,1) - 4*Betasplinejac(Xval2[0],2,1))*(Xval[2] - Xvalp1[2] + 0.166666666666667*hstep*(-Xval[2]*Betaspline(Xval[0], 2) - 4*Xval2[2]*Betaspline(Xval2[0], 2) - Xvalp1[2]*Betaspline(Xvalp1[0], 2) + (-Xval[2] + 1)*Alphaspline(Xval[0], 2) + (-4*Xval2[2] + 4)*Alphaspline(Xval2[0], 2) + (-Xvalp1[2] + 1)*Alphaspline(Xvalp1[0], 2))))/(2*Time+1);
   values[10*(Time+1)+45*(Time)+0+jt] += obj_factor*bounds[2][3]*(// Not supported in C:
// Alphaspline
// Betaspline
0.0555555555555556*pow(hstep, 2)*pow(-4*Alphaspline(Xval2[0], 2) - 4*Betaspline(Xval2[0], 2), 2) + 2)/(2*Time+1);
   values[0*(Time+1)+0*(Time)+0+jt] += obj_factor*bounds[3][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
// Betaspline
// Derivative
// Derivative
// Derivative
// Derivative
0.0868055555555556*pow(hstep, 2)*pow(-Xval[3]*Betasplinejac(Xval[0],3,1) + (-Xval[3] + 1)*Alphasplinejac(Xval[0],3,1), 2) + 0.333333333333333*hstep*(-Xval[3]*Betasplinehes(Xval[0],3,1,1) + (-Xval[3] + 1)*Alphasplinehes(Xval[0],3,1,1))*(0.166666666666667*hstep*(-Xval[3]*Betaspline(Xval[0], 3) - 4*Xval2[3]*Betaspline(Xval2[0], 3) - Xvalp1[3]*Betaspline(Xvalp1[0], 3) + (-Xval[3] + 1)*Alphaspline(Xval[0], 3) + (-4*Xval2[3] + 4)*Alphaspline(Xval2[0], 3) + (-Xvalp1[3] + 1)*Alphaspline(Xvalp1[0], 3)) + Xval[3] - Xvalp1[3]) + 0.25*hstep*(-Xval[3]*Betasplinehes(Xval[0],3,1,1) + (-Xval[3] + 1)*Alphasplinehes(Xval[0],3,1,1))*(0.125*hstep*(-Xval[3]*Betaspline(Xval[0], 3) + Xvalp1[3]*Betaspline(Xvalp1[0], 3) + (-Xval[3] + 1)*Alphaspline(Xval[0], 3) - (-Xvalp1[3] + 1)*Alphaspline(Xvalp1[0], 3)) + 0.5*Xval[3] - Xval2[3] + 0.5*Xvalp1[3]))/(2*Time+1);
   values[1*(Time+1)+0*(Time)+0+jt] += obj_factor*bounds[3][3]*(// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
0.0555555555555556*pow(hstep, 2)*(-Xval[3]*Betasplinejac(Xval[0],3,1) + (-Xval[3] + 1)*Alphasplinejac(Xval[0],3,1))*(-Xvalp1[3]*Betasplinejac(Xvalp1[0],3,1) + (-Xvalp1[3] + 1)*Alphasplinejac(Xvalp1[0],3,1)) + 0.03125*pow(hstep, 2)*(-Xval[3]*Betasplinejac(Xval[0],3,1) + (-Xval[3] + 1)*Alphasplinejac(Xval[0],3,1))*(Xvalp1[3]*Betasplinejac(Xvalp1[0],3,1) - (-Xvalp1[3] + 1)*Alphasplinejac(Xvalp1[0],3,1)))/(2*Time+1);
   values[0*(Time+1)+0*(Time)+0+ 1+jt] += obj_factor*bounds[3][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
// Betaspline
// Derivative
// Derivative
// Derivative
// Derivative
0.0555555555555556*pow(hstep, 2)*pow(-Xvalp1[3]*Betasplinejac(Xvalp1[0],3,1) + (-Xvalp1[3] + 1)*Alphasplinejac(Xvalp1[0],3,1), 2) + 0.03125*pow(hstep, 2)*pow(Xvalp1[3]*Betasplinejac(Xvalp1[0],3,1) - (-Xvalp1[3] + 1)*Alphasplinejac(Xvalp1[0],3,1), 2) + 0.333333333333333*hstep*(-Xvalp1[3]*Betasplinehes(Xvalp1[0],3,1,1) + (-Xvalp1[3] + 1)*Alphasplinehes(Xvalp1[0],3,1,1))*(0.166666666666667*hstep*(-Xval[3]*Betaspline(Xval[0], 3) - 4*Xval2[3]*Betaspline(Xval2[0], 3) - Xvalp1[3]*Betaspline(Xvalp1[0], 3) + (-Xval[3] + 1)*Alphaspline(Xval[0], 3) + (-4*Xval2[3] + 4)*Alphaspline(Xval2[0], 3) + (-Xvalp1[3] + 1)*Alphaspline(Xvalp1[0], 3)) + Xval[3] - Xvalp1[3]) + 0.25*hstep*(Xvalp1[3]*Betasplinehes(Xvalp1[0],3,1,1) - (-Xvalp1[3] + 1)*Alphasplinehes(Xvalp1[0],3,1,1))*(0.125*hstep*(-Xval[3]*Betaspline(Xval[0], 3) + Xvalp1[3]*Betaspline(Xvalp1[0], 3) + (-Xval[3] + 1)*Alphaspline(Xval[0], 3) - (-Xvalp1[3] + 1)*Alphaspline(Xvalp1[0], 3)) + 0.5*Xval[3] - Xval2[3] + 0.5*Xvalp1[3]))/(2*Time+1);
   values[6*(Time+1)+9*(Time)+0+jt] += obj_factor*bounds[3][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
// Betaspline
// Derivative
// Derivative
0.125*hstep*(0.25*hstep*(-Alphaspline(Xval[0], 3) - Betaspline(Xval[0], 3)) + 1.0)*(-Xval[3]*Betasplinejac(Xval[0],3,1) + (-Xval[3] + 1)*Alphasplinejac(Xval[0],3,1)) + 0.166666666666667*hstep*(0.333333333333333*hstep*(-Alphaspline(Xval[0], 3) - Betaspline(Xval[0], 3)) + 2)*(-Xval[3]*Betasplinejac(Xval[0],3,1) + (-Xval[3] + 1)*Alphasplinejac(Xval[0],3,1)) + 0.333333333333333*hstep*(-Alphasplinejac(Xval[0],3,1) - Betasplinejac(Xval[0],3,1))*(0.166666666666667*hstep*(-Xval[3]*Betaspline(Xval[0], 3) - 4*Xval2[3]*Betaspline(Xval2[0], 3) - Xvalp1[3]*Betaspline(Xvalp1[0], 3) + (-Xval[3] + 1)*Alphaspline(Xval[0], 3) + (-4*Xval2[3] + 4)*Alphaspline(Xval2[0], 3) + (-Xvalp1[3] + 1)*Alphaspline(Xvalp1[0], 3)) + Xval[3] - Xvalp1[3]) + 0.25*hstep*(-Alphasplinejac(Xval[0],3,1) - Betasplinejac(Xval[0],3,1))*(0.125*hstep*(-Xval[3]*Betaspline(Xval[0], 3) + Xvalp1[3]*Betaspline(Xvalp1[0], 3) + (-Xval[3] + 1)*Alphaspline(Xval[0], 3) - (-Xvalp1[3] + 1)*Alphaspline(Xvalp1[0], 3)) + 0.5*Xval[3] - Xval2[3] + 0.5*Xvalp1[3]))/(2*Time+1);
   values[7*(Time+1)+9*(Time)+0+jt] += obj_factor*bounds[3][3]*(// Not supported in C:
// Alphaspline
// Betaspline
// Derivative
// Derivative
0.125*hstep*(0.25*hstep*(-Alphaspline(Xval[0], 3) - Betaspline(Xval[0], 3)) + 1.0)*(Xvalp1[3]*Betasplinejac(Xvalp1[0],3,1) - (-Xvalp1[3] + 1)*Alphasplinejac(Xvalp1[0],3,1)) + 0.166666666666667*hstep*(0.333333333333333*hstep*(-Alphaspline(Xval[0], 3) - Betaspline(Xval[0], 3)) + 2)*(-Xvalp1[3]*Betasplinejac(Xvalp1[0],3,1) + (-Xvalp1[3] + 1)*Alphasplinejac(Xvalp1[0],3,1)))/(2*Time+1);
   values[9*(Time+1)+12*(Time)+0+jt] += obj_factor*bounds[3][3]*(// Not supported in C:
// Alphaspline
// Betaspline
(0.125*hstep*(-Alphaspline(Xval[0], 3) - Betaspline(Xval[0], 3)) + 0.5)*(0.25*hstep*(-Alphaspline(Xval[0], 3) - Betaspline(Xval[0], 3)) + 1.0) + (0.166666666666667*hstep*(-Alphaspline(Xval[0], 3) - Betaspline(Xval[0], 3)) + 1)*(0.333333333333333*hstep*(-Alphaspline(Xval[0], 3) - Betaspline(Xval[0], 3)) + 2))/(2*Time+1);
   values[10*(Time+1)+12*(Time)+0+jt] += obj_factor*bounds[3][3]*(// Not supported in C:
// Alphaspline
// Betaspline
// Derivative
// Derivative
0.166666666666667*hstep*(0.333333333333333*hstep*(-Alphaspline(Xvalp1[0], 3) - Betaspline(Xvalp1[0], 3)) - 2)*(-Xval[3]*Betasplinejac(Xval[0],3,1) + (-Xval[3] + 1)*Alphasplinejac(Xval[0],3,1)) + 0.125*hstep*(0.25*hstep*(Alphaspline(Xvalp1[0], 3) + Betaspline(Xvalp1[0], 3)) + 1.0)*(-Xval[3]*Betasplinejac(Xval[0],3,1) + (-Xval[3] + 1)*Alphasplinejac(Xval[0],3,1)))/(2*Time+1);
   values[6*(Time+1)+9*(Time)+0+ 1+jt] += obj_factor*bounds[3][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
// Betaspline
// Derivative
// Derivative
0.166666666666667*hstep*(0.333333333333333*hstep*(-Alphaspline(Xvalp1[0], 3) - Betaspline(Xvalp1[0], 3)) - 2)*(-Xvalp1[3]*Betasplinejac(Xvalp1[0],3,1) + (-Xvalp1[3] + 1)*Alphasplinejac(Xvalp1[0],3,1)) + 0.125*hstep*(0.25*hstep*(Alphaspline(Xvalp1[0], 3) + Betaspline(Xvalp1[0], 3)) + 1.0)*(Xvalp1[3]*Betasplinejac(Xvalp1[0],3,1) - (-Xvalp1[3] + 1)*Alphasplinejac(Xvalp1[0],3,1)) + 0.333333333333333*hstep*(-Alphasplinejac(Xvalp1[0],3,1) - Betasplinejac(Xvalp1[0],3,1))*(0.166666666666667*hstep*(-Xval[3]*Betaspline(Xval[0], 3) - 4*Xval2[3]*Betaspline(Xval2[0], 3) - Xvalp1[3]*Betaspline(Xvalp1[0], 3) + (-Xval[3] + 1)*Alphaspline(Xval[0], 3) + (-4*Xval2[3] + 4)*Alphaspline(Xval2[0], 3) + (-Xvalp1[3] + 1)*Alphaspline(Xvalp1[0], 3)) + Xval[3] - Xvalp1[3]) + 0.25*hstep*(Alphasplinejac(Xvalp1[0],3,1) + Betasplinejac(Xvalp1[0],3,1))*(0.125*hstep*(-Xval[3]*Betaspline(Xval[0], 3) + Xvalp1[3]*Betaspline(Xvalp1[0], 3) + (-Xval[3] + 1)*Alphaspline(Xval[0], 3) - (-Xvalp1[3] + 1)*Alphaspline(Xvalp1[0], 3)) + 0.5*Xval[3] - Xval2[3] + 0.5*Xvalp1[3]))/(2*Time+1);
   values[10*(Time+1)+15*(Time)+0+jt] += obj_factor*bounds[3][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
(0.125*hstep*(-Alphaspline(Xval[0], 3) - Betaspline(Xval[0], 3)) + 0.5)*(0.25*hstep*(Alphaspline(Xvalp1[0], 3) + Betaspline(Xvalp1[0], 3)) + 1.0) + (0.166666666666667*hstep*(-Alphaspline(Xval[0], 3) - Betaspline(Xval[0], 3)) + 1)*(0.333333333333333*hstep*(-Alphaspline(Xvalp1[0], 3) - Betaspline(Xvalp1[0], 3)) - 2))/(2*Time+1);
   values[9*(Time+1)+12*(Time)+0+ 1+jt] += obj_factor*bounds[3][3]*(// Not supported in C:
// Alphaspline
// Betaspline
(0.166666666666667*hstep*(-Alphaspline(Xvalp1[0], 3) - Betaspline(Xvalp1[0], 3)) - 1)*(0.333333333333333*hstep*(-Alphaspline(Xvalp1[0], 3) - Betaspline(Xvalp1[0], 3)) - 2) + (0.125*hstep*(Alphaspline(Xvalp1[0], 3) + Betaspline(Xvalp1[0], 3)) + 0.5)*(0.25*hstep*(Alphaspline(Xvalp1[0], 3) + Betaspline(Xvalp1[0], 3)) + 1.0))/(2*Time+1);
   values[10*(Time+1)+16*(Time)+0+jt] += obj_factor*bounds[3][3]*(// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
0.0555555555555556*pow(hstep, 2)*(-Xval[3]*Betasplinejac(Xval[0],3,1) + (-Xval[3] + 1)*Alphasplinejac(Xval[0],3,1))*(-4*Xval2[3]*Betasplinejac(Xval2[0],3,1) + (-4*Xval2[3] + 4)*Alphasplinejac(Xval2[0],3,1)))/(2*Time+1);
   values[10*(Time+1)+17*(Time)+0+jt] += obj_factor*bounds[3][3]*(// Not supported in C:
// Derivative
// Derivative
// Derivative
// Derivative
0.0555555555555556*pow(hstep, 2)*(-4*Xval2[3]*Betasplinejac(Xval2[0],3,1) + (-4*Xval2[3] + 4)*Alphasplinejac(Xval2[0],3,1))*(-Xvalp1[3]*Betasplinejac(Xvalp1[0],3,1) + (-Xvalp1[3] + 1)*Alphasplinejac(Xvalp1[0],3,1)))/(2*Time+1);
   values[10*(Time+1)+22*(Time)+0+jt] += obj_factor*bounds[3][3]*(// Not supported in C:
// Alphaspline
// Betaspline
// Derivative
// Derivative
0.333333333333333*hstep*(0.166666666666667*hstep*(-Alphaspline(Xval[0], 3) - Betaspline(Xval[0], 3)) + 1)*(-4*Xval2[3]*Betasplinejac(Xval2[0],3,1) + (-4*Xval2[3] + 4)*Alphasplinejac(Xval2[0],3,1)))/(2*Time+1);
   values[10*(Time+1)+23*(Time)+0+jt] += obj_factor*bounds[3][3]*(// Not supported in C:
// Alphaspline
// Betaspline
// Derivative
// Derivative
0.333333333333333*hstep*(0.166666666666667*hstep*(-Alphaspline(Xvalp1[0], 3) - Betaspline(Xvalp1[0], 3)) - 1)*(-4*Xval2[3]*Betasplinejac(Xval2[0],3,1) + (-4*Xval2[3] + 4)*Alphasplinejac(Xval2[0],3,1)))/(2*Time+1);
   values[10*(Time+1)+24*(Time)+0+jt] += obj_factor*bounds[3][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
// Betaspline
// Derivative
// Derivative
// Derivative
// Derivative
0.0555555555555556*pow(hstep, 2)*pow(-4*Xval2[3]*Betasplinejac(Xval2[0],3,1) + (-4*Xval2[3] + 4)*Alphasplinejac(Xval2[0],3,1), 2) + 0.333333333333333*hstep*(-4*Xval2[3]*Betasplinehes(Xval2[0],3,1,1) + (-4*Xval2[3] + 4)*Alphasplinehes(Xval2[0],3,1,1))*(0.166666666666667*hstep*(-Xval[3]*Betaspline(Xval[0], 3) - 4*Xval2[3]*Betaspline(Xval2[0], 3) - Xvalp1[3]*Betaspline(Xvalp1[0], 3) + (-Xval[3] + 1)*Alphaspline(Xval[0], 3) + (-4*Xval2[3] + 4)*Alphaspline(Xval2[0], 3) + (-Xvalp1[3] + 1)*Alphaspline(Xvalp1[0], 3)) + Xval[3] - Xvalp1[3]))/(2*Time+1);
   values[10*(Time+1)+46*(Time)+0+jt] += obj_factor*bounds[3][3]*(// Not supported in C:
// Alphaspline
// Betaspline
// Derivative
// Derivative
0.0555555555555556*pow(hstep, 2)*(-Xval[3]*Betasplinejac(Xval[0],3,1) + (-Xval[3] + 1)*Alphasplinejac(Xval[0],3,1))*(-4*Alphaspline(Xval2[0], 3) - 4*Betaspline(Xval2[0], 3)) - 0.25*hstep*(-Xval[3]*Betasplinejac(Xval[0],3,1) + (-Xval[3] + 1)*Alphasplinejac(Xval[0],3,1)))/(2*Time+1);
   values[10*(Time+1)+47*(Time)+0+jt] += obj_factor*bounds[3][3]*(// Not supported in C:
// Alphaspline
// Betaspline
// Derivative
// Derivative
0.0555555555555556*pow(hstep, 2)*(-Xvalp1[3]*Betasplinejac(Xvalp1[0],3,1) + (-Xvalp1[3] + 1)*Alphasplinejac(Xvalp1[0],3,1))*(-4*Alphaspline(Xval2[0], 3) - 4*Betaspline(Xval2[0], 3)) - 0.25*hstep*(Xvalp1[3]*Betasplinejac(Xvalp1[0],3,1) - (-Xvalp1[3] + 1)*Alphasplinejac(Xvalp1[0],3,1)))/(2*Time+1);
   values[10*(Time+1)+52*(Time)+0+jt] += obj_factor*bounds[3][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
0.333333333333333*hstep*(0.166666666666667*hstep*(-Alphaspline(Xval[0], 3) - Betaspline(Xval[0], 3)) + 1)*(-4*Alphaspline(Xval2[0], 3) - 4*Betaspline(Xval2[0], 3)) - 0.25*hstep*(-Alphaspline(Xval[0], 3) - Betaspline(Xval[0], 3)) - 1.0)/(2*Time+1);
   values[10*(Time+1)+53*(Time)+0+jt] += obj_factor*bounds[3][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
0.333333333333333*hstep*(0.166666666666667*hstep*(-Alphaspline(Xvalp1[0], 3) - Betaspline(Xvalp1[0], 3)) - 1)*(-4*Alphaspline(Xval2[0], 3) - 4*Betaspline(Xval2[0], 3)) - 0.25*hstep*(Alphaspline(Xvalp1[0], 3) + Betaspline(Xvalp1[0], 3)) - 1.0)/(2*Time+1);
   values[10*(Time+1)+54*(Time)+0+jt] += obj_factor*bounds[3][3]*(// Not supported in C:
// Alphaspline
// Alphaspline
// Alphaspline
// Betaspline
// Betaspline
// Betaspline
// Derivative
// Derivative
0.0555555555555556*pow(hstep, 2)*(-4*Xval2[3]*Betasplinejac(Xval2[0],3,1) + (-4*Xval2[3] + 4)*Alphasplinejac(Xval2[0],3,1))*(-4*Alphaspline(Xval2[0], 3) - 4*Betaspline(Xval2[0], 3)) + 0.333333333333333*hstep*(-4*Alphasplinejac(Xval2[0],3,1) - 4*Betasplinejac(Xval2[0],3,1))*(0.166666666666667*hstep*(-Xval[3]*Betaspline(Xval[0], 3) - 4*Xval2[3]*Betaspline(Xval2[0], 3) - Xvalp1[3]*Betaspline(Xvalp1[0], 3) + (-Xval[3] + 1)*Alphaspline(Xval[0], 3) + (-4*Xval2[3] + 4)*Alphaspline(Xval2[0], 3) + (-Xvalp1[3] + 1)*Alphaspline(Xvalp1[0], 3)) + Xval[3] - Xvalp1[3]))/(2*Time+1);
   values[10*(Time+1)+57*(Time)+0+jt] += obj_factor*bounds[3][3]*(// Not supported in C:
// Alphaspline
// Betaspline
0.0555555555555556*pow(hstep, 2)*pow(-4*Alphaspline(Xval2[0], 3) - 4*Betaspline(Xval2[0], 3), 2) + 2)/(2*Time+1);
   } // end for loop 

// Add last element
     for(Index i=0;i<nY;i++) {
        Xval[i] = x[Time + i*(Time+1)];
        Xvalp1[i] = 0;
        Xval2[i] = 0;
     } //end for loop

     for(Index cup=0;cup<nU;cup++) {
        K11val[cup] = x[Time + nY*(Time+1) + cup*(Time+1)];
        K11valp1[cup] = 0;
        K11val2[cup] = 0;
     } //end for loop

     Xdval[0] = VDATA0[2*Time];
     Xdval2[0] = 0;
     Xdvalp1[0] = 0;
     Xdval[1] = MDATA0[2*Time];
     Xdval2[1] = 0;
     Xdvalp1[1] = 0;
     Xdval[2] = HDATA0[2*Time];
     Xdval2[2] = 0;
     Xdvalp1[2] = 0;
     Xdval[3] = NDATA0[2*Time];
     Xdval2[3] = 0;
     Xdvalp1[3] = 0;
     Ival[0] = Iinj[2*Time];
     Ival2[0] = 0;
     Ivalp1[0] = 0;

     for(Index i=0;i<nP;i++) {
        Pval[i] = x[(2*Time+1)*(nY+nU)+i];
     } //end for loop


   values[0*(Time+1)+0*(Time)+0+ Time] += obj_factor*(2)/(2*Time+1);
   values[2*(Time+1)+2*(Time)+0+ Time] += obj_factor*(2)/(2*Time+1);
   values[5*(Time+1)+6*(Time)+0+ Time] += obj_factor*(2)/(2*Time+1);
   values[9*(Time+1)+12*(Time)+0+ Time] += obj_factor*(2)/(2*Time+1);
  } // end else 

   return true;
}


void NEURODYN_NAKL_NLP::finalize_solution(SolverReturn status,
                        Index n, const Number* x, const Number* z_L, const Number* z_U,
                        Index m, const Number* g, const Number* lambda,
                        Number obj_value,
                        const IpoptData* ip_data,
                        IpoptCalculatedQuantities* ip_cq)
{
  FILE *OUTPUT1;
  char filename[20];
  if ((file_fmt < 2) || (modid == 0)){
    sprintf(filename,"D%d_M%d_IC%d.dat", nY,nM,taskid);
  }
  else if ((file_fmt >= 2) && (modid != 0)){
    sprintf(filename,"D%d_M%d_PATH%d_IC%d.dat", nY,nM,pathid,taskid);
  }
  if(beta==0){
    OUTPUT1 = fopen (filename,"w");
  }
  else{
    if ((output_fmt == 1) || (output_fmt < 0)){
      OUTPUT1 = fopen (filename,"w");
    }
    else 
      OUTPUT1 = fopen (filename,"a");
  }
  // Write solution for annealing
  for(Index i=0;i<Ntotal-nP;i++){
    solution[i] = x[i];
  }
  for(Index i=0;i<nP;i++){
    solution[Ntotal-nP+i] = x[(2*Time+1)*(nU+nY)+i];
  }

  // Write to file
  if (output_fmt != 1){
    fprintf(OUTPUT1, "%d %d %e ",beta, (status == SUCCESS), obj_value);
  }
  for (Index i=0;i<Time;i++) {
    for (Index j=0;j<nY;j++) {
      fprintf(OUTPUT1,"%e ", x[j*(Time+1)+i]);
    }
    if (abs(output_fmt) == 2){
      for (Index cup=0;cup<nU;cup++) {
        fprintf(OUTPUT1,"%e ", x[nY*(Time+1)+ cup*(Time+1)+ i]);
      }
    }
    for (Index j=0;j<nY;j++) {
      fprintf(OUTPUT1,"%e ", x[(nY+nU)*(Time+1) + j*Time + i]);
    }
    if (abs(output_fmt) == 2){
      for (Index cup=0;cup<nU;cup++) {
        fprintf(OUTPUT1,"%e ", x[(nY+nU)*(Time+1) + nY*Time + cup*Time + i]);
      }
    }
  }
  for (Index j=0;j<nY;j++) {
     fprintf(OUTPUT1,"%e ", x[j*(Time+1) + Time]);
  }
  if (abs(output_fmt) == 2){
    for (Index cup=0;cup<nU;cup++) {
       fprintf(OUTPUT1,"%e ", x[nY*(Time+1) + cup*(Time+1) + Time]);
    }
  }
  for (Index j=0;j<nP;j++) {
     fprintf(OUTPUT1,"%e ", x[(2*Time+1)*(nY+nU)+j]);
  }
  fprintf(OUTPUT1,"\n");
  fclose (OUTPUT1);
  

}
