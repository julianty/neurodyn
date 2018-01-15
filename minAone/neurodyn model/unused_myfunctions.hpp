#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>

using namespace std;

double sigma(double V, double vBiask, int sign, int V_scale)
{
  double mu = 0.7;
  double Ut = 0.026 * V_scale;
  value = 1 / (1 + exp(sign * mu * (vBiask - V) / Ut)); 
  return value;
}

double sigmajac(double V, double vBias, int sign, int Vscale, int n)
{
  double value = 0.0;
  switch(n) { 
  case 1:
    value = sigma(V, vBias, sign, Vscale) * (1 - sigma(V, vBias, sign, Vscale)); break;
  
  case 2:
    value = 0; break;

  case 3:
    value = 0; break;

  case 4:
    value = 0; break;    

  default:
    cout << "Error in sigma jac"; break;
    value = 0;
  }
  return value;
}

double sigmahes(double V, double vBias, int sign, int Vscale, int n, int m)
{
  double value = 0.0;
  if(m == 1){
    switch(n) {
    case 1:
      value = sigmajac(V, vBias, sign, Vscale, n) * (1 - sigma(V, vBias, sign, Vscale)) + sigma(V, vBias, sign, Vscale) * -sigmajac(V, vBias, sign, Vscale, n); break;
    case 2:
      value = 0; break;
    case 3:
      value = 0; break;
    case 4:
      value = 0; break;
    default:
      value = 0;
      cout << "Error in sigmahes"; break;
    }
    return value;
  }
  else if(m==2){
    return 0;
  }
  else {
    cout << "Error in sigmahes";
    return 0;
  }
}

double alphaspline(double V, int x)
{
  //value of x determines whether m, h, n
  double alpha = 0.0;
  double vBias[7] = {0.6352, 0.7568, 0.8784, 1.0000, 1.1215, 1.2431, 1.3647};
  int alpha_arr[3][7] = {
    {0, 0, 1126, 3756, 7512, 9605, 9605},
    {2225, 751, 0, 0 ,0 ,0 ,0},
    {0, 0, 0, 0, 169, 47, 403}
  };
  for(int i=0; i<=6; i++){
    if(x == 1){
      alpha += alpha_arr[0][i] * sigma(V, vBias[i], 1, 1);
    }
    if(x == 2){
      alpha += alpha_arr[1][i] * sigma(V, vBias[i], -1, 1);
    }
    if(x == 3){
      alpha += alpha_arr[2][i] * sigma(V, vBias[i], 1, 1);
    }
  }
    return alpha;
}

double alphasplinejac(double V, int x, int n)
{
  double alpha = 0.0;
  double vBias[7] = {0.6352, 0.7568, 0.8784, 1.0000, 1.1215, 1.2431, 1.3647};
  int alpha_arr[3][7] = {
    {0, 0, 1126, 3756, 7512, 9605, 9605},
    {2225, 751, 0, 0 ,0 ,0 ,0},
    {0, 0, 0, 0, 169, 47, 403}
  };
  if(n==1){
    for(int i=0; i<=6; i++){
      if(x == 1){
	alpha += alpha_arr[0][i] * sigmajac(V, vBias[i], 1, 1, n);
      }
      if(x == 2){
	alpha += alpha_arr[1][i] * sigmajac(V, vBias[i], -1, 1, n);
      }
      if(x == 3){
	alpha += alpha_arr[2][i] * sigmajac(V, vBias[i], 1, 1, n);
      }
    }
    return alpha;
  }
  else if(n==2){
    return 0;
  }
  else {
    cout << "Error in alphasplinejac";
    return 0;
  }
}

double alphasplinehes(double V, int x, int n, int m)
{
  double alpha = 0.0;
  double vBias[7] = {0.6352, 0.7568, 0.8784, 1.0000, 1.1215, 1.2431, 1.3647};
  int alpha_arr[3][7] = {
    {0, 0, 1126, 3756, 7512, 9605, 9605},
    {2225, 751, 0, 0 ,0 ,0 ,0},
    {0, 0, 0, 0, 169, 47, 403}
  };
  if(m==1){ 
    if(n==1){
      for(int i=0; i<=6; i++){
	if(x == 1){
	  alpha += alpha_arr[0][i] * sigmahes(V, vBias[i], 1, 1, n, m);
	}
	if(x == 2){
	  alpha += alpha_arr[1][i] * sigmahes(V, vBias[i], -1, 1, n, m);
	}
	if(x == 3){
	  alpha += alpha_arr[2][i] * sigmahes(V, vBias[i], 1, 1, n, m);
	}
      }
      return alpha;
    }
    else if(n==2){
      return 0;
    }
    else {
      cout << "Error in alphasplinehes";
      return 0;
    }
  }
  else if(m==2){
    return 0;
  }
  else {
    cout << "Error in alphasplinehes";
    return 0;
  }
}

double betaspline(double V, int x)
{
  double beta = 0.0;
  double vBias[7] = {0.6352, 0.7568, 0.8784, 1.0000, 1.1215, 1.2431, 1.3647};
  int beta_arr[3][7] = {
    {9605, 9605, 9605, 9605, 0, 0, 0},
    {0, 0, 0, 0, 385, 469, 657},
    {9, 0, 0, 0, 0, 0, 9}
  };
  for(int i=0; i<=6; i++){
    if(x == 1){
      beta += beta_arr[0][i] * sigma(V, vBias[i], -1, 1);
    }
    if(x == 2){
      beta += beta_arr[1][i] * sigma(V, vBias[i], 1, 1);
    }
    if(x == 3){
      beta += beta_arr[2][i] * sigma(V, vBias[i], -1, 1);
    }
  }
  return beta;
}

double betasplinejac(double V, int x, int n)
{
  double beta = 0.0;
  double vBias[7] = {0.6352, 0.7568, 0.8784, 1.0000, 1.1215, 1.2431, 1.3647};
  int beta_arr[3][7] = {
    {9605, 9605, 9605, 9605, 0, 0, 0},
    {0, 0, 0, 0, 385, 469, 657},
    {9, 0, 0, 0, 0, 0, 9}
  };
  if(n==1){
    for(int i=0; i<=6; i++){
      if(x == 1){
	beta += beta_arr[0][i] * sigmajac(V, vBias[i], -1, 1, n);
      }
      if(x == 2){
	beta += beta_arr[1][i] * sigmajac(V, vBias[i], 1, 1, n);
      }
      if(x == 3){
	beta += beta_arr[2][i] * sigmajac(V, vBias[i], -1, 1, n);
      }
    }
    return beta;
  }
  else if(n==2){
    return 0;
  }
  else {
    cout << "Error in betasplinejac";
    return 0;
  }
}

double betasplinehes(double V, int x, int n, int m)
{
  double beta = 0.0;
  double vBias[7] = {0.6352, 0.7568, 0.8784, 1.0000, 1.1215, 1.2431, 1.3647};
  int beta_arr[3][7] = {
    {9605, 9605, 9605, 9605, 0, 0, 0},
    {0, 0, 0, 0, 385, 469, 657},
    {9, 0, 0, 0, 0, 0, 9}
  };
  if(m==1){
    if(n==1){
      for(int i=0; i<=6; i++){
	if(x == 1){
	  beta += beta_arr[0][i] * sigmahes(V, vBias[i], -1, 1, n, m);
	}
	if(x == 2){
	  beta += beta_arr[1][i] * sigmahes(V, vBias[i], 1, 1, n, m);
	}
	if(x == 3){
	  beta += beta_arr[2][i] * sigmahes(V, vBias[i], -1, 1, n, m);
	}
      }
      return beta;
    }
    else if(n==2){
      return 0;
    }
    else{
      cout << "Error in betasplinehes";
      return 0;
    }
  }
  else if(m==2){
    return 0;
  }
  else {
    cout <<"Error in betasplinehes";
    return 0;
  }
}
