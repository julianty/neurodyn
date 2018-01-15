#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>

using namespace std;

double Alphaspline(double V, int x)
{
  //value of x determines whether m, h, n
  double alpha = 0.0;
  double vBias[7] = {0.6352, 0.7568, 0.8784, 1.0000, 1.1215, 1.2431, 1.3647};
  int alpha_arr[3][7] = {
    {0, 0, 1126, 3756, 7512, 9605, 9605},
    {2225, 751, 0, 0 ,0 ,0 ,0},
    {0, 0, 0, 0, 169, 47, 403}
  };
  double sigma = 0.0;
  for(int i=0; i<=6; i++){
    if(x == 1){
      sigma = 1 / (1 + exp(-0.7 * (vBias[i] - V) / 0.026));
      alpha += alpha_arr[0][i] * sigma;
    }
    if(x == 2){
      sigma = 1 / (1 + exp(0.7 * (vBias[i] - V) / 0.026));
      alpha += alpha_arr[1][i] * sigma;
    }
    if(x == 3){
      sigma = 1 / (1 + exp(-0.7 * (vBias[i] - V) / 0.026));
      alpha += alpha_arr[2][i] * sigma;
    }xf
  }
    return alpha;
}

double Alphasplinejac(double V, int x, int n)
{
  double alpha = 0.0;
  double vBias[7] = {0.6352, 0.7568, 0.8784, 1.0000, 1.1215, 1.2431, 1.3647};
  int alpha_arr[3][7] = {
    {0, 0, 1126, 3756, 7512, 9605, 9605},
    {2225, 751, 0, 0 ,0 ,0 ,0},
    {0, 0, 0, 0, 169, 47, 403}
  };
  double sigma = 0.0;
  if(n==1){
    for(int i=0; i<=6; i++){
      if(x == 1){
	sigma = 1 / (1 + exp(-0.7 * (vBias[i] - V) / 0.026));
	alpha += alpha_arr[0][i] * sigma * (1 - sigma);
      }
      if(x == 2){
	sigma = 1 / (1 + exp(0.7 * (vBias[i] - V) / 0.026));
	alpha += alpha_arr[1][i] * sigma * (1 - sigma);
      }
      if(x == 3){
	sigma = 1 / (1 + exp(-0.7 * (vBias[i] - V) / 0.026));       
	alpha += alpha_arr[2][i] * sigma * (1 - sigma);
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

double Alphasplinehes(double V, int x, int n, int m)
{
  double alpha = 0.0;
  double vBias[7] = {0.6352, 0.7568, 0.8784, 1.0000, 1.1215, 1.2431, 1.3647};
  int alpha_arr[3][7] = {
    {0, 0, 1126, 3756, 7512, 9605, 9605},
    {2225, 751, 0, 0 ,0 ,0 ,0},
    {0, 0, 0, 0, 169, 47, 403}
  };
  double sigma = 0.0;
  if(m==1){ 
    if(n==1){
      for(int i=0; i<=6; i++){
	if(x == 1){
	  sigma = 1 / (1 + exp(-0.7 * (vBias[i] - V) / 0.026));
	  alpha += alpha_arr[0][i] * sigma - 3 * sigma*sigma + 2 * sigma*sigma*sigma;
	}
	if(x == 2){
	  sigma = 1 / (1 + exp(0.7 * (vBias[i] - V) / 0.026));	  
	  alpha += alpha_arr[1][i] * sigma - 3 * sigma*sigma + 2 * sigma*sigma*sigma;
	}
	if(x == 3){
	  sigma = 1 / (1 + exp(-0.7 * (vBias[i] - V) / 0.026));
	  alpha += alpha_arr[2][i] * sigma - 3 * sigma*sigma + 2 * sigma*sigma*sigma;
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