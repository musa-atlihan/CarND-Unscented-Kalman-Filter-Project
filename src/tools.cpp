#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */

  VectorXd rmse(4);
  rmse << 0,0,0,0;

  //  the estimation vector size should not be zero
  if (estimations.size() == 0)
  {
    cout<<"Error - estimations vector is empty.";
    return rmse;
  }

  // check if vector sizes are equal
  if (estimations.size() != ground_truth.size())
  {
    cout<<"Error - estimations vector size is not equal to ground_truth vector.";
    return rmse;
  }

  // accumulate squared residuals
  for(int i=0; i < estimations.size(); ++i){
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array()*residual.array();
    rmse += residual;
  }

  // calculate the mean
  rmse = rmse/estimations.size();

  // calculate the squared root
  rmse = rmse.array().sqrt();

  return rmse;
}
