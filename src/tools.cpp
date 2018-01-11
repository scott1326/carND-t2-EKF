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

  // Code from the lesson
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  
  // Check validity
  if(estimations.size() != ground_truth.size()
     || estimations.size() == 0) {
    cout << "Invalid estimation or ground_truth data" << endl;

    return rmse;
  }

  // accumulate squared residuals
  for(unsigned int i=0; i < estimations.size(); ++i){
    
    VectorXd residual = estimations[i] - ground_truth[i];

    residual = residual.array()*residual.array();
    rmse += residual;
  }

  // calculate the mean
  rmse = rmse/estimations.size();

  // calculate the square root
  rmse = rmse.array().sqrt();

  // return result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  
  // code from lesson
  MatrixXd Hj(3,4);
  
  // state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // pre-compute some terms
  float c1 = px*px + py*py;
  float c2 = sqrt(c1);
  float c3 = c1*c2;

  // check for divide/0
  if(fabs(c1) < 0.0001){
    cout << "CalculateJacobian() - Error - divide by 0" << endl;
    c1 = 0.0001;
  }
  if (fabs(c2) < 0.0001) {
	  cout << "CalculateJacobian() - Error - divide by 0" << endl;
	  c2 = 0.0001;
  }
  if (fabs(c3) < 0.0001) {
	  cout << "CalculateJacobian() - Error - divide by 0" << endl;
	  c3 = 0.0001;
  }


  // compute
  Hj << (px/c2), (py/c2), 0, 0,
        -(py/c1), (px/c1), 0, 0,
        py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

  return Hj;
}
