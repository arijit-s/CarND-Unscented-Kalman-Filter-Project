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

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	if(estimations.size() == 0 || estimations.size() != ground_truth.size()){
		cout<<"Estimator Vector is size zero"<<endl;
		return rmse;
	}
	//  * the estimation vector size should equal ground truth vector size

	//accumulate squared residuals
	for(int i=0; i < estimations.size(); ++i){
        // ... your code here
		VectorXd c = estimations[i] - ground_truth[i];
		VectorXd c_square = c.array()*c.array();
		rmse = rmse + c_square;
	}

	//calculate the mean
	rmse = rmse/estimations.size();


	//calculate the squared root
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;
}