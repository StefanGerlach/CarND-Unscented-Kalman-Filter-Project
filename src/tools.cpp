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
      * Calculate the RMSE here.
    */
    if (estimations.size() != ground_truth.size() || estimations.size() <= 0) {
      std::cout << "Error computing the Root Mean Squared Error." << std::endl;
    }

    // probe the dimensions. we can access at postion 0 since we checked the length
    // in the statement above
    unsigned int dimensions = estimations[0].rows();
    double tmp_difference = 0.0;

    VectorXd rmse(dimensions);
    rmse.setConstant(0.0);

    // iterate all elements
    for(unsigned int i = 0; i < estimations.size(); i++) {

      // iterate all dimensions
      for(unsigned int d = 0; d < dimensions; d++) {
        tmp_difference = (estimations[i](d) - ground_truth[i](d));
        rmse(d) +=  tmp_difference * tmp_difference;
      }
    }

    // compute mean and square root
    for(unsigned int d = 0; d < dimensions; d++) {
      rmse(d) /= static_cast<double>(estimations.size());
      rmse(d) = std::sqrt(static_cast<double>(rmse(d)));
    }

    return rmse;
}


VectorXd Tools::ConvertPolarToCartesian(const VectorXd &polar_measurement) {

  int expected_dims = 3;
  int destination_dims = 4;

  VectorXd cartesian_measurement(destination_dims);
  cartesian_measurement.setConstant(0.0);

  if(polar_measurement.rows() != expected_dims) {
    throw("ConvertPolarToCartesian() input does not match expected dimensions!");
  }

  float rho = polar_measurement(0); // Distance
  float phi = polar_measurement(1); // Bearing
  float rho_dot = polar_measurement(2); // Range Rate

  cartesian_measurement(0) = std::cos(phi) * rho;
  cartesian_measurement(1) = std::sin(phi) * rho;

  return cartesian_measurement;
}
