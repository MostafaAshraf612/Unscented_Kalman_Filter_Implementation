#include "tools.h"
#include <iostream>
#include <cmath>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

/**
 * Calculate the Root Mean Square Error (RMSE)
 */
VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    // Initialize RMSE vector
    VectorXd rmse(4);
    rmse << 0, 0, 0, 0;

    // Check validity of inputs
    if (estimations.empty() || estimations.size() != ground_truth.size()) {
        std::cout << "Error: Estimation and ground truth vectors are not of the same size or are empty." << std::endl;
        return rmse;
    }

    // Accumulate squared residuals
    for (size_t i = 0; i < estimations.size(); i++) {
        VectorXd residual = estimations[i] - ground_truth[i];
        residual = residual.array() * residual.array();
        rmse += residual;
    }

    // Calculate the mean
    rmse = rmse / estimations.size();

    // Calculate the squared root
    rmse = rmse.array().sqrt();

    return rmse;
}



