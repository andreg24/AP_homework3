#ifndef ODESOLVER_HPP_
#define ODESOLVER_HPP_

#include <iostream>
#include <functional>
#include <cmath>
#include <fstream>
#include <vector>
#include "Eigen/Dense"
#include <chrono>

using namespace Eigen;
using namespace std;

struct ODESolution {
  VectorXd t;
  vector<VectorXd> Y;
};

class ODESolver {
private:
  const function<VectorXd(double,VectorXd)> fun; //Function to implement outside the class
  const double t_start; //Initial value of the interval of the domain of the function
  const double t_end; //Last value of the interval of the domain of the function
  const VectorXd y0; //Initial value of the Cauchy problem

public:
    //Constructor
    ODESolver(function<VectorXd(double, VectorXd)> fun, double t_start, double t_end, VectorXd y0);

    //Runge Kutta 4 method for solving ODEs
    ODESolution RK4(const unsigned int n);

    //Variation of RK4 that saves the results into a csv
    void RK4_csv(const unsigned int n, const string filename);

    //Midpoint method for solving ODEs
    ODESolution midpoint(const unsigned int n);

    //Variation of midpoint method that saves the results into a csv
    void midpoint_csv(const unsigned int n, const string filename);

    // Forward Euler methord for solving ODEs
    ODESolution euler(const unsigned int n);

    // Variation of euler method that saves the results into a csv
    void euler_csv(const unsigned int n,const string filename);

    // Accuracy
    double accuracy(ODESolution& res,const function<VectorXd(double)> analitic);

    //Efficiency
    double efficiency(const string method, const unsigned int n);

    //Stability
    pair<double, double> stability(const string method, ODESolution& res, const double p);

    //Convergence
    double convergence(ODESolution& res,const function<VectorXd(double)> analytic_solution);
 
};


#endif
