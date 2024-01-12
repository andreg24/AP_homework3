#include "odesolver.hpp"

ODESolver::ODESolver(function<VectorXd(double, VectorXd)> fun, double t_start, double t_end, VectorXd y0) : fun(fun), t_start(t_start), t_end(t_end), y0(y0) {};

  ODESolution ODESolver::RK4(const unsigned int n) {
    const double h = (t_end - t_start)/n; //step height. It is fixed
    VectorXd t(n+1); //vector containing t values each step
    t[0] = t_start; 
    vector<VectorXd> Y; //vector containing the values of the solution in t[i] each step.
    VectorXd y = y0; 
    Y.push_back(y0);

    for (unsigned int j=0; j<n; j++) {
      //Runge Kutta increments 
      VectorXd k1 = fun(t[j], y); 
      VectorXd k2 = fun(t[j]+h/2.0, y+0.5*k1*h);
      VectorXd k3 = fun(t[j]+h/2.0, y+0.5*k2*h);
      VectorXd k4 = fun(t[j]+h, y+k3*h);
      // t update
      t[j+1]=t[j]+h;
      // y update using RK4 method 
      y+=+h/6.0*(k1+2*k2+2*k3+k4);
      Y.push_back(y);
      }
    ODESolution res;
    res.t = t;
    res.Y = Y;
    return res;
  }

  void ODESolver::RK4_csv(const unsigned int n,const string filename) {
       const double h = (t_end - t_start)/n; //fixed step height
       double t = t_start;       
       VectorXd y = y0;
       std::ofstream myfile;
      myfile.open(filename); //Opening the file in which we are going to write the results
      if (myfile.is_open()) {

        //header of the csv
        myfile << "t" <<" , ";
        for (unsigned int i=0; i<y.size()-1; i++) {
          myfile<< "y" << i<< " , ";
        };
        myfile << "y"<< y.size()-1<< "\n";
        myfile << t_start << " , "<< y0 << "\n";

        for (unsigned int j=0; j<n; j++) {
           VectorXd k1 = fun(t, y);
           VectorXd k2 = fun(t+h/2.0, y+0.5*k1*h);
           VectorXd k3 = fun(t+h/2.0, y+0.5*k2*h);
           VectorXd k4 = fun(t+h, y+k3*h);
           t+=h;

           y+=+h/6.0*(k1+2*k2+2*k3+k4);
           myfile << t << ",";
           for (unsigned int i=0; i<y.size()-1; i++) {
            myfile << y(i) << ",";
           }
           myfile << y(y.size()-1);
           myfile << "\n";

        }

        myfile.close(); //Closing the file
        std::cout << "Results are now written in resultsRK4.csv" << std::endl;
      } 
      else {
      std::cout << "Unable to open the file 'resultsRK4.csv'" << std::endl;
      }
  }

  ODESolution ODESolver::midpoint(const unsigned int n) {

    const double h = (t_end - t_start) / n; //fixed step height
    VectorXd t(n+1);
    t[0] = t_start;
    vector<VectorXd> Y;
    VectorXd y = y0;
    Y.push_back(y0);

    for (unsigned int j = 0; j < n; j++) {
      VectorXd k1 = fun(t[j], y);
      VectorXd k2 = fun(t[j] + h / 2.0, y + 0.5 * h * k1);
      t[j + 1]=t[j] + h;
      y += h * k2;
      Y.push_back(y);
    }
    ODESolution res;
    res.t = t;
    res.Y = Y;
    return res;
   }

  void ODESolver::midpoint_csv(const unsigned int n,const string filename) {
    const double h = (t_end - t_start) / n;
    double t = t_start;
    VectorXd y = y0;

    std::ofstream myfile;
    myfile.open(filename);
    if (myfile.is_open()) {
      myfile << "t" <<" , ";
      for (unsigned int i=0; i<y.size()-1; i++) {
         myfile<< "y" << i<< " , ";
      };
    myfile << "y"<< y.size()-1<< "\n";
    myfile << t_start << " , "<< y0 << "\n";

    for (unsigned int j = 0; j < n; j++) {
      VectorXd k1 = fun(t, y);
      VectorXd k2 = fun(t + h / 2.0, y + 0.5 * h * k1);
      t+=h;

      y += h * k2;
      myfile << t << ",";
           for (unsigned int i=0; i<y.size()-1; i++) {
            myfile << y(i) << ",";
           }
           myfile << y(y.size()-1);
           myfile << "\n";
           //std::cout<<"ciao"<<y<<endl;
    }

   myfile.close();
        std::cout << "Results are now written in results_midpoint.csv" << std::endl;
    } 
    else {
      std::cout << "Unable to open the file 'results_midpoint.csv'" << std::endl;
    }

}

  ODESolution ODESolver::euler(const unsigned int n) {
    
    const double h = (t_end - t_start) / n;
    VectorXd t(n+1);
    t[0] = t_start;
    vector<VectorXd> Y;
    VectorXd y = y0;
    Y.push_back(y0);

    for (unsigned int j = 0; j < n; j++) {
      VectorXd k1 = fun(t[j], y);
      t[j + 1] = t[j] + h;
      y += h * k1;
      Y.push_back(y);
    }

    ODESolution res;
    res.t = t;
    res.Y = Y;
    return res;
  }

  void ODESolver::euler_csv(unsigned int n,const string filename) {
    const double h = (t_end - t_start) / n;
    double t = t_start;
    VectorXd y = y0;

    std::ofstream myfile;
    myfile.open(filename);
    if (myfile.is_open()) {
      myfile << "t" << " , ";
      for (unsigned int i = 0; i < y.size() - 1; i++) {
        myfile << "y" << i << " , ";
      };
      myfile << "y" << y.size() - 1 << "\n";
      myfile << t_start << " , "<< y0 << "\n";


      for (unsigned int j = 0; j < n; j++) {
        VectorXd k1 = fun(t, y);
        t += h;
        y += h * k1;
        myfile << t << ",";
        for (unsigned int i = 0; i < y.size() - 1; i++) {
          myfile << y(i) << ",";
        }
        myfile << y(y.size() - 1);
        myfile << "\n";
      }

      myfile.close();
      std::cout << "Results are now written in 'results_euler.csv'" << std::endl;
    } 
    else {
      std::cout << "Unable to open the file 'results_euler.csv'" << std::endl;
    }
  }

  double ODESolver::accuracy(ODESolution& res, const function<VectorXd(double)> analitic) {
        double max_error = 0.0;
        const vector<VectorXd>& Y = res.Y;
        const VectorXd& t = res.t;
        for (unsigned int i = 0; i < t.size(); ++i) {
            VectorXd y_analitic = analitic(t[i]); // analitic solution
            VectorXd diff = Y[i] - y_analitic; // Difference between analitic and numerical solutions
            double error = diff.norm(); // Error in each pass
            max_error = std::max(max_error, error); // Update of the max error
        }

        return max_error;
    }

  double ODESolver::efficiency(const string method, const unsigned int n) {
    if (method == "RK4") { // efficiency for RK4 method
      std::chrono::duration<double> sum; 
      for (unsigned int i=0; i<5; i++) { // test RK4 5 times
        auto start = std::chrono::high_resolution_clock::now();
        RK4(n); 
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time = end - start; // time for each computation of RK4
        sum+= time; 
      }
      chrono::duration<double> a = sum/5; //mean of RK4 computations
      return a.count();
    }
    else if (method == "midpoint") { // same procedure for midpoint
      std::chrono::duration<double> sum;
      for (unsigned int i=0; i<5; i++) {
        auto start = std::chrono::high_resolution_clock::now();
        midpoint(n); // 
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time = end - start;
        sum+= time;
      }
      chrono::duration<double> a = sum/5;
      return a.count();
      }
    else if (method == "euler") { // same procedure for euler
      std::chrono::duration<double> sum;
      for (unsigned int i=0; i<5; i++) {
        auto start = std::chrono::high_resolution_clock::now();
        euler(n); 
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time = end - start;
        sum+= time;
      }
      chrono::duration<double> a = sum/5;
      return a.count();
      }
      else { // 
          throw invalid_argument("Invalid input for efficency. You have to choose between RK4, midpoint or euler!");
      }
    };

  pair<double, double> ODESolver::stability(const string method, ODESolution& res,const double p){
    const VectorXd pert = VectorXd::Constant(res.Y[0].size(), p);
    const VectorXd y0_pert=res.Y[0]+pert; // initial value perturbed

    const unsigned int n=res.Y.size(); // number of subintervals
    // Compute the numerical solver with perturbed initial value
    ODESolver perturbed(fun,t_start,t_end,y0_pert); 
    ODESolution res_pert;  // Computer the numerical solution with the choosen initial value

      if (method == "RK4") {
          res_pert = perturbed.RK4(n);
      } else if (method == "midpoint") {
          res_pert = perturbed.midpoint(n);
      } else if (method == "euler") {
          res_pert = perturbed.euler(n);
      } else {
          throw invalid_argument("Invalid input for efficency. You have to choose between RK4, midpoint or euler! ");
      }

    double max_error = 0.0;
    const vector<VectorXd>& Y = res.Y; // numerical solution of the original equation
    const vector<VectorXd>& Y_pert=res_pert.Y; // numerical solution of the perturbed equation
        for (unsigned int i = 0; i < n; ++i) {
            VectorXd diff = Y[i] - Y_pert[i]; // difference between the two solutions of the two equations
            double error = diff.norm(); // norm of the error
            max_error = std::max(max_error, error); // update max error
        }
    pair<double, double> stab; 
    stab.first = p; // perturbation of initial value
    stab.second = max_error; // perturbation of the solution
    return stab;
    };

  double ODESolver::convergence(ODESolution& res, const function<VectorXd(double)> analytic_solution) {
    const unsigned int n = res.Y.size();
    const double h = (t_end - t_start)/n;
    VectorXd h_values(4); // we compute the solver with 4 differnt step height;
    h_values << h, h/2, h/4, h/8; // halved steps

    VectorXd errors(4); // vector containing max error for each h
    for (int i = 0; i < h_values.size(); ++i) {
      auto numerical_solution = euler((t_end - t_start) / h_values(i));

      double max_error = accuracy(numerical_solution, analytic_solution); // max error
      errors(i) = max_error;
    }

    VectorXd convergence_order(3); 
    for (int i = 1; i < errors.size(); ++i) {
      double order = log(errors(i - 1) / errors(i)) / log(2.0); // estimation of convergence error between two computation with two different h
      convergence_order(i - 1) = order;

    }
    // Average convergence error
    double avg_order = convergence_order.mean();
    return avg_order;
    }

    
