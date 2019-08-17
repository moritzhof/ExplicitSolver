//##############################################################################
//##########  Explicit Solver:Periodic and Non-Periodic Interval ###############
//######################### Moritz Travis Hof ##################################
//##############################################################################

#include<iostream>
#include<fstream>
#include<vector>
#include<math.h>
#include<iterator>


void print_vector(std::vector<double>& v){
  for(const auto x : v){
    std::cout<< x << " ";
  }
  std::cout<<std::endl;
}

double function(double& x){
  return std::sin(x);
  //return std::pow(x,2);
}
double function_sec_der(double& x){
  return -1.0*std::sin(x);
  //return 2*x;
}

std::vector<double> Grid(int npts, double a, double b){
  std::vector<double> x;
  double dx =(b-a)/(npts-1.0);
  for(int i = 0; i<npts;++i)
    x.push_back(a+i*dx);
  return x;
}

void periodic_first_derivative(int npts, double dx, std::vector<double>& u, std::vector<double>& u_x){

  double interval = 1.0/(2.0*dx);
  u[0] = (u[0]-u[npts-1])*interval;

  for(int i = 0; i<npts; ++i){
    u_x[i] = (u[i+1]-u[i-1])*interval;
  }

  u[npts-1] = (u[0]-u[npts-2])*interval;

  return;
}

void periodic_second_derivative(int npts, double dx, std::vector<double>& u, std::vector<double>& u_xx){

  double interval = 1/std::pow(dx,2);
  u_xx[0] = (u[1]-2.0*u[0]+u[npts-1])*interval;

  for(int i(1); i<npts-1; ++i){
    u_xx[i] = (u[i+1]-2.0*u[i]+u[i-1])*interval;
  }

  u_xx[npts-1] = (u[0]-2.0*u[npts-2]+u[npts-2])*interval;

  return;
}

void non_periodic_first_derivative(int npts, double dx, std::vector<double>& u, std::vector<double>& u_x){

  double interval = 1.0/(2.0*dx);
  u[0] = (-3.0*u[0]+4.0*u[1]-u[2])*interval;

  for(int i = 0; i<npts; ++i){
    u_x[i] = (u[i+1]-u[i-1])*interval;
  }

  u[npts-1] = (-3.0*u[npts-1]+4.0*u[npts-2]-u[npts-3])*interval;

  return;
}

void non_periodic_second_derivative(int& npts, double& dx, std::vector<double>& u, std::vector<double>& u_xx){


  double interval = 1/std::pow(dx,2);

  u_xx[0] = (2.0*u[0]-5.0*u[1]+4.0*u[2]-u[3])*interval;

  for(int i=1; i<npts-1; ++i){
    u_xx[i] = (u[i+1]-2.0*u[i]+u[i-1])*interval;
  }

  u_xx[npts-1] = (2.0*u[npts-1]-5.0*u[npts-2]+4.0*u[npts-3]-u[npts-4])*interval;

}


int main(int argc, char const *argv[]){
  int choice;
  std::cout << "Please enter if you would like period or non-period on execution" << std::endl;
  std::cout << "Please enter 0 for periodic or 1 for non-periodic" << std::endl;
  std::cin >> choice;
  int npts;
  double dx;
  double a = 0.0, b = 2.0*M_PI;
  std::vector<double> u;
  std::vector<double> u_x;
  std::vector<double> u_xx, u_xx_error;
  std::vector<double> solution;
  std::vector<double> x;

    npts = std::pow(2.0, 11);
    dx = (b-a)/(npts-1);

     x = Grid(npts, a, b);
     u_x.resize(npts);
     u_xx.resize(npts);
     u_xx_error.resize(npts);

    for(std::vector<double>::iterator it = x.begin(); it != x.end(); ++it){
      u.push_back(function(*it));
      solution.push_back(function_sec_der(*it));
    }

    if(choice == 1){
      non_periodic_first_derivative(npts, dx, u, u_x);
      non_periodic_second_derivative(npts, dx, u, u_xx);
    }
    else{
      periodic_first_derivative(npts, dx, u, u_x);
      periodic_second_derivative(npts, dx, u, u_xx);
    }
    for(int j=0; j<npts;++j){
      u_xx_error[j] = std::sqrt(dx*std::pow(u_xx[j]-solution[j], 2));
    }

 if(choice == 1){
   std::cout << "Computations for Non-Periodic Interval Approximations:" << std::endl;
 }
 else{
    std::cout << "Computations for Periodic Interval Approximations:" << std::endl;
 }
 for(std::pair<std::vector<double>::iterator, std::vector<double>::iterator> i(solution.begin(), u_xx.begin()); i.first !=solution.end(); ++i.first, ++i.second){
     std::cout << *i.first <<" \t "<< *i.second <<std::endl;
  }
  std::cout<<std::endl;
  std::cout<< " Error Calculations " << std::endl;
  for(std::vector<double>::iterator it = u_xx_error.begin(); it != u_xx_error.end(); ++it){
    std::cout << *it << " ";
  }
  std::cout<<std::endl;

  //output for plotting in gpuplot
  std::ofstream output_file("./data.txt");
  for(std::pair<std::vector<double>::iterator, std::vector<double>::iterator> i(x.begin(), u_xx.begin()); i.first !=x.end(); ++i.first, ++i.second){
      output_file << *i.first <<" \t "<< *i.second <<std::endl;
  }
  // use the following command in your terminal: gnuplot
  // gnuplot> plot sin(x), "data.txt" xrange [0:2*pi] yrange [-1:1] with line


//########## Print all vectors ... if you would like to ########################
//########################### Just uncomment ###################################
  // print_vector(u);
  // print_vector(u_x);
  // print_vector(u_xx);
  return 0;
}
