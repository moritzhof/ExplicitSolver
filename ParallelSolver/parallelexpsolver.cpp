//##############################################################################
//##########  Explicit Solver:Periodic and Non-Periodic Interval ###############
//######################### Moritz Travis Hof ##################################
//########## Parallel Implentation: First and Second Derivative ################
//##############################################################################
// Please note: The openMPI is writting with the C++ Binding. They have since dis
// continued this with MPI 3 and greater. Will still work with mpic++ compiler

#include<iostream>
#include<fstream>
#include<vector>
#include<math.h>
#include<iterator>
#include<mpi.h>


std::vector<double> Grid(int npts, double a, double b){
  std::vector<double> x;
  double dx =(b-a)/(npts-1.0);
  for(int i = 0; i<npts;++i)
    x.push_back(a+i*dx);
  return x;
}

double function(double x){
  return std::sin(x);
  //return std::pow(x,2);
}
double function_first_der(double x){
  return std::cos(x);
  //return std::pow(x,2);
}
double function_sec_der(double x){
  return -1.0*std::sin(x);
  //return 2*x;
}

void non_periodic_first_derivative(int npts, double dx, std::vector<double>& u, \
  std::vector<double>& u_x, int& totalnodes, int& rank){

    double interval = 1.0/(2.0*dx);
    double mpiTemp;
    MPI::Status status;

    if(rank == 0)
      u_x[0] = (-3.0*u[0] + 4.0*u[1] - u[2])*interval;

    if(rank == (totalnodes-1))
      u_x[npts-1] = (3.0*u[npts-1] - 4.0*u[npts-2] + u[npts-3])*interval;

    for(int i=1;i<npts-1;i++)
      u_x[i] = (u[i+1]-u[i-1])*interval;

    if(rank == 0){

      mpiTemp = u[npts-1];
      MPI::COMM_WORLD.Sendrecv_replace(&mpiTemp,1,MPI_DOUBLE,1,1,1,1, status);
      u_x[npts-1] = (mpiTemp - u[npts-2])*interval;


    }
    else if(rank == (totalnodes-1)){
      mpiTemp = u[0];
      MPI::COMM_WORLD.Sendrecv_replace(&mpiTemp,1,MPI_DOUBLE,rank-1,1, rank-1,1, status);
      u_x[0] = (u[1]-mpiTemp)*interval;
    }
    else{

      mpiTemp = u[0];
      MPI::COMM_WORLD.Sendrecv_replace(&mpiTemp,1,MPI_DOUBLE,rank-1,1, rank-1,1, status);
      u_x[0] = (u[1]-mpiTemp)*interval;

      mpiTemp = u[npts-1];
      MPI::COMM_WORLD.Sendrecv_replace(&mpiTemp,1,MPI_DOUBLE,rank+1,1,rank+1,1, status);
      u_x[npts-1] = (mpiTemp-u[npts-2])*interval;
    }

    return;
  }


void non_periodic_second_derivative(int npts, double dx, std::vector<double>& u, std::vector<double>& u_xx,
         int totalnodes, int rank){
 int i;
 double interval = 1.0/(dx*dx);
 double mpiTemp;
 MPI::Status status;

 if(rank == 0)
   u_xx[0] = (2.0*u[0]-5.0*u[1]+4.0*u[2]-u[3])*interval;

 if(rank == (totalnodes-1))
   u_xx[npts-1] = (2.0*u[npts-1]-5.0*u[npts-2]+4.0*u[npts-3]-u[npts-4])*interval;

 for(i=1;i<npts-1;i++)
   u_xx[i] = (u[i+1]-2.0*u[i]+u[i-1])*interval;

 if(rank == 0){

   mpiTemp = u[npts-1];
   MPI::COMM_WORLD.Sendrecv_replace(&mpiTemp,1,MPI_DOUBLE,1,1,1,1, status);
   u_xx[npts-1] = (mpiTemp - 2.0*u[npts-1] + u[npts-2])*interval;

 }
 else if(rank == (totalnodes-1)){

   mpiTemp = u[0];
   MPI::COMM_WORLD.Sendrecv_replace(&mpiTemp,1,MPI_DOUBLE,rank-1,1,rank-1,1, status);
   u_xx[0] = (u[1] - 2.0*u[0] + mpiTemp)*interval;

 }
 else{

   mpiTemp = u[0];
   MPI::COMM_WORLD.Sendrecv_replace(&mpiTemp,1,MPI_DOUBLE,rank-1,1,
      rank-1,1,status);
   u_xx[0] = (u[1] -2.0*u[0] + mpiTemp)*interval;

   mpiTemp = u[npts-1];
   MPI::COMM_WORLD.Sendrecv_replace(&mpiTemp,1,MPI_DOUBLE,rank+1,1,
      rank+1,1, status);
   u_xx[npts-1] = (mpiTemp -2.0*u[npts-1] + u[npts-2])*interval;

 }

 return;
}

int main(int argc, char *argv[]){

  int totalnodes;
  int rank;

  int local_npts, global_npts;
  double dx;
  double global_a = 0.0, global_b = 2.0*M_PI, local_a, local_b;

  std::vector<double> u;
  std::vector<double> u_x;
  std::vector<double> u_xx;
  std::vector<double> solution;
  std::vector<double> x;



  x = Grid(global_npts, global_a, global_b);
  u_x.resize(local_npts);
  u_xx.resize(local_npts);



     MPI::Init(argc, argv);

     rank = MPI::COMM_WORLD.Get_rank();
     totalnodes = MPI::COMM_WORLD.Get_size();

     global_npts = std::pow(2.0, 11);

     local_npts  = global_npts/totalnodes;
     global_npts = local_npts*totalnodes;
     dx = (global_b-global_a)/global_npts;

     local_a = global_a + dx*local_npts*rank;

     x = Grid(local_npts, global_a, global_b);
     u.resize(local_npts);
     u_x.resize(local_npts);
     u_xx.resize(local_npts);
     solution.resize(local_npts);

     for(int j=0;j<local_npts;j++){
       u[j] = function(local_a + j*dx);
       solution[j] = function_first_der(local_a + j*dx);
     }

     non_periodic_first_derivative(local_npts, dx, u, u_x, totalnodes, rank);
     non_periodic_second_derivative(local_npts,dx, u, u_xx,totalnodes, rank);

     //double approximation;
     //MPI::COMM_WORLD.Reduce(&ux_error,&answer1,1,MPI_DOUBLE,MPI_SUM,0);

     if(rank == 0){
     for(std::pair<std::vector<double>::iterator, std::vector<double>::iterator> i(solution.begin(), u_x.begin()); i.first !=solution.end(); ++i.first, ++i.second){
         std::cout << *i.first <<" \t "<< *i.second <<std::endl;
      }

      //output for plotting in gpuplot
      std::ofstream output_file("./data_first.txt");
      for(std::pair<std::vector<double>::iterator, std::vector<double>::iterator> i(x.begin(), u_x.begin()); i.first !=x.end(); ++i.first, ++i.second){
          output_file << *i.first <<" \t "<< *i.second <<std::endl;
      }

      std::ofstream output_file1("./data_second.txt");
      for(std::pair<std::vector<double>::iterator, std::vector<double>::iterator> i(x.begin(), u_xx.begin()); i.first !=x.end(); ++i.first, ++i.second){
          output_file1 << *i.first <<" \t "<< *i.second <<std::endl;
      }
      // use the following command in your terminal: gnuplot
      // gnuplot> plot sin(x), "data_x.txt" xrange [0:2*pi] yrange [-1:1] with line
    }


     MPI::Finalize();

  return 0;
}
