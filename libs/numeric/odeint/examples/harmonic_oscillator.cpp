/* Boost numeric/odeint/examples/harm_osc.cpp
 
 Copyright 2009 A. Constantine Ashford

 Demonstrates the use of odeint to estimate the (known) solution to a simple,
 undriven harmonic  oscillator:

 y'' = -y - gamma*y'

 This example uses the Euler solver as well as he fourth-order Runge Kutta. 
 The output of each is printed in the output. The accuracy can be compared, or
 the program can be profiled for a speed comparison.
 
 The analytical solution to this ODE is (with initial displacement y and no 
 initial velocity y'):

 Y = exp( -.5*gamma*t ) * cos( sqrt(1-gamma^2)*t )

 This example calculates the analytical solution using the standard math
 library and outputs the values at the same time points for comparison and
 plotting.

 On a unix system, a plot of the output could be created easily using gnuplot
 and the following commands:
 
 $ plot "harm_osc.out" using 1:2 title 'Euler' with dots,\
 $ "harm_osc.out" using 1:5 title 'Runge Kutta 4' with dots,\
 $ "harm_osc.out" using 1:8 title 'Exact' with dots

 The step size, number of steps to integrate, output file name, and gamma can
 be specified on the command line IFF you have the boost library files installed
 *which are at the same version as this example is compiled against*. You can
 enable the command line options by compiling with the preprocessor variable
 BOOST_PROGRAM_OPTIONS_INSTALLED defined.


 Compile with:
 g++ -Wall -O3 -march=k8 -I$BOOST_DIR -I../../../../ harmonic_oscillator.cpp

 Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>

#include <boost/numeric/odeint.hpp>

#ifdef BOOST_PROGRAM_OPTIONS_INSTALLED
#include <boost/program_options.hpp>
namespace opts=boost::program_options;
#endif

namespace odeint=boost::numeric::odeint;

// Namespace for variables specific to this example
namespace harm_osc
{
  double gamma = .15;
  double dt = 0.01; 
  size_t olen = 10000;
  std::ofstream fout;

  // The type of container used to hold the state vector
  typedef std::vector<double> state_type;

  /*
    void harmonic_oscillator(state_type &x, state_type &dxdt, double t)

    Calculates the derivatives of the state vector for the harmonic oscillator
    system. The system equation is described in the opening comments.

    The state vector is:
     x[0] = Y (the output)
     x[1] = Y' (dY/dt)

    Arguments:
     state_type &x: reference to the current state vector [ Y Y' ]
     state_type &dxdt: reference to the vector of calculated derivitives
      [ Y' Y'']
     double t: the current integration time
  */
  void harmonic_oscillator(const state_type &x, state_type &dxdt, const double t)
  {
    dxdt[0] = x[1];
    dxdt[1] = -x[0] - gamma*x[1];
  }

  /*
    Optionally included function to parse the command line. See detailed
    comments and implementation at the bottom of this file
  */
  int parse_command_line(int ac, char ** av);
}

// Main
int main(int argc, char **argv)
{
  using namespace std;

  // Save the original destination of cout in case we redirect to a file
  // as a result of --filename command line option
  streambuf * cout_buf = cout.rdbuf();

  // Parse the command line options iff boost::program_options is installed
  #ifdef BOOST_PROGRAM_OPTIONS_INSTALLED
  int parse_ret = harm_osc::parse_command_line(argc, argv);
  if(parse_ret)
      return parse_ret;
  #endif

  // Declare state vectors for the Euler system, and the Runge Kutta system
  harm_osc::state_type x_euler(2);
  harm_osc::state_type x_rk4(2);

  // Set identical initial conditions for the Euler solver and the Runge Kutta
  x_euler[0] = 1.0;
  x_euler[1] = 0.0;
  x_rk4[0] = 1.0;
  x_rk4[1] = 0.0;

  // Initialize the solvers
  odeint::stepper_euler<harm_osc::state_type> euler;
  odeint::stepper_rk4<harm_osc::state_type> rk4;

  euler.adjust_size( x_euler );
  rk4.adjust_size( x_rk4 );

  
  // Write the output header. The '#' symbol creates a comment in gnuplot
  cout << "# Output from the simple harmonic oscillator example." << endl;
  cout << "#" << endl;
  cout << "# Columns:" << endl;
  cout << "#" << endl;
  cout << "# 	1: time (s)" << endl;
  cout << "#	2: Euler output" << endl;
  cout << "#	3: Euler output dY/dt" << endl;
  cout << "# 	4: Euler output error" << endl;
  cout << "#	5: 4th Order Runge Kutta output" << endl;
  cout << "#	6: 4th Order Runge Kutta dY/dt" << endl;
  cout << "# 	7: 4th Order Runge Kutta output error" << endl;
  cout << "# 	8: Analytical (calculated) solution" <<endl;
  cout << "#" << endl;
  cout << "# A. Constantine Ashford" << endl;
  cout << "#" << endl;

// Initialize the time variable, and the analytical solution variable
  double t = 0.0;
  double reference;
  double e_error = 0, r_error = 0;
  double e_error_max = 0, r_error_max = 0;
  double e_error_rms = 0, r_error_rms = 0;
  int data_column_width = 15;

  // Integrate and write the approximate and analytical results to a table
  for(size_t oi=0; oi < harm_osc::olen; ++oi, t += harm_osc::dt)
    {
      // Calculate the analytical solution for this time
      reference = exp(-.5*harm_osc::gamma*t) *
          cos(sqrt(1-pow(harm_osc::gamma,2))*t);

      // Calculate absolute errors for the data table
      e_error = abs((reference - x_euler[0])/reference);
      r_error = abs((reference - x_rk4[0])/reference);
         
      // Write all the results to a row in the table
      cout << setw(5) << t
	   << setw(data_column_width) << x_euler[0]
	   << setw(data_column_width) << x_euler[1]
	   << setw(data_column_width) << e_error
	   << setw(data_column_width) << x_rk4[0]
	   << setw(data_column_width) << x_rk4[1]
	   << setw(data_column_width) << r_error
	   << setw(data_column_width) << reference
	   << endl;

      // Aggregate statistics
      if( (e_error = abs(e_error)) > e_error_max)
          e_error_max = e_error;
      
      if((r_error = abs(r_error)) > r_error_max)
          r_error_max = r_error;

      e_error_rms += pow(e_error, 2);
      r_error_rms += pow(r_error, 2);
       
      // Advance the solvers
      euler.do_step(harm_osc::harmonic_oscillator,
                      x_euler, t, harm_osc::dt);

      rk4.do_step(harm_osc::harmonic_oscillator,
                    x_rk4, t, harm_osc::dt);
      
    }

  e_error_rms = sqrt(e_error_rms/harm_osc::olen);
  r_error_rms = sqrt(r_error_rms/harm_osc::olen);
  
  // Restore cout to direct to standard out in case we redirected to a file
  cout.rdbuf(cout_buf);
  cout << "#Integration complete." << endl;
  cout << "#Errors for the Euler solver: MAX: " << e_error_max << " RMS: "
       << e_error_rms << endl;
  cout << "#Errors for the Runge Kutta solver: MAX: " << r_error_max << " RMS: "
       << r_error_rms << endl;
  

  // Success
  return 0;
}

#ifdef BOOST_PROGRAM_OPTIONS_INSTALLED
/*
parse_command_line(int ac, char ** av)

Uses boost::program_options to parse the command line.
Allows the example to use runtime options such as:
--step=.001 or 
--step .001

Arguments:
 int ac: the number of command line arguments. Usually just pass in argc
 from main.

 char ** av: an array of c-style strings holding the command line arguments.
 Usually just pass in argv from main.

Return value: 1 if the user chose the --help option, and does *not* wish to
run an integration. Otherwise 0.

*/
int harm_osc::parse_command_line(int ac, char ** av)
{
  using namespace std;

  // Declare supported options
  opts::options_description opts_desc("Allowed options");
  opts_desc.add_options()
    ("help", "Display this message.")
    ("file", opts::value<string>(), "Output file for data.(Default: stdout)")
    ("step", opts::value<double>(), "Set step size (in seconds)")
    ("duration", opts::value<size_t>(), "Set integration duration (in steps)")
    ("gamma", opts::value<double>(), "Set gamma equation parameter");

  // Let Boost Program Options parse the command line
  opts::variables_map opts_map;
  opts::store(opts::parse_command_line(ac, av, opts_desc), opts_map);
  opts::notify(opts_map);

  // Print help message then exit
  if(opts_map.count("help"))
    {
       cout << opts_desc << endl;
       return 1;
    }

  // Alter step size in seconds
  if(opts_map.count("step"))
    {
      harm_osc::dt = opts_map["step"].as<double>();
      cout << "Step size set to " << harm_osc::dt << " seconds" << endl;
    }

  // Change integration duration (in steps)
  if(opts_map.count("duration"))
    {
      harm_osc::olen = opts_map["duration"].as<size_t>();
      cout << "Integration duration set to " << harm_osc::olen << " steps."
	   << endl;
    }

  // Change damping constant
  if(opts_map.count("gamma"))
    {
      harm_osc::gamma = opts_map["gamma"].as<double>();
      cout << "Parameter GAMMA set to " << harm_osc::gamma << "." << endl;
    }

  // Redirect output from standard out to file 
  if(opts_map.count("file"))
    {
      cout << "Writing output to file " << opts_map["file"].as<string>()
	   << endl;
      fout.open(opts_map["file"].as<string>().c_str());
      cout.rdbuf(fout.rdbuf());
    }

  return 0;
}
#endif
