#ifndef PARAM
#define PARAM

#include <deal.II/base/parameter_handler.h>
#include <fstream>
#include <iostream>
#include "Functions_SEM.h"

using namespace dealii;

/* List of parameters
********************
*/

struct code_parameters
{

public:
  //default constructor
  //No other constructor is allowed since this object will contain data read from file
  code_parameters();


 //Medium param
  std::vector<std::string> material_id;
  std::vector<std::string>velocities;
  std::map<unsigned int, double> material_vel;

//Source
  double fo;
  double source_rad;
  double source_pos_x;
  double source_pos_y;

// Receivers
  int n_receivers;
  double  y_pos;
  double  x_init;
  double  x_end;

// Discretization
  std::string input_mesh;
  int refinement;
  unsigned int add_refinement; //around source
  double total_time;
  double time_step;
  double theta;

// finite element
  int set_nq_points;
  int fe_degree;
};

class ParameterReader: public Subscriptor
{
public:

  ParameterReader(const char *parameter_file_);

  //
  code_parameters get_parameters();


private:

  void declare_parameters();
  void read_parameters();
  void pass_parameters ();

  const char* parameter_file;

  //default constructed objects. Note: default constructors should be defined for these class objects
  ParameterHandler  prm;
  code_parameters param_input;

};


#endif
