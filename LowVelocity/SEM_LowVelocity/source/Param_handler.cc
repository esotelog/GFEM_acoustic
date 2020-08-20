
#include <string>
#include "Functions_SEM.h"
#include "Param_handler.h"

using namespace dealii;

/* Code_parameters default constructor
 *****************************
 */
code_parameters::code_parameters():
  // Medium:


  //source
  fo(40.0),
  source_rad(3.125),
  source_pos_x(400),
  source_pos_y(-50),
  //source_pos(0.25,1.875),

  //receivers
  n_receivers (100),
  y_pos (-10),
  x_init (50),
  x_end (750),

  //Discretization: grid & time
  input_mesh("h_mesh.msh"),
  refinement (4),
  add_refinement(0),
  total_time (0.5),
  time_step (0.00165),
  theta (0.51),

   // finite element
  set_nq_points (3),
  fe_degree(1)

{}



/*ParameterReader class
 * ********************



constructor
 *************
 */
ParameterReader::ParameterReader(const char * parameter_file_):
  parameter_file(parameter_file_)
{}

/*Declare Parameters
 * *****************
 */
void ParameterReader::declare_parameters()
{
  //Medium
  prm.enter_subsection("Medium");
  {
      prm.declare_entry("material id",
                        "1",
                        Patterns::Anything(),
                        "material id");

     prm.declare_entry ("velocities",
                        "1000",
                         Patterns::Anything(),
                         "material velocities");


  }
  prm.leave_subsection();



  //Source
  prm. enter_subsection("Source");
  {
    prm.declare_entry("frequency",
              std::to_string(param_input.fo),
              Patterns::Double(0),
              "central frequency of source (Ricker wave), (double)"
            );

    prm.declare_entry ("source radius",
              std::to_string(param_input.source_rad),
              Patterns::Double(0),
              "spatial extent of source, (double)"
              );

    prm.declare_entry ("source_pos_x",
               std::to_string(param_input.source_pos_x),
               Patterns::Double(0),
               "x position of source,  (double)"
               );

    prm.declare_entry ("source_pos_y",
               std::to_string(param_input.source_pos_y),
               Patterns::Double(-10000),
               "y position of source,  (double)"
               );

  }
  prm.leave_subsection();

  //Receivers: line of horizontal receivers

  prm. enter_subsection("Receivers");
  {
    prm.declare_entry ("n_receivers",
               std::to_string(param_input.n_receivers),
               Patterns::Integer(0),
               "number of receivers,  (integer)"
                );

    prm.declare_entry ("y_pos",
              std::to_string(param_input.y_pos),
              Patterns::Double(-10000),
              "y position of receivers, (double)"
              );

    prm.declare_entry ("x_init",
              std::to_string(param_input.x_init),
              Patterns::Double(0),
              "x position of first receiver, (double)"
              );

    prm.declare_entry ("x_end",
              std::to_string(param_input.x_end),
              Patterns::Double(0),
              "x position of last receiver,  (double)"
              );

  }
  prm.leave_subsection();

  prm. enter_subsection("Discretization");
  {
    prm.declare_entry ("mesh file",
                       param_input.input_mesh,
                       Patterns::Anything(),
                       "name of mesh file, include .msh extension"
                        );

    prm.declare_entry ("refinement level",
                        std::to_string(param_input.refinement),
                        Patterns::Integer(0),
                        "mesh refinement,  (integer)"
                        );

    prm.declare_entry  ("add refinement",
                         std::to_string(param_input.add_refinement),
                         Patterns::Integer(0),
                         "add mesh refinement,  (integer)"
                         );

    prm.declare_entry ("simulation_total_time",
                        std::to_string(param_input.total_time),
                        Patterns::Double(0),
                        "total time of simulation, (double)"
                       );

    prm.declare_entry ("time_step",
                        std::to_string(param_input.time_step),
                        Patterns::Double(0),
                        "simulation time_step, (double)"
                       );

    prm.declare_entry ("theta",
                        std::to_string (param_input.theta),
                        Patterns::Double(0,1),
                        "time discretization type: 0 (explicit Euler), 0 (implicit Euler) 0.5 (Crank-Nicolson) (double)"
                       );

  }
  prm.leave_subsection();

  prm. enter_subsection("Finite Element");
  {
    prm.declare_entry ("n_quadrature",
                       std::to_string(param_input.set_nq_points),
                       Patterns::Integer(0),
                       "Number of quadrature points"
                       );

    prm.declare_entry ("fe_degree",
                       std::to_string(param_input.fe_degree),
                       Patterns::Integer(0),
                       "degree of finite element"
                       );

  }
  prm.leave_subsection();

}

/*Read Parameters from file
 * ***********************
 */

void ParameterReader::read_parameters ()
{
  prm.parse_input (parameter_file);
}

/*Pass parameters to structure
 * ***********************
 */

void ParameterReader::pass_parameters()
{
   prm.enter_subsection("Medium");
     param_input.material_id=Utilities::split_string_list( prm.get("material id"), ';');
     param_input.velocities=Utilities::split_string_list( prm.get("velocities"), ';');
   prm.leave_subsection();

   //creating the map material_id->velocity
   for (unsigned int i=0; i< param_input.material_id.size(); i++)
     param_input.material_vel[std::stoul( param_input.material_id[i]) ]
         = std::stold( param_input.velocities[i]);

   prm.enter_subsection("Source");
      param_input.fo = prm.get_double ("frequency");
      param_input.source_rad = prm.get_double ("source radius");
      param_input.source_pos_x = prm.get_double ("source_pos_x");
      param_input.source_pos_y = prm.get_double ("source_pos_y");
   prm.leave_subsection();

   prm.enter_subsection("Receivers");
      param_input.n_receivers= prm.get_integer("n_receivers");
      param_input.y_pos = prm.get_double("y_pos");
      param_input.x_init = prm.get_double("x_init");
      param_input.x_end = prm.get_double("x_end");
   prm.leave_subsection();

   prm.enter_subsection("Discretization");
      param_input.input_mesh = prm.get("mesh file") ;
      param_input.refinement = prm.get_integer("refinement level");
      param_input.add_refinement =prm.get_integer("add refinement");
      param_input.total_time = prm.get_double("simulation_total_time");
      param_input.time_step = prm.get_double ("time_step");
      param_input.theta = prm.get_double ("theta");
   prm.leave_subsection();

   prm. enter_subsection("Finite Element");
      param_input.set_nq_points = prm.get_integer("n_quadrature");
      param_input.fe_degree = prm.get_integer("fe_degree");
   prm.leave_subsection();

}

code_parameters ParameterReader::get_parameters()
{
  declare_parameters();
  read_parameters();
  pass_parameters();
  return param_input;
}
