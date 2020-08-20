/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2006 - 2015 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Wolfgang Bangerth, Texas A&M University, 2006
 */

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/point.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

#include <fstream>
#include <iostream>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include "Functions_SEM.h"
#include "Mesh_model.h"
#include "Param_handler.h"

  using namespace dealii;

  template <int dim>
  class SEM_topography
  {
  public:
    SEM_topography (const struct code_parameters & parameters);
    void run ();

  private:

    void setup_system ();
    void solve_u ();
    void solve_v ();
    void output_results () const;

    const int fe_degree;

    Triangulation<dim>   triangulation;
    FE_Q<dim>            fe;
    DoFHandler<dim>      dof_handler;

   // AffineConstraints<double> constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> laplace_matrix;
    SparseMatrix<double> matrix_old;

    Vector<double>       mass_vector;
    Vector<double>       boundary_vector;
    Vector<double>       vector_u1,vector_un, vector_old,vector_old2;
    Vector<double>       solution_u;
    Vector<double>       old_solution_u;
    Vector<double>       old_solution2_u;

    // Medium
   std::map<unsigned int, double> material_vel;

    //Source
    double fo;
    double source_rad;
    const Point<dim> source_pos;

    //receivers
    const int n_receivers;
    const double  y_pos;
    const double  x_init;
    const double  x_end;
    std::vector<Point<dim> > receiver_locations;

    //grid & time
    const std::string input_mesh;
    const int refinement;
    const unsigned int add_refinement; // refinement around source
    const double total_time;
    const double theta;
    double time, time_step;
    unsigned int timestep_number;

    const int set_nq_points;

    //timer file
    std::ofstream time_dof;
  };

// Constructor
  template <int dim>
  SEM_topography<dim>::SEM_topography (const struct code_parameters & parameters) :
    fe_degree(parameters.fe_degree),
    fe (fe_degree),
    dof_handler (triangulation),

    material_vel(parameters.material_vel),

    fo (parameters.fo),
    source_rad (parameters.source_rad),
    source_pos (parameters.source_pos_x, parameters.source_pos_y),

    n_receivers (parameters.n_receivers),
    y_pos (parameters.y_pos),
    x_init (parameters.x_init),
    x_end (parameters.x_end),

    input_mesh(parameters.input_mesh),
    refinement (parameters.refinement),
    add_refinement (parameters.add_refinement),
    total_time (parameters.total_time),
    theta (parameters.theta),
    time_step (parameters.time_step),

    set_nq_points (parameters.set_nq_points),

    time_dof ("time_dof.dat")
  {
    //**** filling vector receiver_locations****
      double receiver_interval = (x_end-x_init)/(n_receivers-1);

      for (int i=0; i<n_receivers; ++i)
         {
           Point<dim> p  ( x_init + i*receiver_interval, y_pos);
           receiver_locations.push_back(p);
         }

  }


  template <int dim>
  void SEM_topography<dim>::setup_system ()
  {
      Timer   reinit_timer, assembly_timer, factor_timer;
      reinit_timer.start();

    // importing created triangulation
    Mesh_utilities<dim> mesh;
    mesh.Mesh_import(input_mesh, triangulation);

   // const unsigned int cycle=0;
    //mesh.Mesh_info (triangulation, cycle);

    std::cout << "Number of active cells: "
              << triangulation.n_active_cells()
              << std::endl;
    dof_handler.distribute_dofs (fe);
      
    std::cout << "Number of degrees of freedom: "
              << dof_handler.n_dofs()
              << std::endl
              << std::endl;

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, dsp);
    sparsity_pattern.copy_from (dsp);

    // reinit:
    //*** Matrices***
   laplace_matrix.reinit (sparsity_pattern);
   matrix_old.reinit(sparsity_pattern);

    //Vectors:
    mass_vector.reinit(dof_handler.n_dofs());
    boundary_vector.reinit(dof_handler.n_dofs());

    vector_u1.reinit(dof_handler.n_dofs());
    vector_un.reinit(dof_handler.n_dofs());
    vector_old.reinit(dof_handler.n_dofs());
    vector_old2.reinit(dof_handler.n_dofs());

    solution_u.reinit (dof_handler.n_dofs());
    old_solution_u.reinit (dof_handler.n_dofs());
    old_solution2_u.reinit (dof_handler.n_dofs());

    reinit_timer.stop();
    time_dof << reinit_timer.cpu_time()<< " ";
    time_dof <<reinit_timer.wall_time()<< " 0.0" <<std::endl;

    std::cout<<"reinit cpu-time: "<< reinit_timer.cpu_time()<<" seconds"<< std::endl;
    std::cout<<"reinit wall-time: "<< reinit_timer.wall_time()<<" seconds"<< std::endl<< std::endl;

    //***Assembly***
    //**************
    assembly_timer.start();
       {
           QGaussLobatto<dim>  quadrature_formula (fe_degree+1);
           FEValues<dim> fe_values (fe, quadrature_formula,
                                    update_values    |  update_gradients |
                                    update_quadrature_points  |  update_JxW_values);
           const unsigned int   dofs_per_cell = fe.dofs_per_cell;
           const unsigned int   n_q_points    = quadrature_formula.size();

           Vector<double>        cell_vector_mass(dofs_per_cell);
           FullMatrix<double>   cell_matrix_laplace (dofs_per_cell, dofs_per_cell);

           std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
           double vp;

           //boundary
           const  QGaussLobatto<dim-1>  face_quadrature_formula(fe_degree+1);
           FEFaceValues<dim> face_fe_values (fe, face_quadrature_formula,
                                        update_values |  update_JxW_values);

           const unsigned int   face_q_points    = face_quadrature_formula.size();
           Vector<double> cell_boundary_vector(dofs_per_cell);

            //looping through cells
           typename DoFHandler<dim>::active_cell_iterator
           cell = dof_handler.begin_active(),
           endc = dof_handler.end();
           for (; cell!=endc; ++cell)
           {
               cell_vector_mass=0;
               cell_matrix_laplace= 0;

               fe_values.reinit (cell);

               vp = material_vel[cell->material_id()]; //obtaining p velocity for the cell

               for (unsigned int q_index=0; q_index<n_q_points; ++q_index)
               {
                   //std:: cout<<"laplace coef: "<<current_laplaceCoef<<std::endl;
                   for (unsigned int i=0; i<dofs_per_cell; ++i)
                   {
                    cell_vector_mass(i)+= (
                                        fe_values.shape_value(i,q_index) *
                                        fe_values.shape_value(i,q_index) *
                                        fe_values.JxW(q_index));

                       for (unsigned int j=0; j<dofs_per_cell; ++j)
                         {
                           cell_matrix_laplace(i,j) += (vp*vp*
                                                        fe_values.shape_grad(i,q_index) *
                                                        fe_values.shape_grad(j,q_index) *
                                                        fe_values.JxW(q_index));
                        //print local mass matrix
                          //std:: cout<<i <<"; "<<j <<"->"
                            //       <<cell_matrix_mass(i,j)<<std::endl;
                         }
                   }
               }

               cell->get_dof_indices (local_dof_indices);

               mass_vector.add (local_dof_indices,
                                   cell_vector_mass);
               for (unsigned int i=0; i<dofs_per_cell; ++i)
                   for (unsigned int j=0; j<dofs_per_cell; ++j)
                     {
                        laplace_matrix.add(local_dof_indices[i],
                                            local_dof_indices[j],
                                            cell_matrix_laplace(i,j));
                       }

               // query for boundary cell and boundary vector assembly
               for (unsigned int f=0; f< GeometryInfo<dim>::faces_per_cell; ++f)
                 if ((cell->at_boundary(f)) && (cell->face(f)->boundary_id()==10))
                   {
                     cell_boundary_vector=0;
                     face_fe_values.reinit (cell, f);

                     for (unsigned int q_point=0; q_point<face_q_points; ++q_point)
                       for (unsigned int i=0; i<dofs_per_cell; ++i)
                          {
                            cell_boundary_vector(i) += vp*(face_fe_values.shape_value(i,q_point) *
                                                   face_fe_values.shape_value(i,q_point) *
                                                   face_fe_values.JxW(q_point));
                          }

                     //updating global boundary vector
                     boundary_vector.add (local_dof_indices,
                                           cell_boundary_vector);
                   }

               }//end of cell loop
       }//end assembly
    assembly_timer.stop();
    time_dof << assembly_timer.cpu_time()<< " ";
    time_dof <<assembly_timer.wall_time()<< " 0.0" <<std::endl;

    std::cout<<"assembly cpu-time: "<< assembly_timer.cpu_time()<<" seconds"<< std::endl;
    std::cout<<"assembly wall-time: "<< assembly_timer.wall_time()<<" seconds"<< std::endl<< std::endl;

//*** non time -dependent matrices and vectors
 factor_timer.start();
    vector_u1=mass_vector;
    vector_u1 *=4;

    vector_un=mass_vector;
    vector_un *=2;
    vector_un.add(time_step,boundary_vector);

    //inverse of matrix_u1 and matrix_un
    for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
       {
        vector_u1(i) = 1/vector_u1(i);
        vector_un(i) = 1/vector_un(i);
         //std::cout<<i<<"; "<<i<<"->"<<old_solution_u(i)<<std::endl;
        }
    matrix_old.copy_from(laplace_matrix);
    matrix_old *=-2*time_step*time_step;

    vector_old=mass_vector;
    vector_old *=4;

    vector_old2 = mass_vector;
    vector_old2 *=-2;
    vector_old2.add(time_step, boundary_vector);

    factor_timer.stop();
    time_dof << factor_timer.cpu_time()<< " ";
    time_dof << factor_timer.wall_time()<< " 0.0" <<std::endl;

    std::cout<<"factorize cpu-time: "<< factor_timer.cpu_time()<<" seconds"<< std::endl;
    std::cout<<"factorize wall-time: "<< factor_timer.wall_time()<<" seconds"<< std::endl<< std::endl;
  }


  template <int dim>
  void SEM_topography<dim>::output_results () const
  {
    DataOut<dim> data_out;

    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (solution_u, "U");
    //data_out.add_data_vector (solution_v, "V");

    data_out.build_patches ();

    const std::string filename = "solution-" +
                                 Utilities::int_to_string (timestep_number, 3) +
                                 ".vtk";
    std::ofstream output (filename.c_str());
    data_out.write_vtk (output);
  }

  template <int dim>
  void SEM_topography<dim>::run ()
  {
    Timer timer, tot_timer;
    setup_system();

    Vector<double> tmp (solution_u.size());
    Vector<double> forcing_terms (solution_u.size());

    // Creating F vector
    RightHandSide<dim> rhs_function (fo, source_rad, source_pos);

    // open file to store receivers data
    std::ofstream receiver_data("receivers.dat");
    receiver_data.precision(15);

   std::cout<< "total time = " << total_time<<std::endl;

    timer.start();
    tot_timer.start();
    //***time stepping loop***
    for (timestep_number=1, time=time_step;
         time<=total_time;
         time+=time_step, ++timestep_number)
      {
        std::cout << "Time step " << timestep_number
                  << " at t= " << time
                  << std::endl;

   // creating F vector
        rhs_function.set_time (time-time_step);
        VectorTools::create_right_hand_side (dof_handler, QGaussLobatto<dim>(2*fe_degree),
                                             rhs_function, forcing_terms);
        forcing_terms *= 2*time_step*time_step;

        if (timestep_number==1)
          {
            solution_u=forcing_terms;
            solution_u.scale(vector_u1);
           }

         else
          {
            matrix_old.vmult(solution_u,old_solution_u);

            tmp=vector_old;
            tmp.scale(old_solution_u);
            solution_u.add(1,tmp);

            tmp=vector_old2;
            tmp.scale(old_solution2_u);
            solution_u.add(1,tmp);
            solution_u.add(1,forcing_terms);

            solution_u.scale(vector_un);

        }
//*********************************

        //cpu time
        timer.stop();
        time_dof << timer.cpu_time() << " " ;
        time_dof << timer.wall_time() << " " ;
        time_dof  << dof_handler.n_dofs()<< std::endl;

        std::cout<< "cpu time: " <<timer.cpu_time() << " seconds"<<std::endl;
        std::cout<< "wall time: " <<timer.wall_time() << " seconds"<<std::endl;
        std::cout << "dofs: "<<dof_handler.n_dofs()<<std::endl;


        //output_results ();

        // *** solution at receiver_locations **
        receiver_data << time;
                for (unsigned int i=0 ; i<receiver_locations.size(); ++i)
                  receiver_data << " "
                                << VectorTools::point_value (dof_handler,
                                                             solution_u,
                                                             receiver_locations[i])
                                << " ";
               receiver_data << std::endl;

        std::cout << "   Total energy: "
                  << (
                      laplace_matrix.matrix_norm_square (solution_u))
                  << std::endl;

        timer.restart();

        old_solution2_u=old_solution_u;
        old_solution_u = solution_u;
      }

    tot_timer.stop();
    time_dof <<tot_timer.cpu_time() << " ";
    time_dof<< tot_timer.wall_time() <<" 0.0";
    std::cout<<" solution total cpu-time: "<< tot_timer.cpu_time()<<" seconds"<<std::endl;
    std::cout<<" solution total wall-time: "<< tot_timer.wall_time()<<" seconds"<< std::endl<<std::endl;
  }

int main (int argc, char *argv[])
{
  try
    {
      using namespace dealii;
      //argv[1]
      //const char file[]="./../input_sem_topography.prm";
      ParameterReader read_param (argv[1]);
      code_parameters parameters =read_param.get_parameters();

      SEM_topography<2> wave_equation_solver(parameters);
      wave_equation_solver.run ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
