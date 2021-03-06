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
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/sparse_direct.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_enriched.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

#include <math.h>
#include <fstream>
#include <iostream>
# include <utility>
#include <algorithm>

#include "Functions_GFEM.h"
#include "Mesh_model.h"
#include "Param_handler.h"

  using namespace dealii;

  template <int dim>
  class GFEM_LowVelocity
  {
  public:
    GFEM_LowVelocity (const struct code_parameters & parameters);
    ~GFEM_LowVelocity ();
    void run ();

  private:

    void setup_system ();
    void solve_u ();
    void solve_v ();
    void output_results () const;

    Triangulation<dim>   triangulation;
    DoFHandler<dim>      dof_handler;
    std::shared_ptr<FE_Enriched<dim>> fe_ptr; // smart pointer to FE_Enriched

    AffineConstraints<double> constraints;
    SparseDirectUMFPACK  U_direct;
    SparseDirectUMFPACK  V_direct;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> mass_matrix;
    SparseMatrix<double> laplace_matrix;
    SparseMatrix<double> matrix_u;
    SparseMatrix<double> matrix_v;
    SparseMatrix<double> matrix1_old, matrix2_old, matrix3;
    SparseMatrix<double> boundary_matrix;// boundary matrix

    Vector<double>       solution_u, solution_v;
    Vector<double>       old_solution_u, old_solution_v;
    Vector<double>       system_rhs;

    //Medium
    std::map<unsigned int, double> material_vel;

    // Source
    double fo;
    double source_rad;
    const Point<dim> source_pos;

    //receivers
    const int n_receivers;
    const double  y_pos;
    const double  x_init;
    const double  x_end;

    //grid & time
    const std::string input_mesh;
    const int refinement;
    const unsigned int source_refinement;
    const double total_time;
    const double theta;
    double time, time_step;
    unsigned int timestep_number;

    const unsigned int fe_degree;
    const int set_nq_points;
    const unsigned int q_directions;

    std::vector<Point<dim> > receiver_locations;

    //timer-dof file
    std::ofstream time_dof;


  };

// Constructor

  template <int dim>
  GFEM_LowVelocity<dim>::GFEM_LowVelocity (const struct code_parameters & parameters) :

    triangulation(Triangulation<dim>::maximum_smoothing),
    dof_handler (triangulation),

    material_vel(parameters.material_vel),
    fo (parameters.fo),
    source_rad (parameters.source_rad),
    source_pos (parameters.source_pos_x, parameters.source_pos_y),

    n_receivers (parameters.n_receivers),
    y_pos (parameters.y_pos),
    x_init (parameters.x_init),
    x_end (parameters.x_end),

  //Discretization: grid & time
    input_mesh(parameters.input_mesh),
    refinement (parameters.refinement),
    source_refinement(parameters.source_refinement),
    total_time (parameters.total_time),
    theta (parameters.theta),
    time_step (parameters.time_step),

    fe_degree(parameters.fe_degree),
    set_nq_points (parameters.set_nq_points),
    q_directions(parameters.q_directions),

  // timer-dof file
   time_dof ("time_dof.dat")


  {
    //  Enriched Finete element smart pointer (fe_ptr) definition
    //**********************************************************

    // pi is a constant defined in dealII.
    const double pi =numbers::PI;

    //  deference of iterator  (key,value)  corresponding to the min value in the map material_vel
    auto min_vp=  *std::min_element(material_vel.begin(),material_vel.end(),
                                   [](auto i, auto j)->bool {return i.second < j.second;});

    std::vector<double> k_vector({2*pi*fo/(min_vp.second)});

    // Vector of enrichment functions
    std::vector< std::function< const Function<dim> * (const typename Triangulation <dim>::cell_iterator &) >  > functions;

        for (unsigned int n = 0; n < q_directions; ++n)
          for (unsigned int nk = 0; nk < k_vector.size(); ++nk)

          {
            std::shared_ptr<Enrichment<dim>> enrich_ptr(new Enrichment<dim> (k_vector[nk],n,q_directions));

            functions.push_back([=] (const typename Triangulation<dim>::cell_iterator &) -> const Function<dim> * {return enrich_ptr.get();});

           }

        //finite element for standard and enriched part have the same polynomial degree.
       const std::shared_ptr<FE_Q<dim>> feq_ptr(new FE_Q<dim> (fe_degree)) ;

       //FE_Enriched<dim> fe_enriched( feq_ptr.get(), {feq_ptr.get()}, {functions});

       //creating smart pointer for enriched fe
       fe_ptr=std::make_shared<FE_Enriched<dim>,
                  const FiniteElement< dim> *,
                  const std::vector< const FiniteElement< dim > * >  ,
                  const std::vector< std::vector< std::function< const Function< dim > *(const typename Triangulation< dim>::cell_iterator &) > > >
                    >
           ( feq_ptr.get(), {feq_ptr.get()}, {functions});


       // filling vector receiver_locations
       //***********************************
       double receiver_interval = (x_end-x_init)/(n_receivers-1);

       for (int i=0; i<n_receivers; ++i)
          {
            Point<dim> p  ( x_init + i*receiver_interval, y_pos);
            //std::cout <<i<<"->"<<x_init + i*receiver_interval<<std::endl;
            receiver_locations.push_back(p);
          }

  }

  /* Destructor
   * ************************
   */
  template <int dim>
  GFEM_LowVelocity<dim>::~GFEM_LowVelocity()
  {
    dof_handler.clear ();
  }

  template <int dim>
  void GFEM_LowVelocity<dim>::setup_system ()
  {
    //timer start
    Timer   reinit_timer, assembly_timer, factor_timer;
    reinit_timer.start();

    Mesh_utilities<dim> mesh;
    mesh.Mesh_import(input_mesh, triangulation);

    //refine first layer
    const unsigned int  m_id=6;
    mesh.Mesh_refine(triangulation,m_id);

    // refinement around source
    if (source_refinement>0)
      mesh.Source_refine(triangulation, source_pos, source_rad,source_refinement);

     const unsigned int cycle=0;
     mesh.Mesh_info (triangulation, cycle);

//    time_step = GridTools::minimal_cell_diameter(triangulation) /
//               fast_Vp /
//                std::sqrt (1.*dim);

    std::cout << "Number of active cells: "
              << triangulation.n_active_cells()
              << std::endl;
      
    dof_handler.distribute_dofs (*fe_ptr);

    std::cout << "Number of degrees of freedom: "
              << dof_handler.n_dofs()
              << std::endl
              << std::endl;

    constraints.clear ();
    DoFTools::make_hanging_node_constraints (dof_handler,
                                                constraints);
    constraints.close();

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler,  dsp, constraints,true);
    sparsity_pattern.copy_from (dsp);

    //*** Matrices***
    mass_matrix.reinit (sparsity_pattern);
    laplace_matrix.reinit (sparsity_pattern);
    matrix_u.reinit (sparsity_pattern);
    matrix_v.reinit (sparsity_pattern);
    boundary_matrix.reinit (sparsity_pattern);

    matrix1_old.reinit(sparsity_pattern);
    matrix2_old.reinit(sparsity_pattern);
    matrix3.reinit(sparsity_pattern);

    //***Vectors***
    solution_u.reinit (dof_handler.n_dofs());
    solution_v.reinit (dof_handler.n_dofs());
    old_solution_u.reinit (dof_handler.n_dofs());
    old_solution_v.reinit (dof_handler.n_dofs());
    system_rhs.reinit (dof_handler.n_dofs());

    reinit_timer.stop();
    time_dof << reinit_timer.cpu_time()<< " ";
    time_dof <<reinit_timer.wall_time()<< " 0.0" <<std::endl;

    std::cout<<"reinit cpu-time: "<< reinit_timer.cpu_time()<<" seconds"<< std::endl;
    std::cout<<"reinit wall-time: "<< reinit_timer.wall_time()<<" seconds"<< std::endl<< std::endl;

    //***Assembly***
    //**************
    assembly_timer.start();
      {
          QGauss<dim>  quadrature_formula (set_nq_points);
          FEValues<dim> fe_values (*fe_ptr, quadrature_formula,
                                   update_values    |  update_gradients |
                                   update_quadrature_points  |  update_JxW_values);
          const unsigned int   dofs_per_cell = fe_ptr->dofs_per_cell;
          const unsigned int   n_q_points    = quadrature_formula.size();
          FullMatrix<double>   cell_matrix_mass (dofs_per_cell, dofs_per_cell);
          FullMatrix<double>   cell_matrix_laplace (dofs_per_cell, dofs_per_cell);
          
          std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
          double vp;

          //boundary
          const QGauss<dim-1>  face_quadrature_formula(set_nq_points);
          FEFaceValues<dim> face_fe_values (*fe_ptr, face_quadrature_formula,
                                       update_values |  update_JxW_values);
          const unsigned int   face_nq_points = face_quadrature_formula.size();
          FullMatrix<double>   cell_boundary_matrix (dofs_per_cell, dofs_per_cell);

          //looping through cells
          typename DoFHandler<dim>::active_cell_iterator
          cell = dof_handler.begin_active(),
          endc = dof_handler.end();
          for (; cell!=endc; ++cell)
          {
              cell_matrix_mass= 0;
              cell_matrix_laplace= 0;
              
              fe_values.reinit (cell);
              vp = material_vel[cell->material_id()]; //obtaining p velocity for the cell

                  for (unsigned int i=0; i<dofs_per_cell; ++i)               
                      for (unsigned int j=0; j<dofs_per_cell; ++j)                     
                         for (unsigned int q_index=0; q_index<n_q_points; ++q_index)
                           {
                             cell_matrix_mass(i,j) += (
                                               fe_values.shape_value(i,q_index) *
                                               fe_values.shape_value(j,q_index) *
                                               fe_values.JxW(q_index));

                             cell_matrix_laplace(i,j) += (vp*vp*
                                                       fe_values.shape_grad(i,q_index) *
                                                       fe_values.shape_grad(j,q_index) *
                                                       fe_values.JxW(q_index));
                            }


                  // obtaining global dofs
                  cell->get_dof_indices (local_dof_indices);

                  //filling global matrices
                  for (unsigned int i=0; i<dofs_per_cell; ++i)
                      for (unsigned int j=0; j<dofs_per_cell; ++j)
                        {
                           mass_matrix.add (local_dof_indices[i],
                                            local_dof_indices[j],
                                            cell_matrix_mass(i,j));

                           laplace_matrix.add(local_dof_indices[i],
                                             local_dof_indices[j],
                                             cell_matrix_laplace(i,j));
                         }


                  //query if cell at boundary, then compute the local boundary matrix
                   for (unsigned int f=0; f< GeometryInfo<dim>::faces_per_cell; ++f)
                      if (cell->at_boundary(f))
                         {
                          cell_boundary_matrix = 0;
                          face_fe_values.reinit (cell, f);

                          for (unsigned int i=0; i<dofs_per_cell; ++i)
                             for (unsigned int j=0; j<dofs_per_cell; ++j)
                                for (unsigned int q_point=0; q_point<face_nq_points; ++q_point)
                                    cell_boundary_matrix(i,j) += vp*(face_fe_values.shape_value(i,q_point) *
                                                                    face_fe_values.shape_value(j,q_point) *
                                                                    face_fe_values.JxW(q_point));
                          //filling global boundary matrix
                          for (unsigned int i=0; i<dofs_per_cell; ++i)
                             for (unsigned int j=0; j<dofs_per_cell; ++j)
                                 boundary_matrix.add (local_dof_indices[i],
                                                     local_dof_indices[j],
                                                     cell_boundary_matrix(i,j));
                         }

          }//end of cell loop
          
      } //end of assembly

    assembly_timer.stop();
    time_dof << assembly_timer.cpu_time()<< " ";
    time_dof <<assembly_timer.wall_time()<< " 0.0" <<std::endl;

    std::cout<<"assembly cpu-time: "<< assembly_timer.cpu_time()<<" seconds"<< std::endl;
    std::cout<<"assembly wall-time: "<< assembly_timer.wall_time()<<" seconds"<< std::endl<< std::endl;

    // non-time dependent matrices
  factor_timer.start();
        //For solving u:
        matrix1_old.copy_from (mass_matrix);
        matrix1_old.add(time_step*theta, boundary_matrix);
        matrix1_old.add(- time_step*time_step *theta*(1-theta),laplace_matrix);

        //matrix_u
        matrix_u.copy_from (mass_matrix);
        matrix_u.add (theta * theta * time_step * time_step, laplace_matrix);
        matrix_u.add(time_step*theta,boundary_matrix);
        constraints.condense(matrix_u); // adding constraints

        //For solving v:
        matrix2_old.copy_from (boundary_matrix);
        matrix2_old.add(-time_step*(1-theta),laplace_matrix);

        matrix3.copy_from (boundary_matrix);
        matrix3.add(time_step*theta, laplace_matrix);

        //matrix_v
        matrix_v.copy_from(mass_matrix);
        constraints.condense(matrix_v); // adding constraints

       //initializing and factorizing left handside matrices for the direct solver
         U_direct.initialize(matrix_u);
         V_direct.initialize(matrix_v);

       //stop timer and writing to file
          factor_timer.stop();
          time_dof << factor_timer.cpu_time()<< " ";
          time_dof << factor_timer.wall_time()<< " 0.0" <<std::endl;

          std::cout<<"factorize cpu-time: "<< factor_timer.cpu_time()<<" seconds"<< std::endl;
          std::cout<<"factorize wall-time: "<< factor_timer.wall_time()<<" seconds"<< std::endl<< std::endl;

  }

  template <int dim>
  void GFEM_LowVelocity<dim>::solve_u ()
  {
      U_direct.vmult (solution_u, system_rhs);
      constraints.distribute(solution_u);
  }

  template <int dim>
  void GFEM_LowVelocity<dim>::solve_v ()
  {
      V_direct.vmult (solution_v, system_rhs);
      constraints.distribute(solution_v);
  }

  template <int dim>
  void GFEM_LowVelocity<dim>::output_results () const
  {
    DataOut<dim> data_out;

    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (solution_u, "U");
    data_out.add_data_vector (solution_v, "V");

    data_out.build_patches ();

    const std::string filename = "solution-" +
                                 Utilities::int_to_string (timestep_number, 3) +
                                 ".vtk";
    std::ofstream output (filename.c_str());
    data_out.write_vtk (output);
  }


  template <int dim>
  void GFEM_LowVelocity<dim>::run ()
  {
    Timer timer, tot_timer;

    setup_system();

    RightHandSide<dim> rhs_function (fo, source_rad, source_pos);

    Vector<double> tmp (solution_u.size());
    Vector<double> forcing_terms (solution_u.size());

    std::ofstream receiver_data("receivers.dat"); // open file to store receivers data
    receiver_data.precision(15);

    std::cout << "total simulation time = " << total_time<<std::endl;
    timer.start();
    tot_timer.start();

    //***time stepping loop***
    for (timestep_number=1, time=time_step;
         time<=total_time;
         time+=time_step, ++timestep_number)
      {
        std::cout << "Time step " << timestep_number
                  << " at t=" << time
                  << std::endl;

        //***Solving for u**
        // Right hand side system (system_rhs)

        matrix1_old.vmult(system_rhs,old_solution_u);
        mass_matrix.vmult(tmp,old_solution_v);
        system_rhs.add(time_step,tmp);

        rhs_function.set_time (time);
        VectorTools::create_right_hand_side (dof_handler, QGauss<dim>(set_nq_points),
                                                     rhs_function, tmp);
        forcing_terms = tmp;
        forcing_terms *= theta * time_step;

        rhs_function.set_time (time-time_step);
        VectorTools::create_right_hand_side (dof_handler, QGauss<dim>(set_nq_points),
                                             rhs_function, tmp);

        forcing_terms.add ((1-theta) * time_step, tmp);

        system_rhs.add (theta * time_step, forcing_terms);
        constraints.condense(system_rhs);

        solve_u ();

        //**Solving for v**
         // Right hand side system (system_rhs)

        matrix2_old.vmult(system_rhs,old_solution_u);
        matrix3.vmult(tmp,solution_u);
        system_rhs.add(-1,tmp);
        mass_matrix.vmult (tmp, old_solution_v);
        system_rhs.add(1,tmp);

        system_rhs += forcing_terms;
        constraints.condense(system_rhs);

        solve_v ();

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
                  << (mass_matrix.matrix_norm_square (solution_v) +
                      laplace_matrix.matrix_norm_square (solution_u)) / 2
                  << std::endl<<std::endl;
      
        timer.restart();

        old_solution_u = solution_u;
        old_solution_v = solution_v;
      }

    tot_timer.stop();
    time_dof <<tot_timer.cpu_time() << " ";
    time_dof<< tot_timer.wall_time() <<" 0.0";
    std::cout<<"total cpu time: "<< tot_timer.cpu_time()<<" seconds"<< std::endl;
    std::cout<<" solution total wall-time: "<< tot_timer.wall_time()<<" seconds"<< std::endl<<std::endl;

  }

int main (int argc, char *argv[])
{
  try
    {
      using namespace dealii;
      //argv[1]
      //const char file[]="./../input_gfem_LowVelocity.prm";
      ParameterReader read_param (argv[1]);
      code_parameters parameters =read_param.get_parameters();

      GFEM_LowVelocity<2> wave_equation_solver(parameters);
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
