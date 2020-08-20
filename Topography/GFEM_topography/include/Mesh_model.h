#ifndef MESH
#define MESH

#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
//#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/base/logstream.h>

#include <fstream>
#include <iostream>
#include <cmath>



using namespace dealii;


//Mesh utilities
//*****************
template<int dim>
class Mesh_utilities
{
  public:
   // Mesh_utilities();

    void Mesh_import (const std::string & input_mesh,
                      Triangulation<dim> & triangulation );

    void Mesh_refine (Triangulation<dim> & triangulation,
                      const unsigned int & m_id);

    void Source_refine (Triangulation<dim> & triangulation,
                      const Point<dim> & source_pos,
                      const double & source_r,
                      const unsigned int & refine);

    void Mesh_info (const Triangulation<dim> & triangulation,
                   const unsigned int &cycle);
private:

};

#endif
