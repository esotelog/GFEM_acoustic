#include "Mesh_model.h"

using namespace dealii;


/*constructor
 **************************
 */
//template<int dim>
//Mesh_utilities<dim>::Mesh_utilities():
//{}

/* Import mesh
 **************************
 */
template<int dim>
void Mesh_utilities<dim>::Mesh_import (const std::string & input_mesh,
                                       Triangulation<dim> & triangulation)
{
  GridIn<dim> gridin;
  gridin.attach_triangulation(triangulation);
  std::ifstream f("./../" + input_mesh);

  gridin.read_msh(f);
}

/*refine_mesh
 **************
 Explicit single  inclusion refinement. Original mesh needs to grid the inclusion explicitly as well.
 */

template<int dim>
void Mesh_utilities<dim>::Mesh_refine (Triangulation<dim> & triangulation,
                                        const unsigned int & m_id)
{
  for (auto cell : triangulation.active_cell_iterators())
        {
              if ( (cell->material_id()==m_id)  && (cell->center()[1] <-1e-8 ))
                    cell->set_refine_flag();
          }

      triangulation.execute_coarsening_and_refinement();
}

/*Source_refine
 **************
Refinement around source
 */
template<int dim>
void Mesh_utilities<dim>::Source_refine (Triangulation<dim> & triangulation,
                  const Point<dim> & source_pos,
                  const double & source_r,
                  const unsigned int & refine)
{

double current_size;
  for (unsigned int it_refine =0; it_refine < refine; ++it_refine)
    {

        double max_size =GridTools::maximal_cell_diameter(triangulation);

         if (source_r < max_size )
           current_size=max_size;
         else
           current_size=source_r;

       //using C++ range based loops
      for (auto cell : triangulation.active_cell_iterators())
        {
              if ( source_pos.distance(cell->center()) <current_size)
                    cell->set_refine_flag();

          }

      triangulation.execute_coarsening_and_refinement();
    }
}




/*mesh_info
 **************************
 */
// this is a function to provide some basic information about the 
// triangulations. The code inside can be ignored. It produces an 
// eps image of the mesh so we can see the nodes are at the right 
// place.
//
template<int dim>
void Mesh_utilities<dim>::Mesh_info (const Triangulation<dim> & triangulation,
                                const unsigned int & cycle)
{
  {
    std::cout << "Mesh info:" << std::endl
              << " dimension: " << dim << std::endl
              << " no. of cells: " << triangulation.n_active_cells() << std::endl;
    {
      std::map<types::boundary_id, unsigned int> boundary_count;
      for (auto cell : triangulation.active_cell_iterators())
        {
          for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
            {
              if (cell->face(face)->at_boundary())
                boundary_count[cell->face(face)->boundary_id()]++;
            }
        }
      std::cout << " boundary indicators: ";
      for (const std::pair<types::boundary_id, unsigned int> &pair : boundary_count)
        {
          std::cout << pair.first << "(" << pair.second << " times) ";
        }
      std::cout << std::endl;
    }

    std::string filename = "mesh_model_";
    filename += ('0' + cycle);
    filename += ".msh";
    std::ofstream out (filename.c_str());
    GridOut grid_out;
    grid_out.write_msh (triangulation, out);
    std::cout << " written to " << filename
              << std::endl
              << std::endl;
    }

}

// instantiation

template class Mesh_utilities<2>;
