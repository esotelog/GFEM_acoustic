#ifndef FUNCTIONS
#define FUNCTIONS

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/point.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>

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
#include <math.h>

using namespace dealii;



/* class Enrichment()
 * ************************
 */

template <int dim>
class Enrichment : public Function<dim>
{
public:

 Enrichment (const double & _k,
             const unsigned int & _n,
             const unsigned int & _q,
             const double & _phase = 0); // constructor


 virtual double value (const Point<dim>   &p,
                        const unsigned int  component = 0) const; //value

 virtual Tensor<1,dim> gradient (const Point<dim>   &p,
                                    const unsigned int  component = 0) const; //gradient


private:

    double k;//coefficient variable
    unsigned int n,q;
    double phase;
    const double PI;
};


// variable Coefficient (P velocity) for the Laplace Matrix
//**********************************************************
template <int dim>
class LaplaceCoefficient : public Function<dim>
{
public:

    LaplaceCoefficient(const std::vector<double> & Vp_,
                        const std::vector<double> & h_
                        );

    virtual void value_list (const std::vector<Point<dim> > &points,
                             std::vector<double>           &values,
                             const unsigned int            component = 0) const;
    virtual double value (const Point<dim>   &p,
                          const unsigned int  component=0 ) const;
private:
    const std::vector< double> Vp;
    const std::vector< double> h;
};

// source function (Right Hand side)
//****************************************
template <int dim>
class RightHandSide : public Function<dim>
{
public:
    RightHandSide ( const double & fo_,
                   const double & source_rad_,
                   const Point<dim> & source_pos_
                   );
    
    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
private:
    double fo;
    double source_rad;
    const Point<dim> source_pos;
    
};

// Initial value for u and v at time=0
//****************************************
template <int dim>
class InitialValuesU : public Function<dim>
{
public:
    InitialValuesU ();//constructor
    
    virtual double value (const Point<dim>   &p,
                          const unsigned int  component=0 ) const;
};

template <int dim>
class InitialValuesV : public Function<dim>
{
public:
    InitialValuesV ();
    
    virtual double value (const Point<dim>   &p,
                          const unsigned int  component =0) const;
};


#endif
