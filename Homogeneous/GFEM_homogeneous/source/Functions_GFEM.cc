#include "Functions_GFEM.h"

using namespace dealii;


/* Class Enrichment
 ******************
 * Plane wave enrichment
 */

template <int dim>
Enrichment<dim>::Enrichment (const double & _k,
                             const unsigned int & _n,
                             const unsigned int & _q,
                             const double & _phase) :

    Function<dim>(),
    k(_k), n(_n), q(_q),phase(_phase),
    PI(numbers::PI)
{}


template <int dim>
double Enrichment<dim>::value (const Point<dim> &p,
                                   const unsigned int ) const
{

  return cos(k* (p[0] * cos (2*PI*n/q) + p[1] * sin (2*PI*n/q) ) + phase );

}

template <int dim>
Tensor<1,dim> Enrichment<dim>::gradient (const Point<dim>   &p,
                                       const unsigned int) const
{
    double c = -k*sin(k* (p[0] * cos (2*PI*n/q) + p[1] * sin (2*PI*n/q) ) + phase);
    return Tensor<1,dim> ({ c*cos(2*PI*n/q) , c*sin(2*PI*n/q)});
}


// source function (Right Hand side)
//****************************************

// Constructor
template <int dim>
RightHandSide<dim>::RightHandSide ( const double & fo_,
                                    const double & source_rad_,
                                    const Point<dim> & source_pos_
                                    ):
  Function<dim>(),
  fo (fo_),
  source_rad (source_rad_),
  source_pos (source_pos_)
  {}


template <int dim>
double RightHandSide<dim>::value (const Point<dim>  &p,
                                  const unsigned int component) const
{
    Assert (component == 0, ExcInternalError());
    //source period
    const double to=1.0/fo;
    
    
    // source definition:

    double  spatial_term=0.0, temporal_term =0.0;
    double t=this->get_time();
    
    double source_rad2=3.125;

    if ((t<= 2*to) && // source temporal extent
        // source spatial extent (x-xo)^2 +(y-yo)^2 <=source_rad^2
        (source_pos.distance(p) <=source_rad) )
    {
        double c,d,fc,vol;

        c=(numbers::PI*numbers::PI) * (fo*fo) * ((t-to)*(t-to));
        fc=exp(-0.5)/(sqrt(2.0)*numbers::PI*fo);
        temporal_term=(t-to)* exp (-c)/fc;

        d=source_pos.distance(p);
        spatial_term =pow(1.0-1.0/(source_rad*source_rad) * (d*d ), 3);
        vol=numbers::PI/4.0*source_rad2*source_rad2;
        
         return -500000.0*10.0*spatial_term*temporal_term/vol;
    }
    
    //
    else
        return 0.0;
}

// Initial value for u and v at time=0
//****************************************

//constructor for U
template <int dim>
InitialValuesU<dim>::InitialValuesU (): Function<dim>() {}

template <int dim>
double InitialValuesU<dim>::value (const Point<dim>  &/*p*/,
                                   const unsigned int component) const
{
    Assert (component == 0, ExcInternalError());
    return 0;
}

//constructor for V
template <int dim>
InitialValuesV<dim>::InitialValuesV (): Function<dim>() {}

template <int dim>
double InitialValuesV<dim>::value (const Point<dim>  &/*p*/,
                                   const unsigned int component) const
{
    Assert (component == 0, ExcInternalError());
    return 0;
}

// instantiation
template class Enrichment<2>;
template class RightHandSide<2>;
template class InitialValuesU<2>;
template class InitialValuesV <2>;


