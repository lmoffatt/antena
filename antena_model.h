#ifndef ANTENA_MODEL_H
#define ANTENA_MODEL_H


#include "static_string.h"
#include "qm_distribution.h"
#include "qm_distribution_functions.h"
#include "qm_distribution_variables.h"
#include "antena_variables.h"
#include "qm_quimulun.h"
#include "qm_function_standard.h"
#include "qm_function.h"
#include "antena_prior_distributions.h"
#include "qm_units_standard.h"

template <std::size_t J, std::size_t...I0>
auto mymodel_index(std::index_sequence<I0...>,std::index_sequence<>)
{
  return
      quimulun{
          D(i_i<signal,J>{},Normal_Distribution(mean<i_i<signal,J>>{},stddev<i_i<signal,J>>{})),
          F(mean<i_i<signal,J>>{},
              ((i_i<baseline,J>{}
             +(i_i<drift,J>{}*i_i<delay,J>{}))+...+
                (i_i<Amplitude,I0,J>{}*Exp(-i_i<delay,J>{}/i_i<tau,I0,J>{})*
              Cos(N<2>{}*Pi{}*i_i<Frecuency,I0>{}*GHz_F{}*i_i<delay,J>{}*ps_F{}+i_i<phase,I0,J>{}))))
      };


}


template <std::size_t J, std::size_t...I0, std::size_t I00, std::size_t...I1>
auto mymodel_index(std::index_sequence<I0...>,std::index_sequence<I00,I1...>)
{
return quimulun{
          D(i_i<signal,J>{},Normal_Distribution(mean<i_i<signal,J>>{},stddev<i_i<signal,J>>{})),
          F(mean<i_i<signal,J>>{},
            ((((i_i<baseline,J>{}
             +(i_i<drift,J>{}*i_i<delay,J>{}))+...
                +(i_i<Amplitude,I0,J>{}*Exp(-i_i<delay,J>{}/i_i<tau,I0,J>{})*
                Cos(N<2>{}*Pi{}*i_i<Frecuency,I0>{}*GHz_F{}*i_i<delay,J>{}*ps_F{}+i_i<phase,I0,J>{})))
                +(i_i<Amplitude,I00,J>{}*Exp(-i_i<delay,J>{}/i_i<tau,I00,J>{})*
                       Cos(N<2>{}*Pi{}*i_i<Frecuency,I00>{}*GHz_F{}*i_i<delay,J>{}*ps_F{}+i_i<phase,I00,J>{})))+...
                +(i_i<Amplitude,I1,J>{}*Exp(-i_i<delay,J>{}/i_i<tau,I1,J>{})*
            Cos(N<2>{}*Pi{}*i_i<Frecuency,I1>{}*GHz_F{}*i_i<delay,J>{}*ps_F{}+i_i<phase,I1,J>{}))))
};


}

template <std::size_t... J, std::size_t...I0, std::size_t...I1>
auto mymodel_index(std::index_sequence<J...>,std::index_sequence<I0...> i0,std::index_sequence<I1...> i1)
{
  return (mymodel_index<J>(i0,i1)&&...);
}


template <std::size_t... J, std::size_t...I0, std::size_t...I1>
auto totalmodel_index(std::index_sequence<J...> j,std::index_sequence<I0...> i0,std::index_sequence<I1...> i1)
{
  return mymodel_index(j,i0,i1)&&myprior_dist_index(j,i0,i1)
         &&myprior_transf_index(j,i0,i1);
  //+my_common_prior_values();
}

template <std::size_t... J, std::size_t...I0>
auto totalmodel_index(std::index_sequence<J...> j,std::index_sequence<I0...> i0,std::index_sequence<> i1)
{
  return mymodel_index(j,i0,i1)&&
         myprior_dist_index(j,i0,i1)&&
         myprior_transf_index(j,i0,i1);
  //+my_common_prior_values();
}

#endif // ANTENA_MODEL_H
