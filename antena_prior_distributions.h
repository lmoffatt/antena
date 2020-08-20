#ifndef ANTENA_PRIOR_DISTRIBUTIONS_H
#define ANTENA_PRIOR_DISTRIBUTIONS_H

#include "qm_distribution.h"
#include "qm_distribution_variables.h"
#include "qm_quimulun.h"
#include "antena_variables.h"
#include "qm_function_standard.h"




template<std::size_t J>
constexpr auto myprior_dist_stddev()
{
  return quimulun{
      D(i_i<baseline,J>{},Normal_Distribution(mean<baseline>{},stddev<baseline>{})),
      D(i_i<drift,J>{},Normal_Distribution(mean<drift>{},stddev<drift>{})),
      D(Log10_t<stddev<i_i<signal,J>>>{},Normal_Distribution(mean<Log10_t<stddev<signal>>>{},stddev<Log10_t<stddev<signal>>>{}))
  };
}

template<std::size_t I>
constexpr auto myprior_dist_frequency()
{
  if constexpr (I>0)
    return quimulun{
        D(i_i<Log10_t<Frecuency>,I>{},Normal_Distribution(mean<Log10_t<Frecuency>>{},stddev<Log10_t<Frecuency>>{}))};
  else  return quimulun{
        D(i_i<Log10_t<Frecuency>,I>{},Normal_Distribution(mean<Log10_t<i_i<Frecuency,0>>>{},stddev<Log10_t<Frecuency>>{}))};

}

template<std::size_t I, std::size_t J>
constexpr auto myprior_dist_frequency()
{
  return quimulun{
      D(i_i<Log10_t<Frecuency>,I,J>{},Normal_Distribution(mean<Log10_t<Frecuency>>{},stddev<Log10_t<Frecuency>>{}))};
}


template<std::size_t I, std::size_t J>
auto myprior_dist_amplitude()
{
  return quimulun{
      D(i_i<Log10_t<Amplitude>,I,J>{},Normal_Distribution(mean<Log10_t<Amplitude>>{},stddev<Log10_t<Amplitude>>{})),
      D(i_i<phase,I,J>{},Normal_Distribution(mean<phase>{},stddev<phase>{})),
      D(i_i<Log10_t<tau>,I,J>{},Normal_Distribution(mean<Log10_t<tau>>{},stddev<Log10_t<tau>>{}))};
}

template <std::size_t J, std::size_t...I0>
auto myprior_dist_index(std::index_sequence<I0...> ,std::index_sequence<> )
{
  return myprior_dist_stddev<J>()
         &&(myprior_dist_amplitude<I0,J>()&&...);
}

template <std::size_t J, std::size_t...I0, std::size_t...I1>
auto myprior_dist_index(std::index_sequence<I0...> ,std::index_sequence<I1...> )
{
  return myprior_dist_stddev<J>()&&(myprior_dist_frequency<I1,J>()&&...)
         &&(myprior_dist_amplitude<I0,J>()&&...)&&(myprior_dist_amplitude<I1,J>()&&...);
}
template <std::size_t... J, std::size_t...I0, std::size_t...I1>
auto myprior_dist_index(std::index_sequence<J...>,std::index_sequence<I0...> i0,std::index_sequence<I1...> i1)
{
  return (myprior_dist_index<J>(i0,i1)&&...)&&(myprior_dist_frequency<I0>()&&...);
}

template<std::size_t J>
auto myprior_transf_stddev(){return quimulun{
      F(stddev<i_i<signal,J>>{},Exp10(Log10_t<stddev<i_i<signal,J>>>{}))
  };
}

    template<std::size_t J,std::size_t... I>
    auto myprior_transf_Amplitude(std::index_sequence<I...>)
{
  return quimulun{
          F(i_i<Amplitude,I,J>{},Exp10(i_i<Log10_t<Amplitude>,I,J>{}))...,
      F(i_i<tau,I,J>{},Exp10(i_i<Log10_t<tau>,I,J>{}))...,
      };
}

template<std::size_t I>
auto myprior_transf_frequency()
{
  return quimulun{
      F(i_i<Frecuency,I>{},Exp10(i_i<Log10_t<Frecuency>,I>{}))
  };
}
template<std::size_t I, std::size_t J>
auto myprior_transf_frequency()
{
  return quimulun{
      F(i_i<Frecuency,I,J>{},Exp10(i_i<Log10_t<Frecuency>,I,J>{}))
  };
}

template <std::size_t J, std::size_t...I0>
auto myprior_transf_index(std::index_sequence<I0...> i0,std::index_sequence<> )
{
  return myprior_transf_stddev<J>()&&
         myprior_transf_Amplitude<J>(i0);
}


template <std::size_t J, std::size_t...I0, std::size_t I00,std::size_t...I1>
auto myprior_transf_index(std::index_sequence<I0...> i0,std::index_sequence<I00,I1...> i1)
{
  return myprior_transf_stddev<J>()&&(myprior_transf_frequency<I00,J>()&&...&&myprior_transf_frequency<I1,J>())&&
         myprior_transf_Amplitude<J>(i0)&&myprior_transf_Amplitude<J>(i1);
}
template <std::size_t... J, std::size_t...I0, std::size_t...I1>
auto myprior_transf_index(std::index_sequence<J...>,std::index_sequence<I0...> i0,std::index_sequence<I1...> i1)
{
  return (myprior_transf_index<J>(i0,i1)&&...)&&(myprior_transf_frequency<I0>()&&...);
}


#endif // ANTENA_PRIOR_DISTRIBUTIONS_H
