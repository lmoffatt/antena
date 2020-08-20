#ifndef ANTENA_PRIOR_VALUES_H
#define ANTENA_PRIOR_VALUES_H

#include "qm_vector_space.h"
#include "qm_vector_basis.h"
#include "antena_variables.h"
#include "qm_distribution_variables.h"
#include "qm_units_standard.h"
#include "qm_function_standard.h"

#include <cmath>
auto my_common_prior_values(){
  return vs{
      x_i(mean<baseline>{},v(0.0,V_u{})),
      x_i(stddev<baseline>{},v(1e-6,V_u{})),
      x_i(mean<drift>{},v(0.0,V_u{}/ps_u{})),
      x_i(stddev<drift>{},v(1e-7,V_u{}/ps_u{})),
      x_i(mean<Log10_t<stddev<signal>>>{},logv(-6.0,l_u<V_u>{})),
      x_i(stddev<Log10_t<stddev<signal>>>{},v(1.0)),
      x_i(mean<Log10_t<Amplitude>>{},logv(-5.0,l_u<V_u>{})),
      x_i(stddev<Log10_t<Amplitude>>{},v(1.0)),
      x_i(mean<Log10_t<i_i<Frecuency,0>>>{},logv(std::log10(0.0001),l_u<GHz_u>{})),
      x_i(stddev<Log10_t<i_i<Frecuency,0>>>{},v(0.2)),
      x_i(mean<Log10_t<Frecuency>>{},logv(std::log10(9),l_u<GHz_u>{})),
      x_i(stddev<Log10_t<Frecuency>>{},v(0.2)),
      x_i(mean<phase>{},v(PI/4)),
      x_i(stddev<phase>{},v(PI/4)),
      x_i(mean<Log10_t<tau>>{},logv(std::log10(1000),l_u<ps_u>{})),
      x_i(stddev<Log10_t<tau>>{},v(1.0))};
}

#endif // ANTENA_PRIOR_VALUES_H
