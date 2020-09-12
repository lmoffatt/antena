#ifndef ANTENA_DATA_H
#define ANTENA_DATA_H

#include "qm_vector_space.h"
#include "qm_vector_basis.h"
#include "qm_vector_field.h"
#include "qm_vector_index.h"
#include "qm_units_standard.h"
#include "antena_variables.h"
#include "qm_vector_index_variables.h"
#include "qm_mapu.h"
#include "qm_file_IO.h"


template<std::size_t J>
auto get_data_index(){

  using T=vs<
      x_i<i_i<delay,J>,v<double,ps_u>>,
      x_i<i_i<signal,J>,v<double,V_u>>
      >;

  using F=v_f<index_prod<i_i<delay,J>>,T>;

  using I=index_table<index_vector<ind<i_i<delay,J>,ind_size>>>;



  return mapu<I,F>();
}


template<std::size_t J>
auto get_file_data_index(const std::string& fname){

using T=vs<
      x_i<i_i<delay,J>,v<double,ps_u>>,
      x_i<i_i<signal,J>,v<double,V_u>>
      >;

using F=v_f<index_prod<Index<i_i<delay,J>>>,T>;

using I=index_table<index_vector<ind<Index<i_i<delay,J>>,ind_size>>>;



return file_IO<mapu<I,F>>(fname);
}




#endif // ANTENA_DATA_H
