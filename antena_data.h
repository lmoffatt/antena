#ifndef ANTENA_DATA_H
#define ANTENA_DATA_H

#include "qm_vector_space.h"
#include "qm_vector_basis.h"
#include "antena_variables.h"

template<std::size_t J>
auto get_data_index(){return vs{x_i(i_i<delay,J>{}, vec<i_i<delay,J>>{}),
                                            x_i(i_i<signal,J>{},vec<i_i<delay,J>>{})};}

#endif // ANTENA_DATA_H
