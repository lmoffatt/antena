#include "qm_unit.h"
#include "qm_vector_index_basis.h"
#include "qm_vector_field.h"
#include <iostream>
#include "antena_model.h"
#include "antena_data.h"
#include "antena_prior_values.h"




int main(int argc, char **argv)
{
 auto model_8910_012_=totalmodel_index(std::index_sequence<8,9,10>{},std::index_sequence<0,1,2>{},std::index_sequence<>{});
 auto data=extract(get_file_data_index<5>("data"));
 get_file_data_index<5>("data_copy_")<<data;


  std::cerr<<model_8910_012_.name.c_str();
  auto p=my_common_prior_values();

  return 0;
}
