#include "qm_unit.h"
#include "qm_vector_index.h"
#include "qm_vector_field.h"
#include <iostream>
#include "antena_model.h"
#include "antena_prior_values.h"




int main(int argc, char **argv)
{
  auto r=p_t<u<m,1>,u<s,-1>,u<A,2>,u<kg,0>>{};
  std::cerr<<r.name.c_str()<<"\n";
  auto model_8910_012_=totalmodel_index(std::index_sequence<8,9,10>{},std::index_sequence<0,1,2>{},std::index_sequence<>{});
  std::cerr<<model_8910_012_.name.c_str();
  auto p=my_common_prior_values();

  return 0;
}
