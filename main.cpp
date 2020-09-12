#include "qm_unit.h"
#include "qm_vector_index_basis.h"
#include "qm_vector_field.h"
#include <iostream>
#include "antena_model.h"
#include "antena_data.h"
#include "antena_prior_values.h"
#include "qm_distribution_engine.h"
#include "qm_function_engine.h"
#include "qm_quimulun_engine.h"

struct my_engine: public Op_engine, public distribution_engine
{
    using Op_engine::operator[];
    using distribution_engine::operator[];
};




int main(int argc, char **argv)
{
 auto model_8910_012_=totalmodel_index(std::index_sequence<8,9,10>{},std::index_sequence<0,1,2>{},std::index_sequence<>{});
 auto data_5=extract(get_file_data_index<5>("data"))+my_common_prior_values();
 auto data_8=extract(get_file_data_index<8>("data"));
 auto data_9=extract(get_file_data_index<9>("data"));
 auto data_10=extract(get_file_data_index<10>("data"));
 get_file_data_index<5>("data_copy_")<<data_5;
 get_file_data_index<8>("data_copy_")<<data_8;
 get_file_data_index<9>("data_copy_")<<data_9;
 get_file_data_index<10>("data_copy_")<<data_10;

 auto data_8910=join(std::move(data_8),std::move(data_9),std::move(data_10))+my_common_prior_values();

 //std::cerr<<my_common_prior_values();


 auto sample_8910=execute(my_engine{},Replicate{},model_8910_012_,data_8910);
 std::cerr<<sample_8910.name;

  std::cerr<<model_8910_012_.name.c_str();

  return 0;
}
