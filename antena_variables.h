#ifndef ANTENA_VARIABLES_H
#define ANTENA_VARIABLES_H



#include "static_string.h"



struct delay {  constexpr static auto  name=my_static_string("delay"); };




struct signal { constexpr static auto  name=my_static_string("signal"); };

template <class ei, std::size_t... I> struct i_i {
  constexpr static auto  name=(ei::name+...+(my_static_string("_")+to_static_string<I>()));
};



struct baseline{ constexpr static auto name=my_static_string("baseline");};

struct drift{constexpr static auto name=my_static_string("drift");};

struct Amplitude{  constexpr static auto name=my_static_string("Amplitude");};

struct Frecuency{ constexpr static auto name=my_static_string("Frecuency");};


struct phase{ constexpr static auto name=my_static_string("phase");};

struct tau{ constexpr static auto name=my_static_string("tau");};




#endif // ANTENA_VARIABLES_H
