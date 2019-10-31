#include <iostream>
#include <static_string.h>
#include <qm_unit.h>
#include <mytypetraits.h>
#include <qm_tensor_model.h>
#include <fstream>
struct ps{constexpr static auto  className=my_static_string("ps");};
struct V{constexpr static auto  className=my_static_string("V");};
struct GHz{constexpr static auto  className=my_static_string("GHz");};



typedef p<u<ps,1>> ps_u;
typedef p<u<V,1>> V_u;
typedef p<u<GHz,1>> GHz_u;



struct delay
{
  typedef double T;
  typedef ps_u unit;
  constexpr static auto  className=my_static_string("delay");
};
struct signal
{
  typedef double T;
  typedef V_u unit;
  constexpr static auto  className=my_static_string("signal");
};


template <class ei, std::size_t I>
struct rep
{
  using T=typename ei::T;
  using unit=typename ei::unit;
  constexpr static auto  className=ei::className+my_static_string("_")+to_static_string<I>();
};

template <class ei, std::size_t I>
struct ind
{
  using T=typename ei::T;
  using unit=typename ei::unit;
  constexpr static auto  className=ei::className+my_static_string("_")+to_static_string<I>();
};



struct baseline{using T=double; using unit=V_u; constexpr static auto className=my_static_string("baseline");};

struct drift{using T=double; using unit=decltype (V_u{}/ps_u{}); constexpr static auto className=my_static_string("baseline");};

struct Amplitude{using T=double; using unit=V_u; constexpr static auto className=my_static_string("Amplitude");};

struct Frecuency{using T=double; using unit=GHz_u; constexpr static auto className=my_static_string("Frecuency");};


struct phase{using T=double; using unit=dimension_less; constexpr static auto className=my_static_string("phase");};

struct tau{using T=double; using unit=ps_u; constexpr static auto className=my_static_string("tau");};


constexpr const auto GHz_f=v(1E9,p<u<GHz,-1>,u<s,-1>>{});
constexpr const auto ps_f=v(1E-12,p<u<s,1>,u<ps,-1>>{});


struct prior{constexpr static auto className=my_static_string("prior");};

struct model{constexpr static auto className=my_static_string("model");};




int main()
{


  auto myprior=quimulun{prior{},
    D(baseline{},Normal_Distribution{},mean<baseline>{},stddev<baseline>{}),
        D(drift{},Normal_Distribution{},mean<baseline>{},stddev<baseline>{}),
        D(ind<Amplitude,0>{},Normal_Distribution{},mean<ind<Amplitude,0>>{},stddev<ind<Amplitude,0>>{}),
        D(ind<Frecuency,0>{},Normal_Distribution{},mean<ind<Frecuency,0>>{},stddev<ind<Frecuency,0>>{}),
        D(ind<phase,0>{},Normal_Distribution{},mean<ind<phase,0>>{},stddev<ind<phase,0>>{}),
        D(ind<tau,0>{},Normal_Distribution{},mean<ind<tau,0>>{},stddev<ind<tau,0>>{}),
        D(ind<Amplitude,1>{},Normal_Distribution{},mean<ind<Amplitude,1>>{},stddev<ind<Amplitude,1>>{}),
        D(ind<Frecuency,1>{},Normal_Distribution{},mean<ind<Frecuency,1>>{},stddev<ind<Frecuency,1>>{}),
        D(ind<phase,1>{},Normal_Distribution{},mean<ind<phase,1>>{},stddev<ind<phase,1>>{}),
        D(ind<tau,1>{},Normal_Distribution{},mean<ind<tau,1>>{},stddev<ind<tau,1>>{})};

  auto mymodel=
      quimulun{
          model{},
          D(signal{},Normal_Distribution{},mean<signal>{},stddev<signal>{}),
          F(mean<signal>{},
              [](auto t, auto baseline_,auto drift_,
                 auto A0_, auto f0_, auto ph0_, auto tau0_,
                 auto A1_, auto f1_, auto ph1_, auto tau1_ ){
                return baseline_+drift_*t+
                       A0_*exp(-t/tau0_)*cos(2*PI*f0_*GHz_f*t*ps_f+ph0_)+
                       A1_*exp(-t/tau1_)*cos(2*PI*f1_*GHz_f*t*ps_f+ph1_);},
              delay{},baseline{},drift{},
              ind<Amplitude,0>{},ind<Frecuency,0>{},ind<phase,0>{},ind<tau,0>{},
              ind<Amplitude,1>{},ind<Frecuency,1>{},ind<phase,1>{},ind<tau,1>{})
  };

  auto data=vector_space{
      x_i(delay{}, vec<delay>{}),x_i(signal{},vec<delay>{})};

  auto fname="antena_data_1.txt";
  std::ifstream fi(fname);
  from_DataTable(fi,data);

  auto fname2="antena_data_2.txt";
  std::ofstream of(fname2);
  to_DataFrame(of,data);

  auto data_2=decltype (data){};
  of.close();
  std::ifstream if2(fname2);
  from_DataFrame(if2,data_2);

  auto fname3="antena_data_3.txt";
  std::ofstream of3(fname3);
  to_DataTable(of3,data);



  std::cerr<<data;



  return 0;
}
