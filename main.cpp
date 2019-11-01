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


template<class D> struct Distribution_of{constexpr static auto className=D::className+my_static_string("_Distribution");};
template<class D> struct Values_of{constexpr static auto className=D::className+my_static_string("_Values");};

struct prior{constexpr static auto className=my_static_string("prior");};

struct model{constexpr static auto className=my_static_string("model");};






int main()
{
  auto myprior_values=vector_space{
      x_i(mean<baseline>{},v(0.0,V_u{})),
      x_i(stddev<baseline>{},v(1e-6,V_u{})),
      x_i(mean<drift>{},v(0.0,V_u{}/ps_u{})),
      x_i(stddev<drift>{},v(1e-7,V_u{}/ps_u{})),
      x_i(mean<Log10_t<stddev<signal>>>{},logv<double,V_u>(-6.0,{1,V_u{}})),
      x_i(stddev<Log10_t<stddev<signal>>>{},logv<double>(1.0)),
      x_i(mean<ind<Log10_t<Amplitude>,0>>{},logv<double,V_u>(-5.0,{1,V_u{}})),
      x_i(stddev<ind<Log10_t<Amplitude>,0>>{},logv<double>(1.0)),
      x_i(mean<ind<Log10_t<Frecuency>,0>>{},logv<double,GHz_u>(std::log10(9),{1,GHz_u{}})),
      x_i(stddev<ind<Log10_t<Frecuency>,0>>{},logv<double>(0.2)),
      x_i(mean<ind<phase,0>>{},v(PI/4)),
      x_i(stddev<ind<phase,0>>{},v(PI/4)),
      x_i(mean<ind<Log10_t<tau>,0>>{},logv<double,ps_u>(std::log10(1000),{1,ps_u{}})),
      x_i(stddev<ind<Log10_t<tau>,0>>{},logv<double>(1)),
      x_i(mean<ind<Log10_t<Amplitude>,1>>{},logv<double,V_u>(-5.0,{1,V_u{}})),
      x_i(stddev<ind<Log10_t<Amplitude>,1>>{},logv<double>(1.0)),
      x_i(mean<ind<Log10_t<Frecuency>,1>>{},logv<double,GHz_u>(std::log10(9),{1,GHz_u{}})),
      x_i(stddev<ind<Log10_t<Frecuency>,1>>{},logv<double>(0.2)),
      x_i(mean<ind<phase,1>>{},v(PI/4)),
      x_i(stddev<ind<phase,1>>{},v(PI/4)),
      x_i(mean<ind<Log10_t<tau>,1>>{},logv<double,ps_u>(std::log10(1000),{1,ps_u{}})),
      x_i(stddev<ind<Log10_t<tau>,1>>{},logv<double>(1))};

  auto myprior_transf=quimulun{
      F(stddev<signal>{},Log10_rev{},Log10_t<stddev<signal>>{}),

      F(ind<Amplitude,0>{},Log10_rev{},ind<Log10_t<Amplitude>,0>{}),
      F(ind<Frecuency,0>{},Log10_rev{},ind<Log10_t<Frecuency>,0>{}),
      F(ind<tau,0>{},Log10_rev{},ind<Log10_t<tau>,0>{}),

      F(ind<Amplitude,1>{},Log10_rev{},ind<Log10_t<Amplitude>,1>{}),
      F(ind<Frecuency,1>{},Log10_rev{},ind<Log10_t<Frecuency>,1>{}),
      F(ind<tau,1>{},Log10_rev{},ind<Log10_t<tau>,1>{})

  };

  auto myprior_dist=quimulun{
      D(baseline{},Normal_Distribution{},mean<baseline>{},stddev<baseline>{}),
        D(drift{},Normal_Distribution{},mean<drift>{},stddev<drift>{}),
        D(Log10_t<stddev<signal>>{},Normal_Distribution{},mean<Log10_t<stddev<signal>>>{},stddev<Log10_t<stddev<signal>>>{}),
        D(ind<Log10_t<Amplitude>,0>{},Normal_Distribution{},mean<ind<Log10_t<Amplitude>,0>>{},stddev<ind<Log10_t<Amplitude>,0>>{}),
        D(ind<Log10_t<Frecuency>,0>{},Normal_Distribution{},mean<ind<Log10_t<Frecuency>,0>>{},stddev<ind<Log10_t<Frecuency>,0>>{}),
        D(ind<phase,0>{},Normal_Distribution{},mean<ind<phase,0>>{},stddev<ind<phase,0>>{}),
        D(ind<Log10_t<tau>,0>{},Normal_Distribution{},mean<ind<Log10_t<tau>,0>>{},stddev<ind<Log10_t<tau>,0>>{}),
        D(ind<Log10_t<Amplitude>,1>{},Normal_Distribution{},mean<ind<Log10_t<Amplitude>,1>>{},stddev<ind<Log10_t<Amplitude>,1>>{}),
        D(ind<Log10_t<Frecuency>,1>{},Normal_Distribution{},mean<ind<Log10_t<Frecuency>,1>>{},stddev<ind<Log10_t<Frecuency>,1>>{}),
        D(ind<phase,1>{},Normal_Distribution{},mean<ind<phase,1>>{},stddev<ind<phase,1>>{}),
      D(ind<Log10_t<tau>,1>{},Normal_Distribution{},mean<ind<Log10_t<tau>,1>>{},stddev<ind<Log10_t<tau>,1>>{})};

  auto mymodel=
      quimulun{
          D(signal{},Normal_Distribution{},mean<signal>{},stddev<signal>{}),
          F(mean<signal>{},
              [](auto t, auto baseline_,auto drift_,
                 auto A0_, auto f0_, auto ph0_, auto tau0_,
                 auto A1_, auto f1_, auto ph1_, auto tau1_ )
              {

                    return baseline_+(drift_*t)+
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

  auto totalmodel=mymodel+myprior_dist+myprior_transf+myprior_values;

  auto data_prior=data;

  std::random_device rd;
  auto initseed = 0;
  //rd();

  std::mt19937 mt(initseed);

  auto data2=data;
  auto s=sample(totalmodel,std::move(data2),mt);

  std::cerr<<s;

/*  auto fname2="antena_data_2.txt";
  std::ofstream of(fname2);
  to_DataFrame(of,data);

  auto data_2=decltype (data){};
  of.close();
  std::ifstream if2(fname2);
  from_DataFrame(if2,data_2);

  auto fname3="antena_data_3.txt";
  std::ofstream of3(fname3);
  to_DataTable(of3,data);
*/



//  std::cerr<<data;



  return 0;
}
