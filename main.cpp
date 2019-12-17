#include <iostream>
#include <static_string.h>
#include <qm_unit.h>
#include <mytypetraits.h>
#include <qm_tensor_model.h>
#include <fstream>
#include <my_tests.h>
#include <qm_Metropolis_Parallel.h>
#include <qm_data_frame.h>
struct ps{constexpr static auto  className=my_static_string("ps");};
struct V{constexpr static auto  className=my_static_string("V");};
struct GHz{constexpr static auto  className=my_static_string("GHz");};


typedef p<u<ps,1>> ps_u;
typedef p<u<V,1>> V_u;
typedef p<u<GHz,1>> GHz_u;

constexpr const auto GHz_f=v(1E9,p_t<u<GHz,-1>,u<s,-1>>{});
constexpr const auto ps_f=v(1E-12,p_t<u<s,1>,u<ps,-1>>{});

constexpr const auto Gsad=v(1,GHz_u{});
//using r=typename p_t<u<GHz,-1>,u<s,-1>>::uds;



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



template <class ei, std::size_t... I>
struct ind
{
  using T=typename ei::T;
  using unit=typename ei::unit;
  constexpr static auto  className=(ei::className+...+(my_static_string("_")+to_static_string<I>()));
};



struct baseline{using T=double; using unit=V_u; constexpr static auto className=my_static_string("baseline");};

struct drift{using T=double; using unit=decltype (V_u{}/ps_u{}); constexpr static auto className=my_static_string("drift");};

struct Amplitude{using T=double; using unit=V_u; constexpr static auto className=my_static_string("Amplitude");};

struct Frecuency{using T=double; using unit=GHz_u; constexpr static auto className=my_static_string("Frecuency");};


struct phase{using T=double; using unit=dimension_less; constexpr static auto className=my_static_string("phase");};

struct tau{using T=double; using unit=ps_u; constexpr static auto className=my_static_string("tau");};




template<class Id> struct Sorting_index{constexpr static auto className=my_static_string("Sorting_index")+Id::className;};
template<class Id> struct sorted{constexpr static auto className=Id::className+my_static_string("_sorted");};

auto my_common_prior_values(){
  return vector_space{
      x_i(mean<baseline>{},v(0.0,V_u{})),
      x_i(stddev<baseline>{},v(1e-6,V_u{})),
      x_i(mean<drift>{},v(0.0,V_u{}/ps_u{})),
      x_i(stddev<drift>{},v(1e-7,V_u{}/ps_u{})),
      x_i(mean<Log10_t<stddev<signal>>>{},logv(-6.0,V_u{})),
      x_i(stddev<Log10_t<stddev<signal>>>{},v(1.0)),
      x_i(mean<Log10_t<Amplitude>>{},logv(-5.0,V_u{})),
      x_i(stddev<Log10_t<Amplitude>>{},v(1.0)),
      x_i(mean<Log10_t<ind<Frecuency,0>>>{},logv(std::log10(0.0001),GHz_u{})),
      x_i(stddev<Log10_t<ind<Frecuency,0>>>{},v(0.2)),
      x_i(mean<Log10_t<Frecuency>>{},logv(std::log10(9),GHz_u{})),
      x_i(stddev<Log10_t<Frecuency>>{},v(0.2)),
      x_i(mean<phase>{},v(PI/4)),
      x_i(stddev<phase>{},v(PI/4)),
      x_i(mean<Log10_t<tau>>{},logv(std::log10(1000),ps_u{})),
      x_i(stddev<Log10_t<tau>>{},v(1.0))};
}



template<std::size_t J>
auto myprior_dist_stddev()
{
  return quimulun{
      D(ind<baseline,J>{},Normal_Distribution{},mean<baseline>{},stddev<baseline>{}),
      D(ind<drift,J>{},Normal_Distribution{},mean<drift>{},stddev<drift>{}),
      D(Log10_t<stddev<ind<signal,J>>>{},Normal_Distribution{},mean<Log10_t<stddev<signal>>>{},stddev<Log10_t<stddev<signal>>>{})

  };
}

template<std::size_t I>
auto myprior_dist_frequency()
{
  if constexpr (I>0)
    return quimulun{
        D(ind<Log10_t<Frecuency>,I>{},Normal_Distribution{},mean<Log10_t<Frecuency>>{},stddev<Log10_t<Frecuency>>{})};
  else  return quimulun{
        D(ind<Log10_t<Frecuency>,I>{},Normal_Distribution{},mean<Log10_t<ind<Frecuency,0>>>{},stddev<Log10_t<Frecuency>>{})};

}
template<std::size_t I, std::size_t J>
auto myprior_dist_frequency()
{
  return quimulun{
      D(ind<Log10_t<Frecuency>,I,J>{},Normal_Distribution{},mean<Log10_t<Frecuency>>{},stddev<Log10_t<Frecuency>>{})};
}


template<std::size_t I, std::size_t J>
auto myprior_dist_amplitude()
{
  return quimulun{
      D(ind<Log10_t<Amplitude>,I,J>{},Normal_Distribution{},mean<Log10_t<Amplitude>>{},stddev<Log10_t<Amplitude>>{}),
      D(ind<phase,I,J>{},Normal_Distribution{},mean<phase>{},stddev<phase>{}),
      D(ind<Log10_t<tau>,I,J>{},Normal_Distribution{},mean<Log10_t<tau>>{},stddev<Log10_t<tau>>{})};
}

template <std::size_t J, std::size_t...I0>
auto myprior_dist_index(std::index_sequence<I0...> ,std::index_sequence<> )
{
  return myprior_dist_stddev<J>()
         +(myprior_dist_amplitude<I0,J>()+...);
}

template <std::size_t J, std::size_t...I0, std::size_t...I1>
auto myprior_dist_index(std::index_sequence<I0...> ,std::index_sequence<I1...> )
{
  return myprior_dist_stddev<J>()+(myprior_dist_frequency<I1,J>()+...)
         +(myprior_dist_amplitude<I0,J>()+...)+(myprior_dist_amplitude<I1,J>()+...);
}
template <std::size_t... J, std::size_t...I0, std::size_t...I1>
auto myprior_dist_index(std::index_sequence<J...>,std::index_sequence<I0...> i0,std::index_sequence<I1...> i1)
{
  return (myprior_dist_index<J>(i0,i1)+...)+(myprior_dist_frequency<I0>()+...);
}

template<std::size_t... I>
auto myprior_sort_frequency(std::index_sequence<I...>)
{ return quimulun{
      F(Sorting_index<Frecuency>{},[](auto const & x...){
            auto f=std::array{std::pair{x.value(),I}...};
            std::sort(f.begin(),f.end(),[](auto p1, auto p2){return p1.first<p2.first;});
            auto out=std::apply([](auto... p){return std::array{p.second...};},f);
            return v(out,dimension_less{});
          },ind<Log10_t<Frecuency>,I>{}...),
      F(sorted<ind<Log10_t<Frecuency>,I>>{},[](auto const& sorting_index, auto const&... freq){
            auto f=std::array{freq...};
            return f[sorting_index.value()[I]];
          },Sorting_index<Frecuency>{},Cs<ind<Log10_t<Frecuency>,I>...>{})...,
      F(sorted<ind<Frecuency,I>>{},Log10_rev{},sorted<ind<Log10_t<Amplitude>,I>>{})...,

      };
}
template<std::size_t... I, std::size_t J>
auto myprior_sort_Amplitude(std::index_sequence<I...>, std::integral_constant<std::size_t, J>)
{ return quimulun{
      F(sorted<ind<Log10_t<Amplitude>,I,J>>{},[](auto const& sorting_index, auto const&... x){
            auto f=std::array{x...};
            return f[sorting_index.value()[I]];
          },Sorting_index<Frecuency>{},Cs<ind<Log10_t<Amplitude>,I,J>...>{})...,
      F(sorted<ind<Amplitude,I,J>>{},Log10_rev{},sorted<ind<Log10_t<Amplitude>,I,J>>{})...,
      F(sorted<ind<Log10_t<tau>,I,J>>{},[](auto const& sorting_index, auto const&... x){
            auto f=std::array{x...};
            return f[sorting_index.value()[I]];
          },Sorting_index<tau>{},Cs<ind<Log10_t<Amplitude>,I,J>...>{})...,
      F(sorted<ind<tau,I,J>>{},Log10_rev{},sorted<ind<Log10_t<tau>,I,J>>{})...,
      F(sorted<ind<phase,I,J>>{},[](auto const& sorting_index, auto const&... x){
            auto f=std::array{x...};
            return f[sorting_index.value()[I]];
          },Sorting_index<tau>{},Cs<ind<phase,I,J>...>{})...,

      };


}

template<std::size_t J>
auto myprior_transf_stddev(){return quimulun{
      F(stddev<ind<signal,J>>{},Log10_rev{},Log10_t<stddev<ind<signal,J>>>{})
  };
}

    template<std::size_t J,std::size_t... I>
    auto myprior_transf_Amplitude(std::index_sequence<I...>)
{
  return quimulun{
      F(ind<Amplitude,I,J>{},Log10_rev{},ind<Log10_t<Amplitude>,I,J>{})...,
      F(ind<tau,I,J>{},Log10_rev{},ind<Log10_t<tau>,I,J>{})...,
      };
}

template<std::size_t I>
auto myprior_transf_frequency()
{
  return quimulun{
      F(ind<Frecuency,I>{},Log10_rev{},ind<Log10_t<Frecuency>,I>{})
  };
}
template<std::size_t I, std::size_t J>
auto myprior_transf_frequency()
{
  return quimulun{
      F(ind<Frecuency,I,J>{},Log10_rev{},ind<Log10_t<Frecuency>,I,J>{})
  };
}

template<std::size_t I,std::size_t... Is>
auto myprior_transf_Frequency_sum(std::integral_constant<std::size_t,I>,std::index_sequence<Is...>)
{
  return quimulun{
      F(ind<Frecuency,I>{},
          [](auto const& ...f){
            using std::pow;
            return (pow(10.0,f)+...);}

          ,ind<Log10_t<Frecuency>,I-Is>{}...)
  };
}

template<std::size_t I, std::size_t... Is>
auto myprior_transf_Frequency_sum(std::index_sequence<I,Is...>)
{
  return (myprior_transf_Frequency_sum(std::integral_constant<std::size_t,I>{},std::index_sequence<0>{})+
          ...+myprior_transf_Frequency_sum(std::integral_constant<std::size_t,Is>{},std::make_index_sequence<Is-I>{}));
}

template<std::size_t J,std::size_t I,std::size_t... Is>
auto myprior_transf_Frequency_sum(std::index_sequence<Is...>)
{
  return quimulun{
      F(ind<Frecuency,I,J>{},
          [](auto const& ...f){
            using std::pow;
            return (pow(10.0,f)+...);}

          ,ind<Log10_t<Frecuency>,I-Is,J>{}...)
  };
}

template<std::size_t J, std::size_t I,std::size_t... Is>
auto myprior_transf_Frequency_sum(std::index_sequence<I,Is...>)
{
  return (myprior_transf_Frequency_sum<J,I>(std::index_sequence<0>{})+
          ...+myprior_transf_Frequency_sum<J,Is>(std::make_index_sequence<Is-I>{}));
}


template <std::size_t J, std::size_t...I0>
auto myprior_transf_index(std::index_sequence<I0...> i0,std::index_sequence<> )
{
  return myprior_transf_stddev<J>()+
         myprior_transf_Amplitude<J>(i0);
}


template <std::size_t J, std::size_t...I0, std::size_t I00,std::size_t...I1>
auto myprior_transf_index(std::index_sequence<I0...> i0,std::index_sequence<I00,I1...> i1)
{
  return myprior_transf_stddev<J>()+(myprior_transf_frequency<I00,J>()+...+myprior_transf_frequency<I1,J>())+
         myprior_transf_Amplitude<J>(i0)+myprior_transf_Amplitude<J>(i1);
}
template <std::size_t... J, std::size_t...I0, std::size_t...I1>
auto myprior_transf_index(std::index_sequence<J...>,std::index_sequence<I0...> i0,std::index_sequence<I1...> i1)
{
  return (myprior_transf_index<J>(i0,i1)+...)+(myprior_transf_frequency<I0>()+...);
}

template <std::size_t J, std::size_t...I0>
auto mymodel_index(std::index_sequence<I0...>,std::index_sequence<>){return
      quimulun{
          D(ind<signal,J>{},Normal_Distribution{},mean<ind<signal,J>>{},stddev<ind<signal,J>>{}),
          F(mean<ind<signal,J>>{},
              [](auto t, auto baseline_,auto drift_,
                 auto const& A0_, auto const& f0_,  auto const& ph0_, auto const& tau0_)

              {
                auto y0=baseline_+(drift_*t);
                auto y1=((std::get<I0>(A0_.value())()*exp(-t/std::get<I0>(tau0_.value())())*
                            cos(2*PI*std::get<I0>(f0_.value())()*GHz_f*t*ps_f+std::get<I0>(ph0_.value())()))+...);

                return y0+y1;

                //return baseline_;
              },
              ind<delay,J>{},ind<baseline,J>{},ind<drift,J>{},
              std::tuple<ind<Amplitude,I0,J>...>{},
              std::tuple<ind<Frecuency,I0>...>{},
              std::tuple<ind<phase,I0,J>...>{},
              std::tuple<ind<tau,I0,J>...>{})
      };


}


template <std::size_t J, std::size_t...I0, std::size_t I00, std::size_t...I1>
auto mymodel_index(std::index_sequence<I0...>,std::index_sequence<I00,I1...>){return
      quimulun{
          D(ind<signal,J>{},Normal_Distribution{},mean<ind<signal,J>>{},stddev<ind<signal,J>>{}),
          F(mean<ind<signal,J>>{},
              [](auto t, auto baseline_,auto drift_,
                 auto const& A0_, auto const& A1_,auto const& f0_,auto const& f1_,  auto const& ph0_, auto const& ph1_,auto const& tau0_,
                 auto const& tau1_ )

              {
                auto y0=baseline_+(drift_*t);
                auto y1=((std::get<I0>(A0_.value())()*exp(-t/std::get<I0>(tau0_.value())())*
                            cos(2*PI*std::get<I0>(f0_.value())()*GHz_f*t*ps_f+std::get<I0>(ph0_.value())()))+...);
                auto y2=((std::get<0>(A1_.value())()*exp(-t/std::get<0>(tau1_.value())())*
                            cos(2*PI*std::get<0>(f1_.value())()*GHz_f*t*ps_f+std::get<0>(ph1_.value())())));
                if constexpr (sizeof... (I1)==0)

                  return y0+y1+y2;
                else
                {
                  auto y3=((std::get<I1-I00>(A1_.value())()*exp(-t/std::get<I1-I00>(tau1_.value())())*
                              cos(2*PI*std::get<I1-I00>(f1_.value())()*GHz_f*t*ps_f+std::get<I1-I00>(ph1_.value())()))+...);
                  return y0+y1+y2+y3;

                }
                //return baseline_;
              },
              ind<delay,J>{},ind<baseline,J>{},ind<drift,J>{},
              std::tuple<ind<Amplitude,I0,J>...>{},
              std::tuple<ind<Amplitude,I00,J>,ind<Amplitude,I1,J>...>{},
              std::tuple<ind<Frecuency,I0>...>{},
              std::tuple<ind<Frecuency,I00,J>,ind<Frecuency,I1,J>...>{},
              std::tuple<ind<phase,I0,J>...>{},
              std::tuple<ind<phase,I00,J>,ind<phase,I1,J>...>{},
              std::tuple<ind<tau,I0,J>...>{},
              std::tuple<ind<tau,I00,J>,ind<tau,I1,J>...>{})
      };


}

template <std::size_t... J, std::size_t...I0, std::size_t...I1>
auto mymodel_index(std::index_sequence<J...>,std::index_sequence<I0...> i0,std::index_sequence<I1...> i1)
{
  return (mymodel_index<J>(i0,i1)+...);
}


template <std::size_t... J, std::size_t...I0, std::size_t...I1>
auto totalmodel_index(std::index_sequence<J...> j,std::index_sequence<I0...> i0,std::index_sequence<I1...> i1)
{
  return mymodel_index(j,i0,i1)+myprior_dist_index(j,i0,i1)
         +myprior_transf_index(j,i0,i1);
  //+my_common_prior_values();
}

template <std::size_t... J, std::size_t...I0>
auto totalmodel_index(std::index_sequence<J...> j,std::index_sequence<I0...> i0,std::index_sequence<> i1)
{
  return mymodel_index(j,i0,i1)+myprior_dist_index(j,i0,i1)
         +myprior_transf_index(j,i0,i1);
  //+my_common_prior_values();
}


template<std::size_t J>
auto get_data_index(){return vector_space{x_i(ind<delay,J>{}, vec<ind<delay,J>>{}),
                                            x_i(ind<signal,J>{},vec<ind<delay,J>>{})};}

template<std::size_t N>
auto totalmodel_41_(){ return totalmodel_index(std::index_sequence<N>{},std::index_sequence<0,1,2>{},std::index_sequence<>{})+my_common_prior_values();}


int main()
{
  auto myprior_values=vector_space{
      x_i(mean<baseline>{},v(0.0,V_u{})),
      x_i(stddev<baseline>{},v(1e-6,V_u{})),
      x_i(mean<drift>{},v(0.0,V_u{}/ps_u{})),
      x_i(stddev<drift>{},v(1e-7,V_u{}/ps_u{})),
      x_i(mean<Log10_t<stddev<signal>>>{},logv(-6.0,V_u{})),
      x_i(stddev<Log10_t<stddev<signal>>>{},v(1.0)),

      x_i(mean<ind<Log10_t<Amplitude>,0>>{},logv(-5.0,V_u{})),
      x_i(stddev<ind<Log10_t<Amplitude>,0>>{},v(1.0)),
      x_i(mean<ind<Log10_t<Frecuency>,0>>{},logv(std::log10(9),GHz_u{})),
      x_i(stddev<ind<Log10_t<Frecuency>,0>>{},v(0.2)),
      x_i(mean<ind<phase,0>>{},v(PI/4)),
      x_i(stddev<ind<phase,0>>{},v(PI/4)),
      x_i(mean<ind<Log10_t<tau>,0>>{},logv(std::log10(1000),ps_u{})),
      x_i(stddev<ind<Log10_t<tau>,0>>{},v(1.0)),

      x_i(mean<ind<Log10_t<Amplitude>,1>>{},logv(-5.0,V_u{})),
      x_i(stddev<ind<Log10_t<Amplitude>,1>>{},v(1.0)),
      x_i(mean<ind<Log10_t<Frecuency>,1>>{},logv(std::log10(9),GHz_u{})),
      x_i(stddev<ind<Log10_t<Frecuency>,1>>{},v(0.2)),
      x_i(mean<ind<phase,1>>{},v(PI/4)),
      x_i(stddev<ind<phase,1>>{},v(PI/4)),
      x_i(mean<ind<Log10_t<tau>,1>>{},logv(std::log10(1000),ps_u{})),
      x_i(stddev<ind<Log10_t<tau>,1>>{},v(1.0))


  };

      auto myprior_values_3=vector_space{
          x_i(mean<baseline>{},v(0.0,V_u{})),
          x_i(stddev<baseline>{},v(1e-6,V_u{})),
          x_i(mean<drift>{},v(0.0,V_u{}/ps_u{})),
          x_i(stddev<drift>{},v(1e-7,V_u{}/ps_u{})),
          x_i(mean<Log10_t<stddev<signal>>>{},logv(-6.0,V_u{})),
          x_i(stddev<Log10_t<stddev<signal>>>{},v(1.0)),

          x_i(mean<Log10_t<Amplitude>>{},logv(-5.0,V_u{})),
          x_i(stddev<Log10_t<Amplitude>>{},v(1.0)),
          x_i(mean<Log10_t<Frecuency>>{},logv(std::log10(9),GHz_u{})),
          x_i(mean<Log10_t<ind<Frecuency,0>>>{},logv(std::log10(0.0001),GHz_u{})),
          x_i(stddev<Log10_t<Frecuency>>{},v(1.0)),
          x_i(mean<phase>{},v(PI/4)),
          x_i(stddev<phase>{},v(PI/4)),
          x_i(mean<Log10_t<tau>>{},logv(std::log10(1000),ps_u{})),
          x_i(stddev<Log10_t<tau>>{},v(1.0)),



          };

  auto myprior_transf=quimulun{
      F(stddev<signal>{},Log10_rev{},Log10_t<stddev<signal>>{}),

      F(ind<Amplitude,0>{},Log10_rev{},ind<Log10_t<Amplitude>,0>{}),
      F(ind<Frecuency,0>{},Log10_rev{},ind<Log10_t<Frecuency>,0>{}),
      F(ind<tau,0>{},Log10_rev{},ind<Log10_t<tau>,0>{}),

      F(ind<Amplitude,1>{},Log10_rev{},ind<Log10_t<Amplitude>,1>{}),
      F(ind<Frecuency,1>{},Log10_rev{},ind<Log10_t<Frecuency>,1>{}),
      F(ind<tau,1>{},Log10_rev{},ind<Log10_t<tau>,1>{})

  };

  auto myprior_transf_3=quimulun{
      F(stddev<signal>{},Log10_rev{},Log10_t<stddev<signal>>{}),

      F(ind<Amplitude,0>{},Log10_rev{},ind<Log10_t<Amplitude>,0>{}),
      F(ind<Frecuency,0>{},Log10_rev{},ind<Log10_t<Frecuency>,0>{}),
      F(ind<tau,0>{},Log10_rev{},ind<Log10_t<tau>,0>{}),

      F(ind<Amplitude,1>{},Log10_rev{},ind<Log10_t<Amplitude>,1>{}),
      F(ind<Frecuency,1>{},Log10_rev{},ind<Log10_t<Frecuency>,1>{}),
      F(ind<tau,1>{},Log10_rev{},ind<Log10_t<tau>,1>{}),

      F(ind<Amplitude,2>{},Log10_rev{},ind<Log10_t<Amplitude>,1>{}),
      F(ind<Frecuency,2>{},Log10_rev{},ind<Log10_t<Frecuency>,1>{}),
      F(ind<tau,2>{},Log10_rev{},ind<Log10_t<tau>,1>{})

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



      auto myprior_dist_3=quimulun{
          D(baseline{},Normal_Distribution{},mean<baseline>{},stddev<baseline>{}),
          D(drift{},Normal_Distribution{},mean<drift>{},stddev<drift>{}),
          D(Log10_t<stddev<signal>>{},Normal_Distribution{},mean<Log10_t<stddev<signal>>>{},stddev<Log10_t<stddev<signal>>>{}),

          D(ind<Log10_t<Amplitude>,0>{},Normal_Distribution{},mean<Log10_t<Amplitude>>{},stddev<Log10_t<Amplitude>>{}),
          D(ind<Log10_t<Frecuency>,0>{},Normal_Distribution{},mean<Log10_t<ind<Frecuency,0>>>{},stddev<Log10_t<Frecuency>>{}),
          D(ind<phase,0>{},Normal_Distribution{},mean<phase>{},stddev<phase>{}),
          D(ind<Log10_t<tau>,0>{},Normal_Distribution{},mean<Log10_t<tau>>{},stddev<Log10_t<tau>>{}),

          D(ind<Log10_t<Amplitude>,1>{},Normal_Distribution{},mean<Log10_t<Amplitude>>{},stddev<Log10_t<Amplitude>>{}),
          D(ind<Log10_t<Frecuency>,1>{},Normal_Distribution{},mean<Log10_t<Frecuency>>{},stddev<Log10_t<Frecuency>>{}),
          D(ind<phase,1>{},Normal_Distribution{},mean<phase>{},stddev<phase>{}),
          D(ind<Log10_t<tau>,1>{},Normal_Distribution{},mean<Log10_t<tau>>{},stddev<Log10_t<tau>>{}),

          D(ind<Log10_t<Amplitude>,2>{},Normal_Distribution{},mean<Log10_t<Amplitude>>{},stddev<Log10_t<Amplitude>>{}),
          D(ind<Log10_t<Frecuency>,2>{},Normal_Distribution{},mean<Log10_t<Frecuency>>{},stddev<Log10_t<Frecuency>>{}),
          D(ind<phase,2>{},Normal_Distribution{},mean<phase>{},stddev<phase>{}),
          D(ind<Log10_t<tau>,2>{},Normal_Distribution{},mean<Log10_t<tau>>{},stddev<Log10_t<tau>>{}),
          };

  auto myprior_dist_5=quimulun{
      D(baseline{},Normal_Distribution{},mean<baseline>{},stddev<baseline>{}),
      D(drift{},Normal_Distribution{},mean<drift>{},stddev<drift>{}),
      D(Log10_t<stddev<signal>>{},Normal_Distribution{},mean<Log10_t<stddev<signal>>>{},stddev<Log10_t<stddev<signal>>>{}),

      D(ind<Log10_t<Amplitude>,0>{},Normal_Distribution{},mean<ind<Log10_t<Amplitude>,0>>{},stddev<ind<Log10_t<Amplitude>,0>>{}),
      D(ind<Log10_t<Frecuency>,0>{},Normal_Distribution{},mean<ind<Log10_t<Frecuency>,0>>{},stddev<ind<Log10_t<Frecuency>,0>>{}),
      D(ind<phase,0>{},Normal_Distribution{},mean<ind<phase,0>>{},stddev<ind<phase,0>>{}),
      D(ind<Log10_t<tau>,0>{},Normal_Distribution{},mean<ind<Log10_t<tau>,0>>{},stddev<ind<Log10_t<tau>,0>>{}),

      D(ind<Log10_t<Frecuency>,1>{},Normal_Distribution{},mean<Log10_t<Frecuency>>{},stddev<ind<Log10_t<Frecuency>,1>>{}),
      D(ind<Log10_t<Amplitude>,1>{},Normal_Distribution{},mean<Log10_t<Amplitude>>{},stddev<ind<Log10_t<Amplitude>,1>>{}),
      D(ind<phase,1>{},Normal_Distribution{},mean<phase>{},stddev<ind<phase>>{}),
      D(ind<Log10_t<tau>,1>{},Normal_Distribution{},mean<Log10_t<tau>>{},stddev<ind<Log10_t<tau>,1>>{}),

      D(ind<Log10_t<Amplitude>,2>{},Normal_Distribution{},mean<ind<Log10_t<Amplitude>,1>>{},stddev<ind<Log10_t<Amplitude>,1>>{}),
      D(ind<Log10_t<Frecuency>,2>{},Normal_Distribution{},mean<ind<Log10_t<Frecuency>,1>>{},stddev<ind<Log10_t<Frecuency>,1>>{}),
      D(ind<phase,2>{},Normal_Distribution{},mean<ind<phase,1>>{},stddev<ind<phase,1>>{}),
      D(ind<Log10_t<tau>,2>{},Normal_Distribution{},mean<ind<Log10_t<tau>,1>>{},stddev<ind<Log10_t<tau>,1>>{})};




  auto mymodel=
      quimulun{
          D(signal{},Normal_Distribution{},mean<signal>{},stddev<signal>{}),
          F(mean<signal>{},
              [](auto t, auto baseline_,auto drift_,
                 auto A0_, auto f0_, auto ph0_, auto tau0_,
                 auto A1_, auto f1_, auto ph1_, auto tau1_ )
              {
//                using a=typename decltype(t)::t;
//                using b=typename decltype(baseline_)::baseline;
//                using c=typename decltype(drift_)::drift;
 //              using d=typename decltype(A0_)::A0;
//                using e=typename decltype(tau0_)::tau;

                    return baseline_+(drift_*t)+
                       A0_*exp(
                                 -t/tau0_
                                 )*
                           cos(2*PI*f0_*GHz_f*t*ps_f
                               +ph0_)+
                       A1_*exp(-t/tau1_)*cos(2*PI*f1_*GHz_f*t*ps_f+ph1_);},
              delay{},baseline{},drift{},
              ind<Amplitude,0>{},ind<Frecuency,0>{},ind<phase,0>{},ind<tau,0>{},
              ind<Amplitude,1>{},ind<Frecuency,1>{},ind<phase,1>{},ind<tau,1>{})
      };

  auto mymodel_3=
      quimulun{
          D(signal{},Normal_Distribution{},mean<signal>{},stddev<signal>{}),
          F(mean<signal>{},
              [](auto t, auto baseline_,auto drift_,
                 auto A0_f0_, auto ph0_, auto tau0_,
                 auto A1_, auto f1_, auto ph1_, auto tau1_ ,
                 auto A2_, auto f2_, auto ph2_, auto tau2_ )
              {
                //                using a=typename decltype(t)::t;
                //                using b=typename decltype(baseline_)::baseline;
                //                using c=typename decltype(drift_)::drift;
                //              using d=typename decltype(A0_)::A0;
                //                using e=typename decltype(tau0_)::tau;

                return baseline_+(drift_*t)+
                       std::get<0>(A0_f0_.value())()*exp(
                                                           -t/tau0_
                                                           )*
                           cos(2*PI*std::get<1>(A0_f0_.value())()*GHz_f*t*ps_f
                               +ph0_)+
                       A1_*exp(-t/tau1_)*cos(2*PI*f1_*GHz_f*t*ps_f+ph1_)+
                       A2_*exp(-t/tau2_)*cos(2*PI*f2_*GHz_f*t*ps_f+ph2_);
              },
              delay{},baseline{},drift{},
              std::tuple<ind<Amplitude,0>,ind<Frecuency,0>>{},ind<phase,0>{},ind<tau,0>{},
              ind<Amplitude,1>{},ind<Frecuency,1>{},ind<phase,1>{},ind<tau,1>{},
              ind<Amplitude,2>{},ind<Frecuency,2>{},ind<phase,2>{},ind<tau,2>{})
      };


  auto data=vector_space{x_i(delay{}, vec<delay>{}),
                           x_i(signal{},vec<delay>{})};

  using  data_fields=Cs<delay,signal>;
  std::string fname_a="antena_data_1.txt";
  std::ifstream fi(fname_a.c_str());
  from_DataFrame(fi,data);


  auto data_3=get_data_index<3>();
  std::string fname_3="m03.txt";
  std::ifstream fi_3(fname_3.c_str());
  from_DataFrame(fi_3,data_3);

  auto data_4=get_data_index<4>();
  std::string fname_4="m04.txt";
  std::ifstream fi_4(fname_4.c_str());
  from_DataFrame(fi_4,data_4);

  auto data_5=get_data_index<5>();
  std::string fname_5="m05.txt";
  std::ifstream fi_5(fname_5.c_str());
  from_DataFrame(fi_5,data_5);


  auto data_8=get_data_index<8>();
  std::string fname_8="m08.txt";
  std::ifstream fi_8(fname_8.c_str());
  from_DataFrame(fi_8,data_8);

  auto data_9=get_data_index<9>();
  std::string fname_9="m09.txt";
  std::ifstream fi_9(fname_9.c_str());
  from_DataFrame(fi_9,data_9);


  auto data_10=get_data_index<10>();
  std::string fname_10="m10.txt";
  std::ifstream fi_10(fname_10.c_str());
  from_DataFrame(fi_10,data_10);



  auto totalmodel=mymodel+myprior_dist+myprior_transf+myprior_values;
  auto totalmodel_3=mymodel_3+myprior_dist_3+myprior_transf_3+myprior_values_3;

  auto totalmodel_40=totalmodel_index(std::index_sequence<0,1,2>{},std::index_sequence<0,1>{},std::index_sequence<>{})+my_common_prior_values();
  auto totalmodel_41=totalmodel_index(std::index_sequence<0,1,2>{},std::index_sequence<0,1,2>{},std::index_sequence<>{})+my_common_prior_values();

  auto totalmodel_41_1=totalmodel_index(std::index_sequence<1>{},std::index_sequence<0,1,2>{},std::index_sequence<>{})+my_common_prior_values();
  auto totalmodel_41_2=totalmodel_index(std::index_sequence<2>{},std::index_sequence<0,1,2>{},std::index_sequence<>{})+my_common_prior_values();


  auto totalmodel_43=totalmodel_index(std::index_sequence<0,1,2>{},std::index_sequence<0,1,2,3>{},std::index_sequence<>{})+my_common_prior_values();

  auto totalmodel_42=totalmodel_index(std::index_sequence<0,1,2>{},std::index_sequence<0,1,2>{},std::index_sequence<3>{})+my_common_prior_values();
  // using test=typename decltype (myprior_transf_index(std::index_sequence<0,1,2>{},std::index_sequence<0,1,2>{},std::index_sequence<3>{}))::ger;

  auto data_all=data_3+data_4+data_5+data_8+data_9+data_10;

  std::random_device rd;
  auto initseed = 0;
  //rd();

  std::mt19937 mt(initseed);

  auto mtv=v(std::move(mt));
  //auto [par_42,var_42]=sample(totalmodel_42,mtv,data_all);

  auto data_time=data|myselect<Cs<delay>>{};


  std::mt19937 mt2(initseed);
  auto mtv2=v(std::move(mt2));

  auto data2=data_time;


  auto [par, variables, predictions]=simulate(totalmodel,mtv,data);
  auto [par2, variables2]=sample(totalmodel,mtv2,data);

  auto logPriorv=logPrior(totalmodel,data,par2,variables2);
  auto logLikv=logLikelihood(totalmodel,data,par2,variables2);

  std::cerr<<"\nlogPriorv\n"<<logPriorv;
  std::cerr<<"\nlogPriorv\n"<<logLikv;
  auto dpar=Self_Derivative(par);
  auto dvariables3 =calculate(totalmodel,data,dpar);

  auto dlogPriorv=logPrior(totalmodel,data,dpar,dvariables3);
  auto dlogL=logLikelihood(totalmodel,data,dpar,dvariables3);

  auto fimPriorv=fimPrior(totalmodel,data,dpar,dvariables3);
  auto fimLikv=fimLikelihood(totalmodel,data,dpar,dvariables3);


  //  std::cerr<<"\ndlogPriorv\n"<<dlogPriorv;
  //  std::cerr<<"\ndlogLikv\n"<<dlogL;

  //  std::cerr<<"\ndvariables\n"<<dvariables3;


  //    //auto s=sample(totalmodel,std::move(data2),mt);


  // // auto data_sim=s| myselect<data_fields>{};


  //  std::cerr << "parameters \n"<<par <<std::endl;
  //  std::cerr << "parameters 2\n"<<par2 <<std::endl;
  //  std::cerr << "dparameters \n"<<dpar <<std::endl;

  //  std::cerr << "data \n"<<data <<std::endl;
  //  std::cerr << "predictions \n"<<predictions <<std::endl;


  // auto logL=logP(totalmodel,s);
  auto be=std::vector<double>{0,1e-4,3e-4,1e-3,3e-3,1e-2,2e-2,0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
  vector_field<vec<beta_ei>,v<double,dimension_less>> betas;
  auto pbe=betas.begin();
  for (auto &e:be)
  {
    insert_at(betas,pbe,std::move(e));
    ++pbe[beta_ei{}]();

  }

  std::string fname_emcee="emcee";
  //std::ofstream f_emcee(fname_emcee.c_str());

  auto decimate_factor=std::vector<std::size_t>{1ul,1ul,10ul,50ul,100};
  auto decimate_factor_1=std::vector<std::size_t>{1ul,1ul,1ul,1ul,1};

  //  auto mcmc=parallel_emcee(totalmodel_41,data_42,betas,v<std::size_t,dimension_less>(100),initseed,100000,decimate_factor,fname_emcee);
  //auto mcmc=parallel_emcee(totalmodel_42,data_42,betas,v<std::size_t,dimension_less>(100),initseed,100000,decimate_factor,fname_emcee);
  //  auto mcmc=parallel_emcee(totalmodel_43,data_42,betas,v<std::size_t,dimension_less>(100),initseed,100000,decimate_factor,fname_emcee);
  auto mcmc=parallel_emcee_series(totalmodel_41_<5>(),data_all,betas,v<std::size_t,dimension_less>(20),initseed,100000,decimate_factor,fname_emcee);
  auto mcmc_par=parallel_emcee_parallel(totalmodel_41_<5>(),data_all,betas,v<std::size_t,dimension_less>(20),initseed,100000,decimate_factor,fname_emcee);

  // f_emcee.close();
  std::string fname="out.txt";
  std::ofstream f(fname.c_str());
  auto s=dpar;
  // to_DataFrame(f,variables+data);
  // to_DataFrame(f,s);


  auto qui=totalmodel;;
  //auto dlogL=vector_space(logP(qui,data_sim,dpar));Z


  //  std::cerr<<"\n\nFIM\n "<<fimLikv<<"\n";
  auto dpar_new=decltype (dpar){};
  auto fim_new=decltype(fimLikv){};
  auto dlogLve=vector_space(std::move(dlogL));
  auto dlogL_new=decltype(dlogLve){};
  to_DataFrame(f,dlogLve);
  f.close();

  auto s_new=decltype (s){};
  std::ifstream fe;
  fe.open(fname.c_str());
  if (fe.is_open())
  {
    //from_DataFrame(fe,dlogL_new);
    //      from_DataFrame(fe,fim_new);
    //     from_DataFrame(fe,s_new);

  }
  //  std::cerr<<"\ndlogL\n"<<dlogL;
  //  std::cerr<<"\n dlogL_new\n"<<dlogL_new;




  //  assert(dlogLve==dlogL_new);



  return 0;
}
