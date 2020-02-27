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
                auto y1=((std::get<I0>(A0_.value())*exp(-t/std::get<I0>(tau0_.value()))*
                            cos(2*PI*std::get<I0>(f0_.value())*GHz_f*t*ps_f+std::get<I0>(ph0_.value())))+...);

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
                auto y1=((std::get<I0>(A0_.value())*exp(-t/std::get<I0>(tau0_.value()))*
                            cos(2*PI*std::get<I0>(f0_.value())*GHz_f*t*ps_f+std::get<I0>(ph0_.value())))+...);
                auto y2=((std::get<0>(A1_.value())*exp(-t/std::get<0>(tau1_.value()))*
                            cos(2*PI*std::get<0>(f1_.value())*GHz_f*t*ps_f+std::get<0>(ph1_.value()))));
                if constexpr (sizeof... (I1)==0)

                  return y0+y1+y2;
                else
                {
                  auto y3=((std::get<I1-I00>(A1_.value())*exp(-t/std::get<I1-I00>(tau1_.value()))*
                              cos(2*PI*std::get<I1-I00>(f1_.value())*GHz_f*t*ps_f+std::get<I1-I00>(ph1_.value())))+...);
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


int main(int argc, char **argv)
{
  std::string arg;
  std::cerr<<argv[0]<<"\n";

  if (argc>1)
  {
  arg=argv[1];
  std::cerr<<argv[1]<<"\n";
  }
  else
  {
    auto tpid=std::getenv("SLURM_TASK_PID");
    auto jid=std::getenv("SLURM_JOB_ID");

    char const* args[]={"pp_model_345_01_","pp_model_345_012_","pp_model_345_0123_","pp_model_345_012_4","pp_model_345_0_4","pp_model_345_0_45",
                 "pp_model_8910_01_","pp_model_8910_012_","pp_model_8910_0123_","pp_model_8910_012_4","pp_model_8910_0_4","pp_model_8910_0_45"
                 };
    if (tpid!=nullptr)
    {
      std::cerr<<"SLURM_JOB_ID="<<jid<<"\t"<<"SLURM_TASK_PID="<<tpid<<"\n";
      auto tid=std::atol(tpid);
      auto id=tid%12;
      arg=args[id];
      std::cerr<<"\t arg="<<arg<<"\n";
    }
  }


  auto data_3=get_data_index<3>();
//  std::string fname_3="m03.txt";
  std::string fname_3="m03_5.txt";
  std::ifstream fi_3(fname_3.c_str());
  from_DataFrame(fi_3,data_3);

  auto data_4=get_data_index<4>();
//  std::string fname_4="m04.txt";
  std::string fname_4="m04_5.txt";
  std::ifstream fi_4(fname_4.c_str());
  from_DataFrame(fi_4,data_4);

  auto data_5=get_data_index<5>();
//  std::string fname_5="m05.txt";
  std::string fname_5="m05_5.txt";
  std::ifstream fi_5(fname_5.c_str());
  from_DataFrame(fi_5,data_5);


  auto data_8=get_data_index<8>();
//  std::string fname_8="m08.txt";
  std::string fname_8="m08_5.txt";
  std::ifstream fi_8(fname_8.c_str());
  from_DataFrame(fi_8,data_8);

  auto data_9=get_data_index<9>();
//  std::string fname_9="m09.txt";
  std::string fname_9="m09_5.txt";
  std::ifstream fi_9(fname_9.c_str());
  from_DataFrame(fi_9,data_9);


  auto data_10=get_data_index<10>();
//  std::string fname_10="m10.txt";
  std::string fname_10="m10_5.txt";
  std::ifstream fi_10(fname_10.c_str());
  from_DataFrame(fi_10,data_10);




//  auto model_345_01_=totalmodel_index(std::index_sequence<3,4,5>{},std::index_sequence<0,1>{},std::index_sequence<>{})+my_common_prior_values();
//  auto model_345_0_4=totalmodel_index(std::index_sequence<3,4,5>{},std::index_sequence<0>{},std::index_sequence<4>{})+my_common_prior_values();
//  auto model_345_0_45=totalmodel_index(std::index_sequence<3,4,5>{},std::index_sequence<0>{},std::index_sequence<4,5>{})+my_common_prior_values();
//  auto model_345_012_=totalmodel_index(std::index_sequence<3,4,5>{},std::index_sequence<0,1,2>{},std::index_sequence<>{})+my_common_prior_values();
 // auto model_345_012_4=totalmodel_index(std::index_sequence<3,4,5>{},std::index_sequence<0,1,2>{},std::index_sequence<4>{})+my_common_prior_values();
//  auto model_345_0123_=totalmodel_index(std::index_sequence<3,4,5>{},std::index_sequence<0,1,2,3>{},std::index_sequence<>{})+my_common_prior_values();


//  auto model_8910_01_=totalmodel_index(std::index_sequence<8,9,10>{},std::index_sequence<0,1>{},std::index_sequence<>{})+my_common_prior_values();
//  auto model_8910_0_4=totalmodel_index(std::index_sequence<8,9,10>{},std::index_sequence<0>{},std::index_sequence<4>{})+my_common_prior_values();
//  auto model_8910_0_45=totalmodel_index(std::index_sequence<8,9,10>{},std::index_sequence<0>{},std::index_sequence<4,5>{})+my_common_prior_values();
  auto model_8910_012_=totalmodel_index(std::index_sequence<8,9,10>{},std::index_sequence<0,1,2>{},std::index_sequence<>{})+my_common_prior_values();
//  auto model_8910_012_4=totalmodel_index(std::index_sequence<8,9,10>{},std::index_sequence<0,1,2>{},std::index_sequence<4>{})+my_common_prior_values();
//  auto model_8910_0123_=totalmodel_index(std::index_sequence<8,9,10>{},std::index_sequence<0,1,2,3>{},std::index_sequence<>{})+my_common_prior_values();


  auto data_345=data_3+data_4+data_5;
  auto data_8910=data_8+data_9+data_10;

  std::random_device rd;
  auto initseed = rd();
  initseed=0;

  std::mt19937 mt(initseed);

  auto mtv=v(std::move(mt));



  auto be=std::vector<double>{0,1e-4,3e-4,1e-3,3e-3,1e-2,2e-2,0.05,0.1,0.15,0.2,0.3,0.4,0.6,0.8,1.0};
  vector_field<vec<beta_ei>,v<double,dimension_less>> betas;
  auto pbe=betas.begin();
  for (auto &e:be)
  {
    insert_at(betas,pbe,std::move(e));
    ++pbe[beta_ei{}]();
  }

  auto decimate_factor=std::vector<std::size_t>{1000ul,1000ul,10000ul,400000ul};
//  auto decimate_factor=std::vector<std::size_t>{1ul,1ul,1ul,4ul};

//  std::size_t maxiters=20;
  std::size_t maxiters=4000000;
  std::size_t nwalkers=16;

/*  if (arg=="pp_model_345_012_4")
  {
    maxiters*=4;
    parallel_emcee_parallel_parallel_for(model_345_012_4,data_345,betas,v<std::size_t,dimension_less>(nwalkers),initseed,maxiters,decimate_factor,arg);
  }
  else if (arg=="p_model_345_012_4")
  {
    parallel_emcee_parallel(model_345_012_4,data_all,betas,v<std::size_t,dimension_less>(nwalkers),initseed,maxiters,decimate_factor,arg);
  }*/
/*  else if (arg=="sp_model_345_012_4")
  {
    parallel_emcee_series_parallel_for(model_345_012_4,data_345,betas,v<std::size_t,dimension_less>(nwalkers),initseed,maxiters,decimate_factor,arg);
  }
  else if (arg=="ss_model_345_012_4")
  {
    parallel_emcee_series(model_345_012_4,data_all,betas,v<std::size_t,dimension_less>(nwalkers),initseed,maxiters,decimate_factor,arg);
  }*/
  /*
  else if (arg=="pp_model_345_0_4")
  {
    parallel_emcee_parallel_parallel_for(model_345_0_4,data_345,betas,v<std::size_t,dimension_less>(nwalkers),initseed,maxiters,decimate_factor,arg);
  }
  else if (arg=="pp_model_345_01_")
  {
    parallel_emcee_parallel_parallel_for(model_345_01_,data_345,betas,v<std::size_t,dimension_less>(nwalkers),initseed,maxiters,decimate_factor,arg);
  }
  else if (arg=="pp_model_345_012_")
  {
    maxiters*=2;
    parallel_emcee_parallel_parallel_for(model_345_012_,data_345,betas,v<std::size_t,dimension_less>(nwalkers),initseed,maxiters,decimate_factor,arg);
  }
  else if (arg=="pp_model_345_0123_")
  {
    maxiters*=4;
    parallel_emcee_parallel_parallel_for(model_345_0123_,data_345,betas,v<std::size_t,dimension_less>(nwalkers),initseed,maxiters,decimate_factor,arg);
  }
  else if (arg=="pp_model_345_0_45")
  {
    maxiters*=4;
    parallel_emcee_parallel_parallel_for(model_345_0_45,data_345,betas,v<std::size_t,dimension_less>(nwalkers),initseed,maxiters,decimate_factor,arg);
  }
  else if (arg=="pp_model_8910_0_4")
  {
    parallel_emcee_parallel_parallel_for(model_8910_0_4,data_8910,betas,v<std::size_t,dimension_less>(nwalkers),initseed,maxiters,decimate_factor,arg);
  }
  else if (arg=="pp_model_8910_012_4")
  {
    maxiters*=4;
    parallel_emcee_parallel_parallel_for(model_8910_012_4,data_8910,betas,v<std::size_t,dimension_less>(nwalkers),initseed,maxiters,decimate_factor,arg);
  }
  else if (arg=="pp_model_8910_01_")
  {
    parallel_emcee_parallel_parallel_for(model_8910_01_,data_8910,betas,v<std::size_t,dimension_less>(nwalkers),initseed,maxiters,decimate_factor,arg);
  }

  else */if (arg=="p_model_8910_012_")
  {
    maxiters*=2;
    parallel_emcee_parallel(model_8910_012_,data_8910,betas,v<std::size_t,dimension_less>(nwalkers),initseed,maxiters,decimate_factor,arg);
  }
  else if (arg=="pp_model_8910_012_")
  {
    maxiters*=2;
    parallel_emcee_parallel_parallel_for(model_8910_012_,data_8910,betas,v<std::size_t,dimension_less>(nwalkers),initseed,maxiters,decimate_factor,arg);
  }
  else if (arg=="pp_model_8910_012_q")
  {
    maxiters*=2;
    parallel_emcee_parallel_parallel_for_q(model_8910_012_,data_8910,betas,v<std::size_t,dimension_less>(nwalkers),initseed,maxiters,decimate_factor,arg);
  }
  /*  else if (arg=="pp_model_8910_0123_")
  {
    maxiters*=4;
    parallel_emcee_parallel_parallel_for(model_8910_0123_,data_8910,betas,v<std::size_t,dimension_less>(nwalkers),initseed,maxiters,decimate_factor,arg);
  }
  else if (arg=="pp_model_8910_0123_q")
  {
    maxiters*=4;
    parallel_emcee_parallel_parallel_for_q(model_8910_0123_,data_8910,betas,v<std::size_t,dimension_less>(nwalkers),initseed,maxiters,decimate_factor,arg);
  }
 /* else if (arg=="pp_model_8910_0_45")
  {
    maxiters*=4;
    parallel_emcee_parallel_parallel_for(model_8910_0_45,data_8910,betas,v<std::size_t,dimension_less>(nwalkers),initseed,maxiters,decimate_factor,arg);
  }
*/
    return 0;


}
