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

struct drift{using T=double; using unit=decltype (V_u{}/ps_u{}); constexpr static auto className=my_static_string("drift");};

struct Amplitude{using T=double; using unit=V_u; constexpr static auto className=my_static_string("Amplitude");};

struct Frecuency{using T=double; using unit=GHz_u; constexpr static auto className=my_static_string("Frecuency");};


struct phase{using T=double; using unit=dimension_less; constexpr static auto className=my_static_string("phase");};

struct tau{using T=double; using unit=ps_u; constexpr static auto className=my_static_string("tau");};



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
      x_i(stddev<ind<Log10_t<tau>,1>>{},v(1.0))};

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

  auto data=vector_space{x_i(delay{}, vec<delay>{}),
                           x_i(signal{},vec<delay>{})};

  using  data_fields=Cs<delay,signal>;
  auto fname_a="antena_data_1.txt";
  std::ifstream fi(fname_a);
  from_DataTable(fi,data);

  auto totalmodel=mymodel+myprior_dist+myprior_transf+myprior_values;

  auto data_time=data|myselect<Cs<delay>>{};


  std::random_device rd;
  auto initseed = 0;
  //rd();

  std::mt19937 mt(initseed);

  auto mtv=v(std::move(mt));
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


  std::cerr<<"\ndlogPriorv\n"<<dlogPriorv;
  std::cerr<<"\ndlogLikv\n"<<dlogL;

  std::cerr<<"\ndvariables\n"<<dvariables3;


    //auto s=sample(totalmodel,std::move(data2),mt);


 // auto data_sim=s| myselect<data_fields>{};


  std::cerr << "parameters \n"<<par <<std::endl;
  std::cerr << "parameters 2\n"<<par2 <<std::endl;
  std::cerr << "dparameters \n"<<dpar <<std::endl;

  std::cerr << "data \n"<<data <<std::endl;
  std::cerr << "predictions \n"<<predictions <<std::endl;


 // auto logL=logP(totalmodel,s);
  auto be=std::vector<double>{0,1e-4,3e-4,1e-3,3e-3,1e-2,3e-2,0.1,0.2,0.3,0.5,0.8,1.0};
  vector_field<vec<beta_ei>,v<double,dimension_less>> betas;
  auto pbe=betas.begin();
  for (auto &e:be)
  {
    insert_at(betas,pbe,{v<double,dimension_less>(std::move(e))});
    ++pbe[beta_ei{}]();

  }


  auto mcmc=parallel_emcee(totalmodel,data,betas,v<std::size_t,dimension_less>(10),initseed,200,std::cerr);

  std::string fname="out.txt";
  std::ofstream f(fname.c_str());
  auto s=dpar;
  to_DataFrame(f,s);


      auto qui=totalmodel;;
  //auto dlogL=vector_space(logP(qui,data_sim,dpar));


  std::cerr<<"\n\nFIM\n "<<fimLikv<<"\n";
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
    from_DataFrame(fe,dlogL_new);
    //      from_DataFrame(fe,fim_new);
    //     from_DataFrame(fe,s_new);

  }
  std::cerr<<"\ndlogL\n"<<dlogL;
  std::cerr<<"\n dlogL_new\n"<<dlogL_new;




//  assert(dlogLve==dlogL_new);


  assert(test_print_read([](auto& os, auto & x)->auto&{return to_DataFrame(os,x);},
                         [](auto& is, auto& x)->auto&{return from_DataFrame(is,x);},
                         fimLikv,"output2"));

  assert(test_output_extraction_operators(fimLikv,"fimLikv.txt"));


  return 0;
}
