/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__GluSynapse_andras
#define _nrn_initial _nrn_initial__GluSynapse_andras
#define nrn_cur _nrn_cur__GluSynapse_andras
#define _nrn_current _nrn_current__GluSynapse_andras
#define nrn_jacob _nrn_jacob__GluSynapse_andras
#define nrn_state _nrn_state__GluSynapse_andras
#define _net_receive _net_receive__GluSynapse_andras 
#define setRNG setRNG__GluSynapse_andras 
#define state state__GluSynapse_andras 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define tau_d_AMPA _p[0]
#define gmax0_AMPA _p[1]
#define gmax_d_AMPA _p[2]
#define gmax_p_AMPA _p[3]
#define gmax_NMDA _p[4]
#define Use0_TM _p[5]
#define Dep_TM _p[6]
#define Fac_TM _p[7]
#define Nrrp_TM _p[8]
#define Use_d_TM _p[9]
#define Use_p_TM _p[10]
#define volume_CR _p[11]
#define theta_d_GB _p[12]
#define theta_p_GB _p[13]
#define rho0_GB _p[14]
#define synapseID _p[15]
#define verbose _p[16]
#define selected_for_report _p[17]
#define g_AMPA _p[18]
#define g_NMDA _p[19]
#define ica_NMDA _p[20]
#define ica_VDCC _p[21]
#define dep_GB _p[22]
#define pot_GB _p[23]
#define vsyn _p[24]
#define i _p[25]
#define A_AMPA _p[26]
#define B_AMPA _p[27]
#define gmax_AMPA _p[28]
#define A_NMDA _p[29]
#define B_NMDA _p[30]
#define Use_TM _p[31]
#define m_VDCC _p[32]
#define h_VDCC _p[33]
#define cai_CR _p[34]
#define rho_GB _p[35]
#define effcai_GB _p[36]
#define usingR123 _p[37]
#define DA_AMPA _p[38]
#define DB_AMPA _p[39]
#define Dgmax_AMPA _p[40]
#define DA_NMDA _p[41]
#define DB_NMDA _p[42]
#define DUse_TM _p[43]
#define Dm_VDCC _p[44]
#define Dh_VDCC _p[45]
#define Dcai_CR _p[46]
#define Drho_GB _p[47]
#define Deffcai_GB _p[48]
#define v _p[49]
#define _g _p[50]
#define _tsav _p[51]
#define _nd_area  *_ppvar[0]._pval
#define rng_TM	*_ppvar[2]._pval
#define _p_rng_TM	_ppvar[2]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  2;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static double _hoc_bbsavestate();
 static double _hoc_brand();
 static double _hoc_nernst();
 static double _hoc_setRNG();
 static double _hoc_urand();
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern Prop* nrn_point_prop_;
 static int _pointtype;
 static void* _hoc_create_pnt(_ho) Object* _ho; { void* create_point_process();
 return create_point_process(_pointtype, _ho);
}
 static void _hoc_destroy_pnt();
 static double _hoc_loc_pnt(_vptr) void* _vptr; {double loc_point_process();
 return loc_point_process(_pointtype, _vptr);
}
 static double _hoc_has_loc(_vptr) void* _vptr; {double has_loc_point();
 return has_loc_point(_vptr);
}
 static double _hoc_get_loc_pnt(_vptr)void* _vptr; {
 double get_loc_point_process(); return (get_loc_point_process(_vptr));
}
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata(void* _vptr) { Prop* _prop;
 _prop = ((Point_process*)_vptr)->_prop;
   _setdata(_prop);
 }
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 0,0
};
 static Member_func _member_func[] = {
 "loc", _hoc_loc_pnt,
 "has_loc", _hoc_has_loc,
 "get_loc", _hoc_get_loc_pnt,
 "bbsavestate", _hoc_bbsavestate,
 "brand", _hoc_brand,
 "nernst", _hoc_nernst,
 "setRNG", _hoc_setRNG,
 "urand", _hoc_urand,
 0, 0
};
#define bbsavestate bbsavestate_GluSynapse_andras
#define brand brand_GluSynapse_andras
#define nernst nernst_GluSynapse_andras
#define urand urand_GluSynapse_andras
 extern double bbsavestate( _threadargsproto_ );
 extern double brand( _threadargsprotocomma_ double , double );
 extern double nernst( _threadargsprotocomma_ double , double , double );
 extern double urand( _threadargsproto_ );
 /* declare global and static user variables */
#define E_NMDA E_NMDA_GluSynapse_andras
 double E_NMDA = -3;
#define E_AMPA E_AMPA_GluSynapse_andras
 double E_AMPA = 0;
#define cao_CR cao_CR_GluSynapse_andras
 double cao_CR = 2;
#define gamma_p_GB gamma_p_GB_GluSynapse_andras
 double gamma_p_GB = 450;
#define gamma_d_GB gamma_d_GB_GluSynapse_andras
 double gamma_d_GB = 100;
#define gamma_ca_CR gamma_ca_CR_GluSynapse_andras
 double gamma_ca_CR = 0.04;
#define gca_bar_VDCC gca_bar_VDCC_GluSynapse_andras
 double gca_bar_VDCC = 0.0744;
#define htau_VDCC htau_VDCC_GluSynapse_andras
 double htau_VDCC = 27;
#define kh_VDCC kh_VDCC_GluSynapse_andras
 double kh_VDCC = -9.2;
#define km_VDCC km_VDCC_GluSynapse_andras
 double km_VDCC = 9.5;
#define ljp_VDCC ljp_VDCC_GluSynapse_andras
 double ljp_VDCC = 0;
#define min_ca_CR min_ca_CR_GluSynapse_andras
 double min_ca_CR = 7e-005;
#define mtau_VDCC mtau_VDCC_GluSynapse_andras
 double mtau_VDCC = 1;
#define mgo_NMDA mgo_NMDA_GluSynapse_andras
 double mgo_NMDA = 1;
#define rho_star_GB rho_star_GB_GluSynapse_andras
 double rho_star_GB = 0.5;
#define slope_NMDA slope_NMDA_GluSynapse_andras
 double slope_NMDA = 0.072;
#define scale_NMDA scale_NMDA_GluSynapse_andras
 double scale_NMDA = 2.552;
#define tau_effca_GB tau_effca_GB_GluSynapse_andras
 double tau_effca_GB = 200;
#define tau_exp_GB tau_exp_GB_GluSynapse_andras
 double tau_exp_GB = 100;
#define tau_ind_GB tau_ind_GB_GluSynapse_andras
 double tau_ind_GB = 70;
#define tau_ca_CR tau_ca_CR_GluSynapse_andras
 double tau_ca_CR = 12;
#define tau_d_NMDA tau_d_NMDA_GluSynapse_andras
 double tau_d_NMDA = 70;
#define tau_r_NMDA tau_r_NMDA_GluSynapse_andras
 double tau_r_NMDA = 0.29;
#define tau_r_AMPA tau_r_AMPA_GluSynapse_andras
 double tau_r_AMPA = 0.2;
#define vhh_VDCC vhh_VDCC_GluSynapse_andras
 double vhh_VDCC = -39;
#define vhm_VDCC vhm_VDCC_GluSynapse_andras
 double vhm_VDCC = -5.9;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "tau_r_AMPA_GluSynapse_andras", "ms",
 "E_AMPA_GluSynapse_andras", "mV",
 "mgo_NMDA_GluSynapse_andras", "mM",
 "scale_NMDA_GluSynapse_andras", "mM",
 "slope_NMDA_GluSynapse_andras", "/mV",
 "tau_r_NMDA_GluSynapse_andras", "ms",
 "tau_d_NMDA_GluSynapse_andras", "ms",
 "E_NMDA_GluSynapse_andras", "mV",
 "gca_bar_VDCC_GluSynapse_andras", "nS/um2",
 "ljp_VDCC_GluSynapse_andras", "mV",
 "vhm_VDCC_GluSynapse_andras", "mV",
 "km_VDCC_GluSynapse_andras", "mV",
 "vhh_VDCC_GluSynapse_andras", "mV",
 "kh_VDCC_GluSynapse_andras", "mV",
 "mtau_VDCC_GluSynapse_andras", "ms",
 "htau_VDCC_GluSynapse_andras", "ms",
 "gamma_ca_CR_GluSynapse_andras", "1",
 "tau_ca_CR_GluSynapse_andras", "ms",
 "min_ca_CR_GluSynapse_andras", "mM",
 "cao_CR_GluSynapse_andras", "mM",
 "rho_star_GB_GluSynapse_andras", "1",
 "tau_ind_GB_GluSynapse_andras", "s",
 "tau_exp_GB_GluSynapse_andras", "s",
 "tau_effca_GB_GluSynapse_andras", "ms",
 "gamma_d_GB_GluSynapse_andras", "1",
 "gamma_p_GB_GluSynapse_andras", "1",
 "tau_d_AMPA", "ms",
 "gmax0_AMPA", "nS",
 "gmax_d_AMPA", "nS",
 "gmax_p_AMPA", "nS",
 "gmax_NMDA", "nS",
 "Use0_TM", "1",
 "Dep_TM", "ms",
 "Fac_TM", "ms",
 "Nrrp_TM", "1",
 "Use_d_TM", "1",
 "Use_p_TM", "1",
 "volume_CR", "um3",
 "theta_d_GB", "us/liter",
 "theta_p_GB", "us/liter",
 "rho0_GB", "1",
 "A_AMPA", "1",
 "B_AMPA", "1",
 "gmax_AMPA", "nS",
 "A_NMDA", "1",
 "B_NMDA", "1",
 "Use_TM", "1",
 "m_VDCC", "1",
 "h_VDCC", "1",
 "cai_CR", "mM",
 "rho_GB", "1",
 "effcai_GB", "us/liter",
 "g_AMPA", "uS",
 "g_NMDA", "uS",
 "ica_NMDA", "nA",
 "ica_VDCC", "nA",
 "dep_GB", "1",
 "pot_GB", "1",
 "vsyn", "mV",
 "i", "nA",
 0,0
};
 static double A_NMDA0 = 0;
 static double A_AMPA0 = 0;
 static double B_NMDA0 = 0;
 static double B_AMPA0 = 0;
 static double Use_TM0 = 0;
 static double cai_CR0 = 0;
 static double delta_t = 0.01;
 static double effcai_GB0 = 0;
 static double gmax_AMPA0 = 0;
 static double h_VDCC0 = 0;
 static double m_VDCC0 = 0;
 static double rho_GB0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "tau_r_AMPA_GluSynapse_andras", &tau_r_AMPA_GluSynapse_andras,
 "E_AMPA_GluSynapse_andras", &E_AMPA_GluSynapse_andras,
 "mgo_NMDA_GluSynapse_andras", &mgo_NMDA_GluSynapse_andras,
 "scale_NMDA_GluSynapse_andras", &scale_NMDA_GluSynapse_andras,
 "slope_NMDA_GluSynapse_andras", &slope_NMDA_GluSynapse_andras,
 "tau_r_NMDA_GluSynapse_andras", &tau_r_NMDA_GluSynapse_andras,
 "tau_d_NMDA_GluSynapse_andras", &tau_d_NMDA_GluSynapse_andras,
 "E_NMDA_GluSynapse_andras", &E_NMDA_GluSynapse_andras,
 "gca_bar_VDCC_GluSynapse_andras", &gca_bar_VDCC_GluSynapse_andras,
 "ljp_VDCC_GluSynapse_andras", &ljp_VDCC_GluSynapse_andras,
 "vhm_VDCC_GluSynapse_andras", &vhm_VDCC_GluSynapse_andras,
 "km_VDCC_GluSynapse_andras", &km_VDCC_GluSynapse_andras,
 "vhh_VDCC_GluSynapse_andras", &vhh_VDCC_GluSynapse_andras,
 "kh_VDCC_GluSynapse_andras", &kh_VDCC_GluSynapse_andras,
 "mtau_VDCC_GluSynapse_andras", &mtau_VDCC_GluSynapse_andras,
 "htau_VDCC_GluSynapse_andras", &htau_VDCC_GluSynapse_andras,
 "gamma_ca_CR_GluSynapse_andras", &gamma_ca_CR_GluSynapse_andras,
 "tau_ca_CR_GluSynapse_andras", &tau_ca_CR_GluSynapse_andras,
 "min_ca_CR_GluSynapse_andras", &min_ca_CR_GluSynapse_andras,
 "cao_CR_GluSynapse_andras", &cao_CR_GluSynapse_andras,
 "rho_star_GB_GluSynapse_andras", &rho_star_GB_GluSynapse_andras,
 "tau_ind_GB_GluSynapse_andras", &tau_ind_GB_GluSynapse_andras,
 "tau_exp_GB_GluSynapse_andras", &tau_exp_GB_GluSynapse_andras,
 "tau_effca_GB_GluSynapse_andras", &tau_effca_GB_GluSynapse_andras,
 "gamma_d_GB_GluSynapse_andras", &gamma_d_GB_GluSynapse_andras,
 "gamma_p_GB_GluSynapse_andras", &gamma_p_GB_GluSynapse_andras,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 
#define _watch_array _ppvar + 4 
 static void _hoc_destroy_pnt(_vptr) void* _vptr; {
   Prop* _prop = ((Point_process*)_vptr)->_prop;
   if (_prop) { _nrn_free_watch(_prop->dparam, 4, 5);}
   destroy_point_process(_vptr);
}
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[9]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"GluSynapse_andras",
 "tau_d_AMPA",
 "gmax0_AMPA",
 "gmax_d_AMPA",
 "gmax_p_AMPA",
 "gmax_NMDA",
 "Use0_TM",
 "Dep_TM",
 "Fac_TM",
 "Nrrp_TM",
 "Use_d_TM",
 "Use_p_TM",
 "volume_CR",
 "theta_d_GB",
 "theta_p_GB",
 "rho0_GB",
 "synapseID",
 "verbose",
 "selected_for_report",
 0,
 "g_AMPA",
 "g_NMDA",
 "ica_NMDA",
 "ica_VDCC",
 "dep_GB",
 "pot_GB",
 "vsyn",
 "i",
 0,
 "A_AMPA",
 "B_AMPA",
 "gmax_AMPA",
 "A_NMDA",
 "B_NMDA",
 "Use_TM",
 "m_VDCC",
 "h_VDCC",
 "cai_CR",
 "rho_GB",
 "effcai_GB",
 0,
 "rng_TM",
 0};
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
  if (nrn_point_prop_) {
	_prop->_alloc_seq = nrn_point_prop_->_alloc_seq;
	_p = nrn_point_prop_->param;
	_ppvar = nrn_point_prop_->dparam;
 }else{
 	_p = nrn_prop_data_alloc(_mechtype, 52, _prop);
 	/*initialize range parameters*/
 	tau_d_AMPA = 1.7;
 	gmax0_AMPA = 1;
 	gmax_d_AMPA = 1;
 	gmax_p_AMPA = 2;
 	gmax_NMDA = 0.55;
 	Use0_TM = 0.5;
 	Dep_TM = 100;
 	Fac_TM = 10;
 	Nrrp_TM = 1;
 	Use_d_TM = 0.2;
 	Use_p_TM = 0.8;
 	volume_CR = 0.087;
 	theta_d_GB = 0.006;
 	theta_p_GB = 0.001;
 	rho0_GB = 0;
 	synapseID = 0;
 	verbose = 0;
 	selected_for_report = 0;
  }
 	_prop->param = _p;
 	_prop->param_size = 52;
  if (!nrn_point_prop_) {
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 10, _prop);
  }
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 "cai_CR", 1e-006,
 "effcai_GB", 0.001,
 0,0
};
 
#define _tqitem &(_ppvar[3]._pvoid)
 static void _net_receive(Point_process*, double*, double);
 static void _net_init(Point_process*, double*, double);
 static void _thread_mem_init(Datum*);
 static void _thread_cleanup(Datum*);
 static void bbcore_write(double*, int*, int*, int*, _threadargsproto_);
 extern void hoc_reg_bbcore_write(int, void(*)(double*, int*, int*, int*, _threadargsproto_));
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _GluSynapse_andras_reg() {
	int _vectorized = 1;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex, 5,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
  _extcall_thread = (Datum*)ecalloc(4, sizeof(Datum));
  _thread_mem_init(_extcall_thread);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 1, _thread_mem_init);
     _nrn_thread_reg(_mechtype, 0, _thread_cleanup);
   hoc_reg_bbcore_write(_mechtype, bbcore_write);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 52, 10);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "bbcorepointer");
  hoc_register_dparam_semantics(_mechtype, 3, "netsend");
  hoc_register_dparam_semantics(_mechtype, 4, "watch");
  hoc_register_dparam_semantics(_mechtype, 5, "watch");
  hoc_register_dparam_semantics(_mechtype, 6, "watch");
  hoc_register_dparam_semantics(_mechtype, 7, "watch");
  hoc_register_dparam_semantics(_mechtype, 8, "watch");
  hoc_register_dparam_semantics(_mechtype, 9, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_init[_mechtype] = _net_init;
 pnt_receive_size[_mechtype] = 5;
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 GluSynapse_andras C:/Users/mkbli/Dropbox/deep_dendrite/mods/GluSynapse_andras.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double FARADAY = 96485.3;
 static double PI = 3.14159;
 static double R = 8.3145;
static int _reset;
static char *modelname = "Glutamatergic synapse";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int setRNG(_threadargsproto_);
 
#define _deriv1_advance _thread[0]._i
#define _dith1 1
#define _recurse _thread[2]._i
#define _newtonspace1 _thread[3]._pvoid
extern void* nrn_cons_newtonspace(int);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist2[11];
  static int _slist1[11], _dlist1[11];
 static int state(_threadargsproto_);
 
/*VERBATIM*/
/**
 * This Verbatim block is needed to generate random numbers from a uniform
 * distribution U(0, 1).
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "nrnran123.h"
double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   double _lminf_VDCC , _lhinf_VDCC ;
 DA_AMPA = - A_AMPA / tau_r_AMPA ;
   DB_AMPA = - B_AMPA / tau_d_AMPA ;
   Dgmax_AMPA = ( gmax_d_AMPA + rho_GB * ( gmax_p_AMPA - gmax_d_AMPA ) - gmax_AMPA ) / ( ( 1e3 ) * tau_exp_GB ) ;
   DA_NMDA = - A_NMDA / tau_r_NMDA ;
   DB_NMDA = - B_NMDA / tau_d_NMDA ;
   DUse_TM = ( Use_d_TM + rho_GB * ( Use_p_TM - Use_d_TM ) - Use_TM ) / ( ( 1e3 ) * tau_exp_GB ) ;
   _lminf_VDCC = 1.0 / ( 1.0 + exp ( ( ( vhm_VDCC - ljp_VDCC ) - v ) / km_VDCC ) ) ;
   _lhinf_VDCC = 1.0 / ( 1.0 + exp ( ( ( vhh_VDCC - ljp_VDCC ) - v ) / kh_VDCC ) ) ;
   Dm_VDCC = ( _lminf_VDCC - m_VDCC ) / mtau_VDCC ;
   Dh_VDCC = ( _lhinf_VDCC - h_VDCC ) / htau_VDCC ;
   Dcai_CR = - ( 1e-9 ) * ( ica_NMDA + ica_VDCC ) * gamma_ca_CR / ( ( 1e-15 ) * volume_CR * 2.0 * FARADAY ) - ( cai_CR - min_ca_CR ) / tau_ca_CR ;
   Deffcai_GB = - effcai_GB / tau_effca_GB + ( cai_CR - min_ca_CR ) ;
   Drho_GB = ( - rho_GB * ( 1.0 - rho_GB ) * ( rho_star_GB - rho_GB ) + pot_GB * gamma_p_GB * ( 1.0 - rho_GB ) - dep_GB * gamma_d_GB * rho_GB ) / ( ( 1e3 ) * tau_ind_GB ) ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 double _lminf_VDCC , _lhinf_VDCC ;
 DA_AMPA = DA_AMPA  / (1. - dt*( ( - 1.0 ) / tau_r_AMPA )) ;
 DB_AMPA = DB_AMPA  / (1. - dt*( ( - 1.0 ) / tau_d_AMPA )) ;
 Dgmax_AMPA = Dgmax_AMPA  / (1. - dt*( ( ( ( - 1.0 ) ) ) / ( ( 1e3 ) * tau_exp_GB ) )) ;
 DA_NMDA = DA_NMDA  / (1. - dt*( ( - 1.0 ) / tau_r_NMDA )) ;
 DB_NMDA = DB_NMDA  / (1. - dt*( ( - 1.0 ) / tau_d_NMDA )) ;
 DUse_TM = DUse_TM  / (1. - dt*( ( ( ( - 1.0 ) ) ) / ( ( 1e3 ) * tau_exp_GB ) )) ;
 _lminf_VDCC = 1.0 / ( 1.0 + exp ( ( ( vhm_VDCC - ljp_VDCC ) - v ) / km_VDCC ) ) ;
 _lhinf_VDCC = 1.0 / ( 1.0 + exp ( ( ( vhh_VDCC - ljp_VDCC ) - v ) / kh_VDCC ) ) ;
 Dm_VDCC = Dm_VDCC  / (1. - dt*( ( ( ( - 1.0 ) ) ) / mtau_VDCC )) ;
 Dh_VDCC = Dh_VDCC  / (1. - dt*( ( ( ( - 1.0 ) ) ) / htau_VDCC )) ;
 Dcai_CR = Dcai_CR  / (1. - dt*( ( - ( ( 1.0 ) ) / tau_ca_CR ) )) ;
 Deffcai_GB = Deffcai_GB  / (1. - dt*( ( - 1.0 ) / tau_effca_GB )) ;
 Drho_GB = Drho_GB  / (1. - dt*( ( ( (( (( - 1.0 )*( ( 1.0 - rho_GB ) ) + ( - rho_GB )*( ( ( - 1.0 ) ) )) )*( ( rho_star_GB - rho_GB ) ) + ( - rho_GB * ( 1.0 - rho_GB ) )*( ( ( - 1.0 ) ) )) + ( pot_GB * gamma_p_GB )*( ( ( - 1.0 ) ) ) - ( dep_GB * gamma_d_GB )*( 1.0 ) ) ) / ( ( 1e3 ) * tau_ind_GB ) )) ;
  return 0;
}
 /*END CVODE*/
 
static int state (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset=0; int error = 0;
 { double* _savstate1 = _thread[_dith1]._pval;
 double* _dlist2 = _thread[_dith1]._pval + 11;
 int _counte = -1;
 if (!_recurse) {
 _recurse = 1;
 {int _id; for(_id=0; _id < 11; _id++) { _savstate1[_id] = _p[_slist1[_id]];}}
 error = nrn_newton_thread(_newtonspace1, 11,_slist2, _p, state, _dlist2, _ppvar, _thread, _nt);
 _recurse = 0; if(error) {abort_run(error);}}
 {
   double _lminf_VDCC , _lhinf_VDCC ;
 DA_AMPA = - A_AMPA / tau_r_AMPA ;
   DB_AMPA = - B_AMPA / tau_d_AMPA ;
   Dgmax_AMPA = ( gmax_d_AMPA + rho_GB * ( gmax_p_AMPA - gmax_d_AMPA ) - gmax_AMPA ) / ( ( 1e3 ) * tau_exp_GB ) ;
   DA_NMDA = - A_NMDA / tau_r_NMDA ;
   DB_NMDA = - B_NMDA / tau_d_NMDA ;
   DUse_TM = ( Use_d_TM + rho_GB * ( Use_p_TM - Use_d_TM ) - Use_TM ) / ( ( 1e3 ) * tau_exp_GB ) ;
   _lminf_VDCC = 1.0 / ( 1.0 + exp ( ( ( vhm_VDCC - ljp_VDCC ) - v ) / km_VDCC ) ) ;
   _lhinf_VDCC = 1.0 / ( 1.0 + exp ( ( ( vhh_VDCC - ljp_VDCC ) - v ) / kh_VDCC ) ) ;
   Dm_VDCC = ( _lminf_VDCC - m_VDCC ) / mtau_VDCC ;
   Dh_VDCC = ( _lhinf_VDCC - h_VDCC ) / htau_VDCC ;
   Dcai_CR = - ( 1e-9 ) * ( ica_NMDA + ica_VDCC ) * gamma_ca_CR / ( ( 1e-15 ) * volume_CR * 2.0 * FARADAY ) - ( cai_CR - min_ca_CR ) / tau_ca_CR ;
   Deffcai_GB = - effcai_GB / tau_effca_GB + ( cai_CR - min_ca_CR ) ;
   Drho_GB = ( - rho_GB * ( 1.0 - rho_GB ) * ( rho_star_GB - rho_GB ) + pot_GB * gamma_p_GB * ( 1.0 - rho_GB ) - dep_GB * gamma_d_GB * rho_GB ) / ( ( 1e3 ) * tau_ind_GB ) ;
   {int _id; for(_id=0; _id < 11; _id++) {
if (_deriv1_advance) {
 _dlist2[++_counte] = _p[_dlist1[_id]] - (_p[_slist1[_id]] - _savstate1[_id])/dt;
 }else{
_dlist2[++_counte] = _p[_slist1[_id]] - _savstate1[_id];}}}
 } }
 return _reset;}
 
static double _watch1_cond(_pnt) Point_process* _pnt; {
 	double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
	_thread= (Datum*)0; _nt = (_NrnThread*)_pnt->_vnt;
 	_p = _pnt->_prop->param; _ppvar = _pnt->_prop->dparam;
	v = NODEV(_pnt->node);
	return  ( effcai_GB ) - ( theta_d_GB ) ;
}
 
static double _watch2_cond(_pnt) Point_process* _pnt; {
 	double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
	_thread= (Datum*)0; _nt = (_NrnThread*)_pnt->_vnt;
 	_p = _pnt->_prop->param; _ppvar = _pnt->_prop->dparam;
	v = NODEV(_pnt->node);
	return  -( ( effcai_GB ) - ( theta_d_GB ) ) ;
}
 
static double _watch3_cond(_pnt) Point_process* _pnt; {
 	double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
	_thread= (Datum*)0; _nt = (_NrnThread*)_pnt->_vnt;
 	_p = _pnt->_prop->param; _ppvar = _pnt->_prop->dparam;
	v = NODEV(_pnt->node);
	return  ( effcai_GB ) - ( theta_p_GB ) ;
}
 
static double _watch4_cond(_pnt) Point_process* _pnt; {
 	double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
	_thread= (Datum*)0; _nt = (_NrnThread*)_pnt->_vnt;
 	_p = _pnt->_prop->param; _ppvar = _pnt->_prop->dparam;
	v = NODEV(_pnt->node);
	return  -( ( effcai_GB ) - ( theta_p_GB ) ) ;
}
 
static void _net_receive (_pnt, _args, _lflag) Point_process* _pnt; double* _args; double _lflag; 
{  double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   int _watch_rm = 0;
   _thread = (Datum*)0; _nt = (_NrnThread*)_pnt->_vnt;   _p = _pnt->_prop->param; _ppvar = _pnt->_prop->dparam;
  if (_tsav > t){ extern char* hoc_object_name(); hoc_execerror(hoc_object_name(_pnt->ob), ":Event arrived out of order. Must call ParallelContext.set_maxstep AFTER assigning minimum NetCon.delay");}
 _tsav = t;   if (_lflag == 1. ) {*(_tqitem) = 0;}
 {
   double _lp_rec , _lreleased , _ltp , _lfactor ;
 if ( verbose > 0.0 ) {
      printf ( "Time = %g ms, incoming spike at synapse %g\n" , t , synapseID ) ;
      }
   if ( _lflag  == 0.0 ) {
     if ( _args[0] <= 0.0 ) {
       if ( verbose > 0.0 ) {
         printf ( "Inactive synapse, weight = %g\n" , _args[0] ) ;
         }
       }
     else {
       if ( verbose > 0.0 ) {
         printf ( "Flag 0, Regular spike\n" ) ;
         }
       _args[1] = Use_TM + _args[1] * ( 1.0 - Use_TM ) * exp ( - ( t - _args[2] ) / Fac_TM ) ;
       if ( verbose > 0.0 ) {
         printf ( "\tVesicle release probability = %g\n" , _args[1] ) ;
         }
       _lp_rec = 1.0 - exp ( - ( t - _args[2] ) / Dep_TM ) ;
       if ( verbose > 0.0 ) {
         printf ( "\tVesicle recovery probability = %g\n" , _lp_rec ) ;
         }
       if ( verbose > 0.0 ) {
         printf ( "\tVesicle available before recovery = %g\n" , _args[3] ) ;
         }
       _args[3] = _args[3] + brand ( _threadargscomma_ _args[4] , _lp_rec ) ;
       if ( verbose > 0.0 ) {
         printf ( "\tVesicles available after recovery = %g\n" , _args[3] ) ;
         }
       _lreleased = brand ( _threadargscomma_ _args[3] , _args[1] ) ;
       if ( verbose > 0.0 ) {
         printf ( "\tReleased %g vesicles out of %g\n" , _lreleased , _args[3] ) ;
         }
       _ltp = ( tau_r_AMPA * tau_d_AMPA ) / ( tau_d_AMPA - tau_r_AMPA ) * log ( tau_d_AMPA / tau_r_AMPA ) ;
       _lfactor = 1.0 / ( - exp ( - _ltp / tau_r_AMPA ) + exp ( - _ltp / tau_d_AMPA ) ) ;
         if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for general derivimplicit and KINETIC case */
    int __i, __neq = 11;
    double __state = A_AMPA;
    double __primary_delta = (A_AMPA + _lreleased / Nrrp_TM * _lfactor) - __state;
    double __dtsav = dt;
    for (__i = 0; __i < __neq; ++__i) {
      _p[_dlist1[__i]] = 0.0;
    }
    _p[_dlist1[0]] = __primary_delta;
    dt *= 0.5;
    v = NODEV(_pnt->node);
#if NRN_VECTORIZED
    _thread = _nt->_ml_list[_mechtype]->_thread;
#endif
    _ode_matsol_instance1(_threadargs_);
    dt = __dtsav;
    for (__i = 0; __i < __neq; ++__i) {
      _p[_slist1[__i]] += _p[_dlist1[__i]];
    }
  } else {
 A_AMPA = A_AMPA + _lreleased / Nrrp_TM * _lfactor ;
         }
   if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for general derivimplicit and KINETIC case */
    int __i, __neq = 11;
    double __state = B_AMPA;
    double __primary_delta = (B_AMPA + _lreleased / Nrrp_TM * _lfactor) - __state;
    double __dtsav = dt;
    for (__i = 0; __i < __neq; ++__i) {
      _p[_dlist1[__i]] = 0.0;
    }
    _p[_dlist1[1]] = __primary_delta;
    dt *= 0.5;
    v = NODEV(_pnt->node);
#if NRN_VECTORIZED
    _thread = _nt->_ml_list[_mechtype]->_thread;
#endif
    _ode_matsol_instance1(_threadargs_);
    dt = __dtsav;
    for (__i = 0; __i < __neq; ++__i) {
      _p[_slist1[__i]] += _p[_dlist1[__i]];
    }
  } else {
 B_AMPA = B_AMPA + _lreleased / Nrrp_TM * _lfactor ;
         }
 _ltp = ( tau_r_NMDA * tau_d_NMDA ) / ( tau_d_NMDA - tau_r_NMDA ) * log ( tau_d_NMDA / tau_r_NMDA ) ;
       _lfactor = 1.0 / ( - exp ( - _ltp / tau_r_NMDA ) + exp ( - _ltp / tau_d_NMDA ) ) ;
         if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for general derivimplicit and KINETIC case */
    int __i, __neq = 11;
    double __state = A_NMDA;
    double __primary_delta = (A_NMDA + _lreleased / Nrrp_TM * _lfactor) - __state;
    double __dtsav = dt;
    for (__i = 0; __i < __neq; ++__i) {
      _p[_dlist1[__i]] = 0.0;
    }
    _p[_dlist1[3]] = __primary_delta;
    dt *= 0.5;
    v = NODEV(_pnt->node);
#if NRN_VECTORIZED
    _thread = _nt->_ml_list[_mechtype]->_thread;
#endif
    _ode_matsol_instance1(_threadargs_);
    dt = __dtsav;
    for (__i = 0; __i < __neq; ++__i) {
      _p[_slist1[__i]] += _p[_dlist1[__i]];
    }
  } else {
 A_NMDA = A_NMDA + _lreleased / Nrrp_TM * _lfactor ;
         }
   if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for general derivimplicit and KINETIC case */
    int __i, __neq = 11;
    double __state = B_NMDA;
    double __primary_delta = (B_NMDA + _lreleased / Nrrp_TM * _lfactor) - __state;
    double __dtsav = dt;
    for (__i = 0; __i < __neq; ++__i) {
      _p[_dlist1[__i]] = 0.0;
    }
    _p[_dlist1[4]] = __primary_delta;
    dt *= 0.5;
    v = NODEV(_pnt->node);
#if NRN_VECTORIZED
    _thread = _nt->_ml_list[_mechtype]->_thread;
#endif
    _ode_matsol_instance1(_threadargs_);
    dt = __dtsav;
    for (__i = 0; __i < __neq; ++__i) {
      _p[_slist1[__i]] += _p[_dlist1[__i]];
    }
  } else {
 B_NMDA = B_NMDA + _lreleased / Nrrp_TM * _lfactor ;
         }
 _args[3] = _args[3] - _lreleased ;
       _args[4] = Nrrp_TM - _args[3] ;
       if ( verbose > 0.0 ) {
         printf ( "\tFinal vesicle count, Recovered = %g, Unrecovered = %g, Nrrp = %g\n" , _args[3] , _args[4] , Nrrp_TM ) ;
         }
       _args[2] = t ;
       }
     }
   else if ( _lflag  == 1.0 ) {
     if ( verbose > 0.0 ) {
       printf ( "Flag 1, Initialize watchers\n" ) ;
       }
       _nrn_watch_activate(_watch_array, _watch1_cond, 1, _pnt, _watch_rm++, 2.0);
   _nrn_watch_activate(_watch_array, _watch2_cond, 2, _pnt, _watch_rm++, 3.0);
   _nrn_watch_activate(_watch_array, _watch3_cond, 3, _pnt, _watch_rm++, 4.0);
   _nrn_watch_activate(_watch_array, _watch4_cond, 4, _pnt, _watch_rm++, 5.0);
 }
   else if ( _lflag  == 2.0 ) {
     if ( verbose > 0.0 ) {
       printf ( "Flag 2, Activate depression mechanisms\n" ) ;
       }
     dep_GB = 1.0 ;
     }
   else if ( _lflag  == 3.0 ) {
     if ( verbose > 0.0 ) {
       printf ( "Flag 3, Deactivate depression mechanisms\n" ) ;
       }
     dep_GB = 0.0 ;
     }
   else if ( _lflag  == 4.0 ) {
     if ( verbose > 0.0 ) {
       printf ( "Flag 4, Activate potentiation mechanisms\n" ) ;
       }
     pot_GB = 1.0 ;
     }
   else if ( _lflag  == 5.0 ) {
     if ( verbose > 0.0 ) {
       printf ( "Flag 5, Deactivate potentiation mechanisms\n" ) ;
       }
     pot_GB = 0.0 ;
     }
   } }
 
static void _net_init(Point_process* _pnt, double* _args, double _lflag) {
       double* _p = _pnt->_prop->param;
    Datum* _ppvar = _pnt->_prop->dparam;
    Datum* _thread = (Datum*)0;
    _NrnThread* _nt = (_NrnThread*)_pnt->_vnt;
 _args[0] = 1.0 ;
   _args[1] = 0.0 ;
   _args[2] = 0.0 ;
   _args[3] = Nrrp_TM ;
   _args[4] = 0.0 ;
   }
 
double nernst ( _threadargsprotocomma_ double _lci , double _lco , double _lz ) {
   double _lnernst;
 _lnernst = ( 1000.0 ) * R * ( celsius + 273.15 ) / ( _lz * FARADAY ) * log ( _lco / _lci ) ;
   if ( verbose > 1.0 ) {
      printf ( "nernst:%g R:%g temperature (c):%g \n" , _lnernst , R , celsius ) ;
      }
   
return _lnernst;
 }
 
static double _hoc_nernst(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (_NrnThread*)((Point_process*)_vptr)->_vnt;
 _r =  nernst ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) , *getarg(3) );
 return(_r);
}
 
static int  setRNG ( _threadargsproto_ ) {
   
/*VERBATIM*/
    #ifndef CORENEURON_BUILD
    // For compatibility, allow for either MCellRan4 or Random123
    // Distinguish by the arg types
    // Object => MCellRan4, seeds (double) => Random123
    usingR123 = 0;
    if( ifarg(1) && hoc_is_double_arg(1) ) {
        nrnran123_State** pv = (nrnran123_State**)(&_p_rng_TM);
        uint32_t a2 = 0;
        uint32_t a3 = 0;
        if (*pv) {
            nrnran123_deletestream(*pv);
            *pv = (nrnran123_State*)0;
        }
        if (ifarg(2)) {
            a2 = (uint32_t)*getarg(2);
        }
        if (ifarg(3)) {
            a3 = (uint32_t)*getarg(3);
        }
        *pv = nrnran123_newstream3((uint32_t)*getarg(1), a2, a3);
        usingR123 = 1;
    } else if( ifarg(1) ) {   // not a double, so assume hoc object type
        void** pv = (void**)(&_p_rng_TM);
        *pv = nrn_random_arg(1);
    } else {  // no arg, so clear pointer
        void** pv = (void**)(&_p_rng_TM);
        *pv = (void*)0;
    }
    #endif
  return 0; }
 
static double _hoc_setRNG(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (_NrnThread*)((Point_process*)_vptr)->_vnt;
 _r = 1.;
 setRNG ( _p, _ppvar, _thread, _nt );
 return(_r);
}
 
double urand ( _threadargsproto_ ) {
   double _lurand;
 
/*VERBATIM*/
    double value;
    if ( usingR123 ) {
        value = nrnran123_dblpick((nrnran123_State*)_p_rng_TM);
    } else if (_p_rng_TM) {
        #ifndef CORENEURON_BUILD
        value = nrn_random_pick(_p_rng_TM);
        #endif
    } else {
        value = 0.0;
    }
    _lurand = value;
 
return _lurand;
 }
 
static double _hoc_urand(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (_NrnThread*)((Point_process*)_vptr)->_vnt;
 _r =  urand ( _p, _ppvar, _thread, _nt );
 return(_r);
}
 
double brand ( _threadargsprotocomma_ double _ln , double _lp ) {
   double _lbrand;
 double _lresult , _lcount , _lsuccess ;
 _lsuccess = 0.0 ;
   {int  _lcount ;for ( _lcount = 0 ; _lcount <= ( ((int) _ln ) - 1 ) ; _lcount ++ ) {
     _lresult = urand ( _threadargs_ ) ;
     if ( _lresult <= _lp ) {
       _lsuccess = _lsuccess + 1.0 ;
       }
     } }
   _lbrand = _lsuccess ;
   
return _lbrand;
 }
 
static double _hoc_brand(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (_NrnThread*)((Point_process*)_vptr)->_vnt;
 _r =  brand ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) );
 return(_r);
}
 
double bbsavestate ( _threadargsproto_ ) {
   double _lbbsavestate;
 _lbbsavestate = 0.0 ;
   
/*VERBATIM*/
    #ifndef CORENEURON_BUILD
        /* first arg is direction (0 save, 1 restore), second is array*/
        /* if first arg is -1, fill xdir with the size of the array */
        double *xdir, *xval, *hoc_pgetarg();
        long nrn_get_random_sequence(void* r);
        void nrn_set_random_sequence(void* r, int val);
        xdir = hoc_pgetarg(1);
        xval = hoc_pgetarg(2);
        if (_p_rng_TM) {
            // tell how many items need saving
            if (*xdir == -1) {  // count items
                if( usingR123 ) {
                    *xdir = 2.0;
                } else {
                    *xdir = 1.0;
                }
                return 0.0;
            } else if(*xdir ==0 ) {  // save
                if( usingR123 ) {
                    uint32_t seq;
                    unsigned char which;
                    nrnran123_getseq( (nrnran123_State*)_p_rng_TM, &seq, &which );
                    xval[0] = (double) seq;
                    xval[1] = (double) which;
                } else {
                    xval[0] = (double)nrn_get_random_sequence(_p_rng_TM);
                }
            } else {  // restore
                if( usingR123 ) {
                    nrnran123_setseq( (nrnran123_State*)_p_rng_TM, (uint32_t)xval[0], (char)xval[1] );
                } else {
                    nrn_set_random_sequence(_p_rng_TM, (long)(xval[0]));
                }
            }
        }
    #endif
 
return _lbbsavestate;
 }
 
static double _hoc_bbsavestate(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (_NrnThread*)((Point_process*)_vptr)->_vnt;
 _r =  bbsavestate ( _p, _ppvar, _thread, _nt );
 return(_r);
}
 
/*VERBATIM*/
static void bbcore_write(double* dArray, int* iArray, int* doffset, int* ioffset, _threadargsproto_) {
    // make sure offset array non-null
    if (iArray) {
        // get handle to random123 instance
        nrnran123_State** pv = (nrnran123_State**)(&_p_rng_TM);
        // get location for storing ids
        uint32_t* ia = ((uint32_t*)iArray) + *ioffset;
        // retrieve/store identifier seeds
        nrnran123_getids3(*pv, ia, ia+1, ia+2);
        // retrieve/store stream sequence
        unsigned char which;
        nrnran123_getseq(*pv, ia+3, &which);
        ia[4] = (int)which;
    }

    // increment integer offset (2 identifier), no double data
    *ioffset += 5;
    *doffset += 0;
}

static void bbcore_read(double* dArray, int* iArray, int* doffset, int* ioffset, _threadargsproto_) {
    // make sure it's not previously set
    assert(!_p_rng_TM);
    uint32_t* ia = ((uint32_t*)iArray) + *ioffset;
    // make sure non-zero identifier seeds
    if (ia[0] != 0 || ia[1] != 0 || ia[2] != 0) {
        nrnran123_State** pv = (nrnran123_State**)(&_p_rng_TM);
        // get new stream
        *pv = nrnran123_newstream3(ia[0], ia[1], ia[2]);
        // restore sequence
        nrnran123_setseq(*pv, ia[3], (char)ia[4]);
    }
    // increment intger offset (2 identifiers), no double data
    *ioffset += 5;
    *doffset += 0;
}
 
static int _ode_count(int _type){ return 11;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
     _ode_spec1 (_p, _ppvar, _thread, _nt);
 }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 11; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
 _ode_matsol_instance1(_threadargs_);
 }}
 
static void _thread_mem_init(Datum* _thread) {
   _thread[_dith1]._pval = (double*)ecalloc(22, sizeof(double));
   _newtonspace1 = nrn_cons_newtonspace(11);
 }
 
static void _thread_cleanup(Datum* _thread) {
   free((void*)(_thread[_dith1]._pval));
   nrn_destroy_newtonspace(_newtonspace1);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  A_NMDA = A_NMDA0;
  A_AMPA = A_AMPA0;
  B_NMDA = B_NMDA0;
  B_AMPA = B_AMPA0;
  Use_TM = Use_TM0;
  cai_CR = cai_CR0;
  effcai_GB = effcai_GB0;
  gmax_AMPA = gmax_AMPA0;
  h_VDCC = h_VDCC0;
  m_VDCC = m_VDCC0;
  rho_GB = rho_GB0;
 {
   A_AMPA = 0.0 ;
   B_AMPA = 0.0 ;
   gmax_AMPA = gmax0_AMPA ;
   A_NMDA = 0.0 ;
   B_NMDA = 0.0 ;
   Use_TM = Use0_TM ;
   cai_CR = min_ca_CR ;
   rho_GB = rho0_GB ;
   effcai_GB = 0.0 ;
   dep_GB = 0.0 ;
   pot_GB = 0.0 ;
   net_send ( _tqitem, (double*)0, _ppvar[1]._pvoid, t +  0.0 , 1.0 ) ;
   }
 
}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _tsav = -1e20;
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
 initmodel(_p, _ppvar, _thread, _nt);
}
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   double _lEca_syn , _lmggate , _li_AMPA , _li_NMDA , _lPf_NMDA , _lgca_bar_abs_VDCC , _lgca_VDCC ;
 g_AMPA = ( 1e-3 ) * gmax_AMPA * ( B_AMPA - A_AMPA ) ;
   _li_AMPA = g_AMPA * ( v - E_AMPA ) ;
   _lmggate = 1.0 / ( 1.0 + exp ( - slope_NMDA * v ) * ( mgo_NMDA / scale_NMDA ) ) ;
   g_NMDA = ( 1e-3 ) * gmax_NMDA * _lmggate * ( B_NMDA - A_NMDA ) ;
   _li_NMDA = g_NMDA * ( v - E_NMDA ) ;
   _lPf_NMDA = ( 4.0 * cao_CR ) / ( 4.0 * cao_CR + ( 1.0 / 1.38 ) * 120.0 ) * 0.6 ;
   ica_NMDA = _lPf_NMDA * g_NMDA * ( v - 40.0 ) ;
   _lgca_bar_abs_VDCC = gca_bar_VDCC * 4.0 * PI * pow( ( 3.0 / 4.0 * volume_CR * 1.0 / PI ) , ( 2.0 / 3.0 ) ) ;
   _lgca_VDCC = ( 1e-3 ) * _lgca_bar_abs_VDCC * m_VDCC * m_VDCC * h_VDCC ;
   _lEca_syn = nernst ( _threadargscomma_ cai_CR , cao_CR , 2.0 ) ;
   ica_VDCC = _lgca_VDCC * ( v - _lEca_syn ) ;
   vsyn = v ;
   i = _li_AMPA + _li_NMDA + ica_VDCC ;
   }
 _current += i;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
 	}
 _g = (_g - _rhs)/.001;
 _g *=  1.e2/(_nd_area);
 _rhs *= 1.e2/(_nd_area);
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}
 
}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}
 
}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
double _dtsav = dt;
if (secondorder) { dt *= 0.5; }
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
 {  _deriv1_advance = 1;
 derivimplicit_thread(11, _slist1, _dlist1, _p, state, _ppvar, _thread, _nt);
_deriv1_advance = 0;
     if (secondorder) {
    int _i;
    for (_i = 0; _i < 11; ++_i) {
      _p[_slist1[_i]] += dt*_p[_dlist1[_i]];
    }}
 }}}
 dt = _dtsav;
}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(A_AMPA) - _p;  _dlist1[0] = &(DA_AMPA) - _p;
 _slist1[1] = &(B_AMPA) - _p;  _dlist1[1] = &(DB_AMPA) - _p;
 _slist1[2] = &(gmax_AMPA) - _p;  _dlist1[2] = &(Dgmax_AMPA) - _p;
 _slist1[3] = &(A_NMDA) - _p;  _dlist1[3] = &(DA_NMDA) - _p;
 _slist1[4] = &(B_NMDA) - _p;  _dlist1[4] = &(DB_NMDA) - _p;
 _slist1[5] = &(Use_TM) - _p;  _dlist1[5] = &(DUse_TM) - _p;
 _slist1[6] = &(m_VDCC) - _p;  _dlist1[6] = &(Dm_VDCC) - _p;
 _slist1[7] = &(h_VDCC) - _p;  _dlist1[7] = &(Dh_VDCC) - _p;
 _slist1[8] = &(cai_CR) - _p;  _dlist1[8] = &(Dcai_CR) - _p;
 _slist1[9] = &(effcai_GB) - _p;  _dlist1[9] = &(Deffcai_GB) - _p;
 _slist1[10] = &(rho_GB) - _p;  _dlist1[10] = &(Drho_GB) - _p;
 _slist2[0] = &(A_NMDA) - _p;
 _slist2[1] = &(A_AMPA) - _p;
 _slist2[2] = &(B_NMDA) - _p;
 _slist2[3] = &(B_AMPA) - _p;
 _slist2[4] = &(Use_TM) - _p;
 _slist2[5] = &(cai_CR) - _p;
 _slist2[6] = &(effcai_GB) - _p;
 _slist2[7] = &(gmax_AMPA) - _p;
 _slist2[8] = &(h_VDCC) - _p;
 _slist2[9] = &(m_VDCC) - _p;
 _slist2[10] = &(rho_GB) - _p;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "GluSynapse_andras.mod";
static const char* nmodl_file_text = 
  "COMMENT\n"
  "/**\n"
  " * @file GluSynapse.mod\n"
  " * @brief Probabilistic synapse featuring long-term plasticity\n"
  " * @author king, chindemi, rossert\n"
  " * @date 2019-09-20\n"
  " * @version 1.0.0\n"
  " * @remark Copyright BBP/EPFL 2005-2019; All rights reserved.\n"
  "           Do not distribute without further notice.\n"
  " */\n"
  " Glutamatergic synapse model featuring:\n"
  "1) AMPA receptor with a dual-exponential conductance profile.\n"
  "2) NMDA receptor  with a dual-exponential conductance profile and magnesium\n"
  "   block as described in Jahr and Stevens 1990.\n"
  "3) Tsodyks-Markram presynaptic short-term plasticity as Barros et al. 2019.\n"
  "   Implementation based on the work of Eilif Muller, Michael Reimann and\n"
  "   Srikanth Ramaswamy (Blue Brain Project, August 2011), who introduced the\n"
  "   2-state Markov model of vesicle release. The new model is an extension of\n"
  "   Fuhrmann et al. 2002, motivated by the following constraints:\n"
  "        a) No consumption on failure\n"
  "        b) No release until recovery\n"
  "        c) Same ensemble averaged trace as canonical Tsodyks-Markram using same\n"
  "           parameters determined from experiment.\n"
  "   For a pre-synaptic spike or external spontaneous release trigger event, the\n"
  "   synapse will only release if it is in the recovered state, and with\n"
  "   probability u (which follows facilitation dynamics). If it releases, it will\n"
  "   transition to the unrecovered state. Recovery is as a Poisson process with\n"
  "   rate 1/Dep.\n"
  "   John Rahmon and Giuseppe Chindemi introduced multi-vesicular release as an\n"
  "   extension of the 2-state Markov model of vesicle release described above\n"
  "   (Blue Brain Project, February 2017).\n"
  "4) NMDAR-mediated calcium current. Fractional calcium current Pf_NMDA from\n"
  "   Schneggenburger et al. 1993. Fractional NMDAR conductance treated as a\n"
  "   calcium-only permeable channel with Erev = 40 mV independent of extracellular\n"
  "   calcium concentration (see Jahr and Stevens 1993). Implemented by Christian\n"
  "   Rossert and Giuseppe Chindemi (Blue Brain Project, 2016).\n"
  "5) Spine volume.\n"
  "6) VDCC.\n"
  "7) Postsynaptic calcium dynamics.\n"
  "8) Long-term synaptic plasticity. Calcium-based STDP model based on Graupner and\n"
  "   Brunel 2012.\n"
  "Model implementation, optimization and simulation curated by James King (Blue\n"
  "Brain Project, 2019).\n"
  "ENDCOMMENT\n"
  "\n"
  "TITLE Glutamatergic synapse\n"
  "\n"
  "NEURON {\n"
  "    THREADSAFE\n"
  "    POINT_PROCESS GluSynapse_andras\n"
  "    : AMPA Receptor\n"
  "    GLOBAL tau_r_AMPA, E_AMPA\n"
  "    RANGE tau_d_AMPA, gmax0_AMPA, gmax_d_AMPA, gmax_p_AMPA, g_AMPA\n"
  "    : NMDA Receptor\n"
  "    GLOBAL mgo_NMDA, scale_NMDA, slope_NMDA\n"
  "    GLOBAL tau_r_NMDA, tau_d_NMDA, E_NMDA\n"
  "    RANGE gmax_NMDA, g_NMDA\n"
  "    : Stochastic Tsodyks-Markram Multi-Vesicular Release\n"
  "    RANGE Use0_TM, Dep_TM, Fac_TM, Nrrp_TM\n"
  "    RANGE Use_d_TM, Use_p_TM\n"
  "    BBCOREPOINTER rng_TM\n"
  "    : NMDAR-mediated calcium current\n"
  "    RANGE ica_NMDA\n"
  "    : Spine\n"
  "    RANGE volume_CR\n"
  "    : VDCC (R-type)\n"
  "    GLOBAL ljp_VDCC, vhm_VDCC, km_VDCC, mtau_VDCC, vhh_VDCC, kh_VDCC, htau_VDCC, gca_bar_VDCC\n"
  "    RANGE ica_VDCC\n"
  "    : Postsynaptic Ca2+ dynamics\n"
  "    GLOBAL gamma_ca_CR, tau_ca_CR, min_ca_CR, cao_CR\n"
  "    : Long-term synaptic plasticity\n"
  "    GLOBAL rho_star_GB, tau_ind_GB, tau_exp_GB, tau_effca_GB\n"
  "    GLOBAL gamma_d_GB, gamma_p_GB\n"
  "    RANGE theta_d_GB, theta_p_GB, rho0_GB, dep_GB, pot_GB\n"
  "    : Misc\n"
  "    RANGE vsyn, synapseID, selected_for_report, verbose\n"
  "    NONSPECIFIC_CURRENT i\n"
  "}\n"
  "\n"
  "\n"
  "UNITS {\n"
  "    (nA)    = (nanoamp)\n"
  "    (mV)    = (millivolt)\n"
  "    (uS)    = (microsiemens)\n"
  "    (nS)    = (nanosiemens)\n"
  "    (pS)    = (picosiemens)\n"
  "    (umho)  = (micromho)\n"
  "    (um)    = (micrometers)\n"
  "    (mM)    = (milli/liter)\n"
  "    (uM)    = (micro/liter)\n"
  "    FARADAY = (faraday) (coulomb)\n"
  "    PI      = (pi)      (1)\n"
  "    R       = (k-mole)  (joule/degC)\n"
  "}\n"
  "\n"
  "\n"
  "PARAMETER {\n"
  "    celsius                     (degC)\n"
  "    : AMPA Receptor\n"
  "    tau_r_AMPA      = 0.2       (ms)        : Tau rise, dual-exponential conductance profile\n"
  "    tau_d_AMPA      = 1.7       (ms)        : Tau decay, IMPORTANT: tau_r < tau_d\n"
  "    E_AMPA          = 0         (mV)        : Reversal potential\n"
  "    gmax0_AMPA      = 1.0       (nS)        : Initial peak conductance\n"
  "    gmax_d_AMPA     = 1.0       (nS)        : Peak conductance in the depressed state\n"
  "    gmax_p_AMPA     = 2.0       (nS)        : Peak conductance in the potentitated state\n"
  "    : NMDA Receptor\n"
  "    mgo_NMDA        = 1         (mM)        : Extracellular magnesium concentration\n"
  "    scale_NMDA      = 2.552     (mM)        : Scale of the mg block (Vargas-Caballero and Robinson 2003)\n"
  "    slope_NMDA      = 0.072     (/mV)       : Slope of the ma block (Vargas-Caballero and Robinson 2003)\n"
  "    tau_r_NMDA      = 0.29      (ms)        : Tau rise, dual-exponential conductance profile\n"
  "    tau_d_NMDA      = 70        (ms)        : Tau decay, IMPORTANT: tau_r < tau_d\n"
  "    E_NMDA          = -3        (mV)        : Reversal potential (Vargas-Caballero and Robinson 2003)\n"
  "    gmax_NMDA       = 0.55      (nS)        : Peak conductance\n"
  "    : Stochastic Tsodyks-Markram Multi-Vesicular Release\n"
  "    Use0_TM         = 0.5       (1)         : Initial utilization of synaptic efficacy\n"
  "    Dep_TM          = 100       (ms)        : Relaxation time constant from depression\n"
  "    Fac_TM          = 10        (ms)        : Relaxation time constant from facilitation\n"
  "    Nrrp_TM         = 1         (1)         : Number of release sites for given contact\n"
  "    Use_d_TM        = 0.2       (1)         : Depressed Use\n"
  "    Use_p_TM        = 0.8       (1)         : Potentiated Use\n"
  "    : Spine\n"
  "    volume_CR       = 0.087     (um3)       : From spine data by Ruth Benavides-Piccione (unpublished)\n"
  "    : VDCC (R-type)\n"
  "    gca_bar_VDCC    = 0.0744    (nS/um2)    : Density spines: 20 um-2 (Sabatini 2000), unitary conductance VGCC 3.72 pS (Bartol 2015)\n"
  "    ljp_VDCC        = 0         (mV)\n"
  "    vhm_VDCC        = -5.9      (mV)        : v 1/2 for act, Magee and Johnston 1995 (corrected for m*m)\n"
  "    km_VDCC         = 9.5       (mV)        : act slope, Magee and Johnston 1995 (corrected for m*m)\n"
  "    vhh_VDCC        = -39       (mV)        : v 1/2 for inact, Magee and Johnston 1995\n"
  "    kh_VDCC         = -9.2      (mV)        : inact, Magee and Johnston 1995\n"
  "    mtau_VDCC       = 1         (ms)        : max time constant (guess)\n"
  "    htau_VDCC       = 27        (ms)        : max time constant 100*0.27\n"
  "    : Postsynaptic Ca2+ dynamics\n"
  "    gamma_ca_CR     = 0.04      (1)         : Percent of free calcium (not buffered), Sabatini et al 2002: kappa_e = 24+-11 (also 14 (2-31) or 22 (18-33))\n"
  "    tau_ca_CR       = 12        (ms)        : Rate of removal of calcium, Sabatini et al 2002: 14ms (12-20ms)\n"
  "    min_ca_CR       = 70e-6     (mM)        : Sabatini et al 2002: 70+-29 nM, per AP: 1.1 (0.6-8.2) uM = 1100 e-6 mM = 1100 nM\n"
  "    cao_CR          = 2.0       (mM)        : Extracellular calcium concentration in slices\n"
  "    : Long-term synaptic plasticity\n"
  "    rho_star_GB     = 0.5       (1)\n"
  "    tau_ind_GB      = 70        (s)\n"
  "    tau_exp_GB      = 100       (s)\n"
  "    tau_effca_GB    = 200       (ms)\n"
  "    gamma_d_GB      = 100       (1)\n"
  "    gamma_p_GB      = 450       (1)\n"
  "    theta_d_GB      = 0.006     (us/liter)\n"
  "    theta_p_GB      = 0.001     (us/liter)\n"
  "    rho0_GB         = 0         (1)\n"
  "    : Misc\n"
  "    synapseID       = 0\n"
  "    verbose         = 0\n"
  "    selected_for_report = 0\n"
  "}\n"
  "\n"
  "VERBATIM\n"
  "/**\n"
  " * This Verbatim block is needed to generate random numbers from a uniform\n"
  " * distribution U(0, 1).\n"
  " */\n"
  "#include <stdlib.h>\n"
  "#include <stdio.h>\n"
  "#include <math.h>\n"
  "#include \"nrnran123.h\"\n"
  "double nrn_random_pick(void* r);\n"
  "void* nrn_random_arg(int argpos);\n"
  "ENDVERBATIM\n"
  "\n"
  "ASSIGNED {\n"
  "    : AMPA Receptor\n"
  "    g_AMPA          (uS)\n"
  "    : NMDA Receptor\n"
  "    g_NMDA          (uS)\n"
  "    : Stochastic Tsodyks-Markram Multi-Vesicular Release\n"
  "    rng_TM                  : Random Number Generator\n"
  "    usingR123               : TEMPORARY until mcellran4 completely deprecated\n"
  "    : NMDAR-mediated calcium current\n"
  "    ica_NMDA        (nA)\n"
  "    : VDCC (R-type)\n"
  "    ica_VDCC        (nA)\n"
  "    : Long-term synaptic plasticity\n"
  "    dep_GB          (1)\n"
  "    pot_GB          (1)\n"
  "    : Misc\n"
  "    v               (mV)\n"
  "    vsyn            (mV)\n"
  "    i               (nA)\n"
  "}\n"
  "\n"
  "STATE {\n"
  "    : AMPA Receptor\n"
  "    A_AMPA      (1)\n"
  "    B_AMPA      (1)\n"
  "    gmax_AMPA   (nS)\n"
  "    : NMDA Receptor\n"
  "    A_NMDA      (1)\n"
  "    B_NMDA      (1)\n"
  "    : Stochastic Tsodyks-Markram Multi-Vesicular Release\n"
  "    Use_TM      (1)\n"
  "    : VDCC (R-type)\n"
  "    m_VDCC      (1)\n"
  "    h_VDCC      (1)\n"
  "    : Postsynaptic Ca2+ dynamics\n"
  "    cai_CR      (mM)        <1e-6>\n"
  "    : Long-term synaptic plasticity\n"
  "    rho_GB      (1)\n"
  "    effcai_GB   (us/liter)  <1e-3>\n"
  "}\n"
  "\n"
  "INITIAL{\n"
  "    : AMPA Receptor\n"
  "    A_AMPA      = 0\n"
  "    B_AMPA      = 0\n"
  "    gmax_AMPA   = gmax0_AMPA\n"
  "    : NMDA Receptor\n"
  "    A_NMDA      = 0\n"
  "    B_NMDA      = 0\n"
  "    : Stochastic Tsodyks-Markram Multi-Vesicular Release\n"
  "    Use_TM      = Use0_TM\n"
  "    : Postsynaptic Ca2+ dynamics\n"
  "    cai_CR      = min_ca_CR\n"
  "    : Long-term synaptic plasticity\n"
  "    rho_GB      = rho0_GB\n"
  "    effcai_GB   = 0\n"
  "    dep_GB      = 0\n"
  "    pot_GB      = 0\n"
  "    : Initialize watchers\n"
  "    net_send(0, 1)\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "    LOCAL Eca_syn, mggate, i_AMPA, i_NMDA, Pf_NMDA, gca_bar_abs_VDCC, gca_VDCC\n"
  "    SOLVE state METHOD derivimplicit\n"
  "    : AMPA Receptor\n"
  "    g_AMPA = (1e-3)*gmax_AMPA*(B_AMPA - A_AMPA)\n"
  "    i_AMPA = g_AMPA*(v-E_AMPA)\n"
  "    : NMDA Receptor\n"
  "    mggate = 1 / (1 + exp(-slope_NMDA*v) * (mgo_NMDA/scale_NMDA))\n"
  "    g_NMDA = (1e-3)*gmax_NMDA*mggate*(B_NMDA - A_NMDA)\n"
  "    i_NMDA = g_NMDA*(v - E_NMDA)\n"
  "    : NMDAR-mediated calcium current\n"
  "    Pf_NMDA  = (4*cao_CR) / (4*cao_CR + (1/1.38) * 120 (mM)) * 0.6\n"
  "    ica_NMDA = Pf_NMDA*g_NMDA*(v-40.0)\n"
  "    : VDCC (R-type), assuming sphere for spine head\n"
  "    gca_bar_abs_VDCC = gca_bar_VDCC * 4(um2)*PI*(3(1/um3)/4*volume_CR*1/PI)^(2/3)\n"
  "    gca_VDCC = (1e-3) * gca_bar_abs_VDCC * m_VDCC * m_VDCC * h_VDCC\n"
  "    Eca_syn = nernst(cai_CR, cao_CR, 2)\n"
  "    ica_VDCC = gca_VDCC*(v-Eca_syn)\n"
  "    : Update synaptic voltage (for recording convenience)\n"
  "    vsyn = v\n"
  "    : Update current\n"
  "    i = i_AMPA + i_NMDA + ica_VDCC\n"
  "}\n"
  "\n"
  "DERIVATIVE state {\n"
  "    LOCAL minf_VDCC, hinf_VDCC\n"
  "    : AMPA Receptor\n"
  "    A_AMPA'      = - A_AMPA/tau_r_AMPA\n"
  "    B_AMPA'      = - B_AMPA/tau_d_AMPA\n"
  "    gmax_AMPA'   = (gmax_d_AMPA + rho_GB*(gmax_p_AMPA - gmax_d_AMPA) - gmax_AMPA) / ((1e3)*tau_exp_GB)\n"
  "    : NMDA Receptor\n"
  "    A_NMDA'      = - A_NMDA/tau_r_NMDA\n"
  "    B_NMDA'      = - B_NMDA/tau_d_NMDA\n"
  "    : Stochastic Tsodyks-Markram Multi-Vesicular Release\n"
  "    Use_TM'      = (Use_d_TM + rho_GB*(Use_p_TM - Use_d_TM) - Use_TM) / ((1e3)*tau_exp_GB)\n"
  "    : VDCC (R-type)\n"
  "    minf_VDCC    = 1 / (1 + exp(((vhm_VDCC - ljp_VDCC) - v) / km_VDCC))\n"
  "    hinf_VDCC    = 1 / (1 + exp(((vhh_VDCC - ljp_VDCC) - v) / kh_VDCC))\n"
  "    m_VDCC'      = (minf_VDCC-m_VDCC)/mtau_VDCC\n"
  "    h_VDCC'      = (hinf_VDCC-h_VDCC)/htau_VDCC\n"
  "    : Postsynaptic Ca2+ dynamics\n"
  "    cai_CR'      = - (1e-9)*(ica_NMDA + ica_VDCC)*gamma_ca_CR/((1e-15)*volume_CR*2*FARADAY)\n"
  "                   - (cai_CR - min_ca_CR)/tau_ca_CR\n"
  "    : Long-term synaptic plasticity\n"
  "    effcai_GB'   = - effcai_GB/tau_effca_GB + (cai_CR - min_ca_CR)\n"
  "    rho_GB'      = ( - rho_GB*(1 - rho_GB)*(rho_star_GB - rho_GB)\n"
  "                     + pot_GB*gamma_p_GB*(1 - rho_GB)\n"
  "                     - dep_GB*gamma_d_GB*rho_GB ) / ((1e3)*tau_ind_GB)\n"
  "}\n"
  "\n"
  "NET_RECEIVE (weight, u, tsyn (ms), recovered, unrecovered) {\n"
  "    LOCAL p_rec, released, tp, factor\n"
  "    INITIAL {\n"
  "        weight = 1\n"
  "        u = 0\n"
  "        tsyn = 0 (ms)\n"
  "        recovered = Nrrp_TM\n"
  "        unrecovered = 0\n"
  "    }\n"
  "    if(verbose > 0){ UNITSOFF printf(\"Time = %g ms, incoming spike at synapse %g\\n\", t, synapseID) UNITSON }\n"
  "    if(flag == 0) {\n"
  "        if(weight <= 0){\n"
  "            : Do not perform any calculations if the synapse (netcon) is deactivated.\n"
  "            : This avoids drawing from the random stream\n"
  "            : WARNING In this model *weight* is only used to activate/deactivate the\n"
  "            :         synapse. The conductance is stored in gmax_AMPA and gmax_NMDA.\n"
  "            if(verbose > 0){ printf(\"Inactive synapse, weight = %g\\n\", weight) }\n"
  "        } else {\n"
  "            : Flag 0: Regular spike\n"
  "            if(verbose > 0){ printf(\"Flag 0, Regular spike\\n\") }\n"
  "            : Update facilitation variable as Eq. 2 in Fuhrmann et al. 2002\n"
  "            u = Use_TM + u*(1 - Use_TM)*exp(-(t - tsyn)/Fac_TM)\n"
  "            if ( verbose > 0 ) { printf(\"\\tVesicle release probability = %g\\n\", u) }\n"
  "            : Recovery\n"
  "            p_rec = 1 - exp(-(t - tsyn)/Dep_TM)\n"
  "            if ( verbose > 0 ) { printf(\"\\tVesicle recovery probability = %g\\n\", p_rec) }\n"
  "            if ( verbose > 0 ) { printf(\"\\tVesicle available before recovery = %g\\n\", recovered) }\n"
  "            recovered = recovered + brand(unrecovered, p_rec)\n"
  "            if ( verbose > 0 ) { printf(\"\\tVesicles available after recovery = %g\\n\", recovered) }\n"
  "            : Release\n"
  "            released = brand(recovered, u)\n"
  "            if ( verbose > 0 ) { printf(\"\\tReleased %g vesicles out of %g\\n\", released, recovered) }\n"
  "            : Update AMPA variables\n"
  "            tp = (tau_r_AMPA*tau_d_AMPA)/(tau_d_AMPA-tau_r_AMPA)*log(tau_d_AMPA/tau_r_AMPA)  : Time to peak\n"
  "            factor = 1 / (-exp(-tp/tau_r_AMPA)+exp(-tp/tau_d_AMPA))  : Normalization factor\n"
  "            A_AMPA = A_AMPA + released/Nrrp_TM*factor\n"
  "            B_AMPA = B_AMPA + released/Nrrp_TM*factor\n"
  "            : Update NMDA variables\n"
  "            tp = (tau_r_NMDA*tau_d_NMDA)/(tau_d_NMDA-tau_r_NMDA)*log(tau_d_NMDA/tau_r_NMDA)  : Time to peak\n"
  "            factor = 1 / (-exp(-tp/tau_r_NMDA)+exp(-tp/tau_d_NMDA))  : Normalization factor\n"
  "            A_NMDA = A_NMDA + released/Nrrp_TM*factor\n"
  "            B_NMDA = B_NMDA + released/Nrrp_TM*factor\n"
  "            : Update vesicle pool\n"
  "            recovered = recovered - released\n"
  "            unrecovered = Nrrp_TM - recovered\n"
  "            if ( verbose > 0 ) { printf(\"\\tFinal vesicle count, Recovered = %g, Unrecovered = %g, Nrrp = %g\\n\", recovered, unrecovered, Nrrp_TM) }\n"
  "            : Update tsyn\n"
  "            : tsyn knows about all spikes, not only those that released\n"
  "            : i.e. each spike can increase the u, regardless of recovered state\n"
  "            :      and each spike trigger an evaluation of recovery\n"
  "            tsyn = t\n"
  "        }\n"
  "    } else if(flag == 1) {\n"
  "        : Flag 1, Initialize watchers\n"
  "        if(verbose > 0){ printf(\"Flag 1, Initialize watchers\\n\") }\n"
  "        WATCH (effcai_GB > theta_d_GB) 2\n"
  "        WATCH (effcai_GB < theta_d_GB) 3\n"
  "        WATCH (effcai_GB > theta_p_GB) 4\n"
  "        WATCH (effcai_GB < theta_p_GB) 5\n"
  "    } else if(flag == 2) {\n"
  "        : Flag 2, Activate depression mechanisms\n"
  "        if(verbose > 0){ printf(\"Flag 2, Activate depression mechanisms\\n\") }\n"
  "        dep_GB = 1\n"
  "    } else if(flag == 3) {\n"
  "        : Flag 3, Deactivate depression mechanisms\n"
  "        if(verbose > 0){ printf(\"Flag 3, Deactivate depression mechanisms\\n\") }\n"
  "        dep_GB = 0\n"
  "    } else if(flag == 4) {\n"
  "        : Flag 4, Activate potentiation mechanisms\n"
  "        if(verbose > 0){ printf(\"Flag 4, Activate potentiation mechanisms\\n\") }\n"
  "        pot_GB = 1\n"
  "    } else if(flag == 5) {\n"
  "        : Flag 5, Deactivate potentiation mechanisms\n"
  "        if(verbose > 0){ printf(\"Flag 5, Deactivate potentiation mechanisms\\n\") }\n"
  "        pot_GB = 0\n"
  "    }\n"
  "}\n"
  "\n"
  "FUNCTION nernst(ci(mM), co(mM), z) (mV) {\n"
  "    nernst = (1000) * R * (celsius + 273.15) / (z*FARADAY) * log(co/ci)\n"
  "    if(verbose > 1) { UNITSOFF printf(\"nernst:%g R:%g temperature (c):%g \\n\", nernst, R, celsius) UNITSON }\n"
  "}\n"
  "\n"
  "PROCEDURE setRNG() {\n"
  "    VERBATIM\n"
  "    #ifndef CORENEURON_BUILD\n"
  "    // For compatibility, allow for either MCellRan4 or Random123\n"
  "    // Distinguish by the arg types\n"
  "    // Object => MCellRan4, seeds (double) => Random123\n"
  "    usingR123 = 0;\n"
  "    if( ifarg(1) && hoc_is_double_arg(1) ) {\n"
  "        nrnran123_State** pv = (nrnran123_State**)(&_p_rng_TM);\n"
  "        uint32_t a2 = 0;\n"
  "        uint32_t a3 = 0;\n"
  "        if (*pv) {\n"
  "            nrnran123_deletestream(*pv);\n"
  "            *pv = (nrnran123_State*)0;\n"
  "        }\n"
  "        if (ifarg(2)) {\n"
  "            a2 = (uint32_t)*getarg(2);\n"
  "        }\n"
  "        if (ifarg(3)) {\n"
  "            a3 = (uint32_t)*getarg(3);\n"
  "        }\n"
  "        *pv = nrnran123_newstream3((uint32_t)*getarg(1), a2, a3);\n"
  "        usingR123 = 1;\n"
  "    } else if( ifarg(1) ) {   // not a double, so assume hoc object type\n"
  "        void** pv = (void**)(&_p_rng_TM);\n"
  "        *pv = nrn_random_arg(1);\n"
  "    } else {  // no arg, so clear pointer\n"
  "        void** pv = (void**)(&_p_rng_TM);\n"
  "        *pv = (void*)0;\n"
  "    }\n"
  "    #endif\n"
  "    ENDVERBATIM\n"
  "}\n"
  "\n"
  "FUNCTION urand() {\n"
  "    VERBATIM\n"
  "    double value;\n"
  "    if ( usingR123 ) {\n"
  "        value = nrnran123_dblpick((nrnran123_State*)_p_rng_TM);\n"
  "    } else if (_p_rng_TM) {\n"
  "        #ifndef CORENEURON_BUILD\n"
  "        value = nrn_random_pick(_p_rng_TM);\n"
  "        #endif\n"
  "    } else {\n"
  "        value = 0.0;\n"
  "    }\n"
  "    _lurand = value;\n"
  "    ENDVERBATIM\n"
  "}\n"
  "\n"
  "FUNCTION brand(n, p) {\n"
  "    LOCAL result, count, success\n"
  "    success = 0\n"
  "    FROM count = 0 TO (n - 1) {\n"
  "        result = urand()\n"
  "        if(result <= p) {\n"
  "            success = success + 1\n"
  "        }\n"
  "    }\n"
  "    brand = success\n"
  "}\n"
  "\n"
  "FUNCTION bbsavestate() {\n"
  "    bbsavestate = 0\n"
  "    VERBATIM\n"
  "    #ifndef CORENEURON_BUILD\n"
  "        /* first arg is direction (0 save, 1 restore), second is array*/\n"
  "        /* if first arg is -1, fill xdir with the size of the array */\n"
  "        double *xdir, *xval, *hoc_pgetarg();\n"
  "        long nrn_get_random_sequence(void* r);\n"
  "        void nrn_set_random_sequence(void* r, int val);\n"
  "        xdir = hoc_pgetarg(1);\n"
  "        xval = hoc_pgetarg(2);\n"
  "        if (_p_rng_TM) {\n"
  "            // tell how many items need saving\n"
  "            if (*xdir == -1) {  // count items\n"
  "                if( usingR123 ) {\n"
  "                    *xdir = 2.0;\n"
  "                } else {\n"
  "                    *xdir = 1.0;\n"
  "                }\n"
  "                return 0.0;\n"
  "            } else if(*xdir ==0 ) {  // save\n"
  "                if( usingR123 ) {\n"
  "                    uint32_t seq;\n"
  "                    unsigned char which;\n"
  "                    nrnran123_getseq( (nrnran123_State*)_p_rng_TM, &seq, &which );\n"
  "                    xval[0] = (double) seq;\n"
  "                    xval[1] = (double) which;\n"
  "                } else {\n"
  "                    xval[0] = (double)nrn_get_random_sequence(_p_rng_TM);\n"
  "                }\n"
  "            } else {  // restore\n"
  "                if( usingR123 ) {\n"
  "                    nrnran123_setseq( (nrnran123_State*)_p_rng_TM, (uint32_t)xval[0], (char)xval[1] );\n"
  "                } else {\n"
  "                    nrn_set_random_sequence(_p_rng_TM, (long)(xval[0]));\n"
  "                }\n"
  "            }\n"
  "        }\n"
  "    #endif\n"
  "    ENDVERBATIM\n"
  "}\n"
  "\n"
  "VERBATIM\n"
  "static void bbcore_write(double* dArray, int* iArray, int* doffset, int* ioffset, _threadargsproto_) {\n"
  "    // make sure offset array non-null\n"
  "    if (iArray) {\n"
  "        // get handle to random123 instance\n"
  "        nrnran123_State** pv = (nrnran123_State**)(&_p_rng_TM);\n"
  "        // get location for storing ids\n"
  "        uint32_t* ia = ((uint32_t*)iArray) + *ioffset;\n"
  "        // retrieve/store identifier seeds\n"
  "        nrnran123_getids3(*pv, ia, ia+1, ia+2);\n"
  "        // retrieve/store stream sequence\n"
  "        unsigned char which;\n"
  "        nrnran123_getseq(*pv, ia+3, &which);\n"
  "        ia[4] = (int)which;\n"
  "    }\n"
  "\n"
  "    // increment integer offset (2 identifier), no double data\n"
  "    *ioffset += 5;\n"
  "    *doffset += 0;\n"
  "}\n"
  "\n"
  "static void bbcore_read(double* dArray, int* iArray, int* doffset, int* ioffset, _threadargsproto_) {\n"
  "    // make sure it's not previously set\n"
  "    assert(!_p_rng_TM);\n"
  "    uint32_t* ia = ((uint32_t*)iArray) + *ioffset;\n"
  "    // make sure non-zero identifier seeds\n"
  "    if (ia[0] != 0 || ia[1] != 0 || ia[2] != 0) {\n"
  "        nrnran123_State** pv = (nrnran123_State**)(&_p_rng_TM);\n"
  "        // get new stream\n"
  "        *pv = nrnran123_newstream3(ia[0], ia[1], ia[2]);\n"
  "        // restore sequence\n"
  "        nrnran123_setseq(*pv, ia[3], (char)ia[4]);\n"
  "    }\n"
  "    // increment intger offset (2 identifiers), no double data\n"
  "    *ioffset += 5;\n"
  "    *doffset += 0;\n"
  "}\n"
  "ENDVERBATIM\n"
  "\n"
  ;
#endif
