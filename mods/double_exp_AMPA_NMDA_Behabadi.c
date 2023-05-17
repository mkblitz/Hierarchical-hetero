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
 
#define nrn_init _nrn_init__double_exp_AMPA_NMDA_Behabadi
#define _nrn_initial _nrn_initial__double_exp_AMPA_NMDA_Behabadi
#define nrn_cur _nrn_cur__double_exp_AMPA_NMDA_Behabadi
#define _nrn_current _nrn_current__double_exp_AMPA_NMDA_Behabadi
#define nrn_jacob _nrn_jacob__double_exp_AMPA_NMDA_Behabadi
#define nrn_state _nrn_state__double_exp_AMPA_NMDA_Behabadi
#define _net_receive _net_receive__double_exp_AMPA_NMDA_Behabadi 
#define _f_mgblock _f_mgblock__double_exp_AMPA_NMDA_Behabadi 
#define betadyn betadyn__double_exp_AMPA_NMDA_Behabadi 
#define mgblock mgblock__double_exp_AMPA_NMDA_Behabadi 
 
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
#define tau1 _p[0]
#define tau2 _p[1]
#define tau1_NMDA _p[2]
#define tau2_NMDA _p[3]
#define e _p[4]
#define mg _p[5]
#define NMDA_ratio _p[6]
#define i _p[7]
#define i_NMDA _p[8]
#define i_AMPA _p[9]
#define g_NMDA _p[10]
#define g_AMPA _p[11]
#define vsyn _p[12]
#define A_AMPA _p[13]
#define B_AMPA _p[14]
#define A_NMDA _p[15]
#define B_NMDA _p[16]
#define f_NMDA _p[17]
#define factor_AMPA _p[18]
#define factor_NMDA _p[19]
#define DA_AMPA _p[20]
#define DB_AMPA _p[21]
#define DA_NMDA _p[22]
#define DB_NMDA _p[23]
#define v _p[24]
#define _g _p[25]
#define _tsav _p[26]
#define _nd_area  *_ppvar[0]._pval
 
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
 static int hoc_nrnpointerindex =  -1;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 /* declaration of user functions */
 static double _hoc_mgblock();
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
 "mgblock", _hoc_mgblock,
 0, 0
};
 
static void _check_mgblock(double*, Datum*, Datum*, _NrnThread*); 
static void _check_table_thread(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, int _type) {
   _check_mgblock(_p, _ppvar, _thread, _nt);
 }
 /* declare global and static user variables */
#define usetable usetable_double_exp_AMPA_NMDA_Behabadi
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_double_exp_AMPA_NMDA_Behabadi", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "tau1", "ms",
 "tau2", "ms",
 "tau1_NMDA", "ms",
 "tau2_NMDA", "ms",
 "e", "mV",
 "mg", "mM",
 "A_AMPA", "uS",
 "B_AMPA", "uS",
 "A_NMDA", "uS",
 "B_NMDA", "uS",
 "i", "nA",
 "i_NMDA", "nA",
 "i_AMPA", "nA",
 "g_NMDA", "uS",
 "g_AMPA", "uS",
 "vsyn", "mV",
 0,0
};
 static double A_NMDA0 = 0;
 static double A_AMPA0 = 0;
 static double B_NMDA0 = 0;
 static double B_AMPA0 = 0;
 static double delta_t = 0.01;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "usetable_double_exp_AMPA_NMDA_Behabadi", &usetable_double_exp_AMPA_NMDA_Behabadi,
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
 static void _hoc_destroy_pnt(_vptr) void* _vptr; {
   destroy_point_process(_vptr);
}
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[2]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"double_exp_AMPA_NMDA_Behabadi",
 "tau1",
 "tau2",
 "tau1_NMDA",
 "tau2_NMDA",
 "e",
 "mg",
 "NMDA_ratio",
 0,
 "i",
 "i_NMDA",
 "i_AMPA",
 "g_NMDA",
 "g_AMPA",
 "vsyn",
 0,
 "A_AMPA",
 "B_AMPA",
 "A_NMDA",
 "B_NMDA",
 0,
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
 	_p = nrn_prop_data_alloc(_mechtype, 27, _prop);
 	/*initialize range parameters*/
 	tau1 = 0.05;
 	tau2 = 0.5;
 	tau1_NMDA = 2.1;
 	tau2_NMDA = 18.8;
 	e = 0;
 	mg = 1;
 	NMDA_ratio = 0;
  }
 	_prop->param = _p;
 	_prop->param_size = 27;
  if (!nrn_point_prop_) {
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
  }
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _net_receive(Point_process*, double*, double);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _double_exp_AMPA_NMDA_Behabadi_reg() {
	int _vectorized = 1;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex, 1,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_table_reg(_mechtype, _check_table_thread);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 27, 3);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_size[_mechtype] = 1;
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 double_exp_AMPA_NMDA_Behabadi C:/Users/mkbli/Dropbox/deep_dendrite/mods/double_exp_AMPA_NMDA_Behabadi.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double *_t_f_NMDA;
static int _reset;
static char *modelname = "double Exp NMDA synapse";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int _f_mgblock(_threadargsprotocomma_ double);
static int mgblock(_threadargsprotocomma_ double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static void _n_mgblock(_threadargsprotocomma_ double _lv);
 static int _slist1[4], _dlist1[4];
 static int betadyn(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   DA_AMPA = - A_AMPA / tau1 ;
   DB_AMPA = - B_AMPA / tau2 ;
   DA_NMDA = - A_NMDA / tau1_NMDA ;
   DB_NMDA = - B_NMDA / tau2_NMDA ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 DA_AMPA = DA_AMPA  / (1. - dt*( ( - 1.0 ) / tau1 )) ;
 DB_AMPA = DB_AMPA  / (1. - dt*( ( - 1.0 ) / tau2 )) ;
 DA_NMDA = DA_NMDA  / (1. - dt*( ( - 1.0 ) / tau1_NMDA )) ;
 DB_NMDA = DB_NMDA  / (1. - dt*( ( - 1.0 ) / tau2_NMDA )) ;
  return 0;
}
 /*END CVODE*/
 static int betadyn (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) { {
    A_AMPA = A_AMPA + (1. - exp(dt*(( - 1.0 ) / tau1)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau1 ) - A_AMPA) ;
    B_AMPA = B_AMPA + (1. - exp(dt*(( - 1.0 ) / tau2)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau2 ) - B_AMPA) ;
    A_NMDA = A_NMDA + (1. - exp(dt*(( - 1.0 ) / tau1_NMDA)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau1_NMDA ) - A_NMDA) ;
    B_NMDA = B_NMDA + (1. - exp(dt*(( - 1.0 ) / tau2_NMDA)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau2_NMDA ) - B_NMDA) ;
   }
  return 0;
}
 
static void _net_receive (_pnt, _args, _lflag) Point_process* _pnt; double* _args; double _lflag; 
{  double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   _thread = (Datum*)0; _nt = (_NrnThread*)_pnt->_vnt;   _p = _pnt->_prop->param; _ppvar = _pnt->_prop->dparam;
  if (_tsav > t){ extern char* hoc_object_name(); hoc_execerror(hoc_object_name(_pnt->ob), ":Event arrived out of order. Must call ParallelContext.set_maxstep AFTER assigning minimum NetCon.delay");}
 _tsav = t; {
   _args[0] = _args[0] * .001 ;
     if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = A_AMPA;
    double __primary = (A_AMPA - _args[0] * factor_AMPA) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - 1.0 ) / tau1 ) ) )*( - ( 0.0 ) / ( ( - 1.0 ) / tau1 ) - __primary );
    A_AMPA += __primary;
  } else {
 A_AMPA = A_AMPA - _args[0] * factor_AMPA ;
     }
   if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = B_AMPA;
    double __primary = (B_AMPA + _args[0] * factor_AMPA) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - 1.0 ) / tau2 ) ) )*( - ( 0.0 ) / ( ( - 1.0 ) / tau2 ) - __primary );
    B_AMPA += __primary;
  } else {
 B_AMPA = B_AMPA + _args[0] * factor_AMPA ;
     }
   if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = A_NMDA;
    double __primary = (A_NMDA - _args[0] * NMDA_ratio * factor_NMDA) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - 1.0 ) / tau1_NMDA ) ) )*( - ( 0.0 ) / ( ( - 1.0 ) / tau1_NMDA ) - __primary );
    A_NMDA += __primary;
  } else {
 A_NMDA = A_NMDA - _args[0] * NMDA_ratio * factor_NMDA ;
     }
   if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = B_NMDA;
    double __primary = (B_NMDA + _args[0] * NMDA_ratio * factor_NMDA) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - 1.0 ) / tau2_NMDA ) ) )*( - ( 0.0 ) / ( ( - 1.0 ) / tau2_NMDA ) - __primary );
    B_NMDA += __primary;
  } else {
 B_NMDA = B_NMDA + _args[0] * NMDA_ratio * factor_NMDA ;
     }
 } }
 static double _mfac_mgblock, _tmin_mgblock;
  static void _check_mgblock(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  static double _sav_mg;
  if (!usetable) {return;}
  if (_sav_mg != mg) { _maktable = 1;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_mgblock =  - 100.0 ;
   _tmax =  100.0 ;
   _dx = (_tmax - _tmin_mgblock)/201.; _mfac_mgblock = 1./_dx;
   for (_i=0, _x=_tmin_mgblock; _i < 202; _x += _dx, _i++) {
    _f_mgblock(_p, _ppvar, _thread, _nt, _x);
    _t_f_NMDA[_i] = f_NMDA;
   }
   _sav_mg = mg;
  }
 }

 static int mgblock(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _lv) { 
#if 0
_check_mgblock(_p, _ppvar, _thread, _nt);
#endif
 _n_mgblock(_p, _ppvar, _thread, _nt, _lv);
 return 0;
 }

 static void _n_mgblock(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_mgblock(_p, _ppvar, _thread, _nt, _lv); return; 
}
 _xi = _mfac_mgblock * (_lv - _tmin_mgblock);
 if (isnan(_xi)) {
  f_NMDA = _xi;
  return;
 }
 if (_xi <= 0.) {
 f_NMDA = _t_f_NMDA[0];
 return; }
 if (_xi >= 201.) {
 f_NMDA = _t_f_NMDA[201];
 return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 f_NMDA = _t_f_NMDA[_i] + _theta*(_t_f_NMDA[_i+1] - _t_f_NMDA[_i]);
 }

 
static int  _f_mgblock ( _threadargsprotocomma_ double _lv ) {
   f_NMDA = 1.0 / ( 1.0 + 0.3 * exp ( - 0.1 * ( _lv ) ) ) ;
    return 0; }
 
static double _hoc_mgblock(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (_NrnThread*)((Point_process*)_vptr)->_vnt;
 
#if 1
 _check_mgblock(_p, _ppvar, _thread, _nt);
#endif
 _r = 1.;
 mgblock ( _p, _ppvar, _thread, _nt, *getarg(1) );
 return(_r);
}
 
static int _ode_count(int _type){ return 4;}
 
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
	for (_i=0; _i < 4; ++_i) {
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

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  A_NMDA = A_NMDA0;
  A_AMPA = A_AMPA0;
  B_NMDA = B_NMDA0;
  B_AMPA = B_AMPA0;
 {
   double _ltp_AMPA , _ltp_NMDA ;
 A_AMPA = 0. ;
   B_AMPA = 0. ;
   A_NMDA = 0. ;
   B_NMDA = 0. ;
   _ltp_AMPA = ( tau1 * tau2 ) / ( tau2 - tau1 ) * log ( tau2 / tau1 ) ;
   factor_AMPA = - exp ( - _ltp_AMPA / tau1 ) + exp ( - _ltp_AMPA / tau2 ) ;
   factor_AMPA = 1. / factor_AMPA ;
   _ltp_NMDA = ( tau1_NMDA * tau2_NMDA ) / ( tau2_NMDA - tau1_NMDA ) * log ( tau2_NMDA / tau1_NMDA ) ;
   factor_NMDA = - exp ( - _ltp_NMDA / tau1_NMDA ) + exp ( - _ltp_NMDA / tau2_NMDA ) ;
   factor_NMDA = 1. / factor_NMDA ;
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

#if 0
 _check_mgblock(_p, _ppvar, _thread, _nt);
#endif
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
   mgblock ( _threadargscomma_ v ) ;
   g_AMPA = ( A_AMPA + B_AMPA ) ;
   g_NMDA = ( A_NMDA + B_NMDA ) ;
   i_NMDA = f_NMDA * ( g_NMDA ) * ( v - e ) ;
   i_AMPA = ( g_AMPA ) * ( v - e ) ;
   i = i_AMPA + i_NMDA ;
   vsyn = v ;
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
 {   betadyn(_p, _ppvar, _thread, _nt);
  }}}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(A_AMPA) - _p;  _dlist1[0] = &(DA_AMPA) - _p;
 _slist1[1] = &(B_AMPA) - _p;  _dlist1[1] = &(DB_AMPA) - _p;
 _slist1[2] = &(A_NMDA) - _p;  _dlist1[2] = &(DA_NMDA) - _p;
 _slist1[3] = &(B_NMDA) - _p;  _dlist1[3] = &(DB_NMDA) - _p;
   _t_f_NMDA = makevector(202*sizeof(double));
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "double_exp_AMPA_NMDA_Behabadi.mod";
static const char* nmodl_file_text = 
  "TITLE   double Exp NMDA synapse\n"
  "\n"
  "NEURON {\n"
  "	THREADSAFE\n"
  "		POINT_PROCESS double_exp_AMPA_NMDA_Behabadi\n"
  "		RANGE tau1, tau2, tau1_NMDA, tau2_NMDA, mg, i, e, NMDA_ratio, i_NMDA, i_AMPA, g_NMDA, g_AMPA, vsyn\n"
  "		NONSPECIFIC_CURRENT i\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(nA) = (nanoamp)\n"
  "	(mV) = (millivolt)\n"
  "	(uS) = (microsiemens)\n"
  "	(mM) = (milli/liter)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	tau1 = .05 (ms)\n"
  "	tau2 = .5 (ms)\n"
  "	tau1_NMDA = 2.1 (ms)\n"
  "    tau2_NMDA = 18.8 (ms)\n"
  "	e  = 0    (mV)   : reversal potential, Dalby 2003\n"
  "	mg = 1      (mM)    : external magnesium concentration\n"
  "    NMDA_ratio = 0\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	v (mV)   		: postsynaptic voltage\n"
  "	i (nA)   		: nonspecific current = g*(v - Erev) \n"
  "	i_NMDA (nA)\n"
  "	i_AMPA (nA)\n"
  "	g_NMDA (uS)\n"
  "	g_AMPA (uS)\n"
  "	vsyn   (mV)\n"
  "\n"
  "	f_NMDA				: voltage dependendent magnesium blockade\n"
  "\n"
  "	factor_AMPA\n"
  "	factor_NMDA\n"
  "}\n"
  "\n"
  "\n"
  "STATE { \n"
  "	A_AMPA (uS)\n"
  "	B_AMPA (uS)\n"
  "    A_NMDA (uS)\n"
  "    B_NMDA (uS)\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "	LOCAL tp_AMPA, tp_NMDA\n"
  "\n"
  "    A_AMPA = 0.\n"
  "	B_AMPA = 0.\n"
  "    A_NMDA = 0.\n"
  "    B_NMDA = 0.\n"
  "\n"
  "    tp_AMPA = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)\n"
  "    factor_AMPA = -exp(-tp_AMPA/tau1) + exp(-tp_AMPA/tau2)\n"
  "    factor_AMPA = 1./factor_AMPA\n"
  "\n"
  "    tp_NMDA = (tau1_NMDA*tau2_NMDA)/(tau2_NMDA - tau1_NMDA) * log(tau2_NMDA/tau1_NMDA)\n"
  "    factor_NMDA = -exp(-tp_NMDA/tau1_NMDA) + exp(-tp_NMDA/tau2_NMDA)\n"
  "    factor_NMDA = 1./factor_NMDA\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE betadyn METHOD cnexp\n"
  "  	mgblock(v)\n"
  "	:i = ((A_AMPA + B_AMPA) + f_NMDA * (A_NMDA + B_NMDA)) * (v - e)\n"
  "	g_AMPA = (A_AMPA + B_AMPA)\n"
  "	g_NMDA = (A_NMDA + B_NMDA)\n"
  "	i_NMDA = f_NMDA * (g_NMDA) * (v - e)\n"
  "	i_AMPA = (g_AMPA) * (v - e)\n"
  "	i = i_AMPA + i_NMDA\n"
  "	vsyn = v\n"
  "}\n"
  "\n"
  "DERIVATIVE betadyn {\n"
  "    A_AMPA' = -A_AMPA/tau1\n"
  "	B_AMPA' = -B_AMPA/tau2\n"
  "    A_NMDA' = -A_NMDA/tau1_NMDA\n"
  "    B_NMDA' = -B_NMDA/tau2_NMDA\n"
  "}\n"
  "\n"
  "NET_RECEIVE(weight (uS)) {\n"
  "	weight = weight *.001 \n"
  "	A_AMPA = A_AMPA - weight*factor_AMPA\n"
  "	B_AMPA = B_AMPA + weight*factor_AMPA\n"
  "    A_NMDA = A_NMDA - weight*NMDA_ratio*factor_NMDA\n"
  "    B_NMDA = B_NMDA + weight*NMDA_ratio*factor_NMDA\n"
  "}\n"
  "\n"
  "\n"
  "PROCEDURE mgblock( v(mV) ) {\n"
  "	: from Jahr & Stevens\n"
  "\n"
  "	TABLE f_NMDA DEPEND mg\n"
  "		FROM -100 TO 100 WITH 201\n"
  "\n"
  "	:f_NMDA = 1 / (1 + exp(0.062 (/mV) * -v) * (mg / 3.57 (mM)))\n"
  "	f_NMDA = 1/(1 + 0.3*exp(-0.1 (/mV) *(v)))\n"
  "}\n"
  "\n"
  "\n"
  "\n"
  "COMMENT\n"
  "Author Johan Hake (c) spring 2004\n"
  ":     Summate input from many presynaptic sources and saturate \n"
  ":     each one of them during heavy presynaptic firing\n"
  "\n"
  ": [1] Destexhe, A., Z. F. Mainen and T. J. Sejnowski (1998)\n"
  ":     Kinetic models of synaptic transmission\n"
  ":     In C. Koch and I. Segev (Eds.), Methods in Neuronal Modeling\n"
  "\n"
  ": [2] Rotter, S. and M. Diesmann (1999) Biol. Cybern. 81, 381-402\n"
  ":     Exact digital simulation of time-invariant linear systems with application \n"
  ":     to neural modeling\n"
  "\n"
  "Mainen ZF, Malinow R, Svoboda K (1999) Nature. 399, 151-155.\n"
  "Synaptic calcium transients in single spines indicate that NMDA\n"
  "receptors are not saturated.\n"
  "\n"
  "Chapman DE, Keefe KA, Wilcox KS (2003) J Neurophys. 89: 69-80.\n"
  "Evidence for functionally distinct synaptic nmda receptors in ventromedial\n"
  "vs. dorsolateral striatum.\n"
  "\n"
  "Dalby, N. O., and Mody, I. (2003). Activation of NMDA receptors in rat\n"
  "dentate gyrus granule cells by spontaneous and evoked transmitter\n"
  "release. J Neurophysiol 90, 786-797.\n"
  "\n"
  "Jahr CE, Stevens CF. (1990) Voltage dependence of NMDA activated\n"
  "macroscopic conductances predicted by single channel kinetics. J\n"
  "Neurosci 10: 3178, 1990.\n"
  "\n"
  "Gutfreund H, Kinetics for the Life Sciences, Cambridge University Press, 1995, pg 234.\n"
  "(suggested by Ted Carnevale)\n"
  "ENDCOMMENT\n"
  ;
#endif
