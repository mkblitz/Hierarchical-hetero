#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _CaDynamics_E2_reg();
extern void _Ca_HVA_reg();
extern void _Ca_LVAst_reg();
extern void _GluSynapse_reg();
extern void _GluSynapse_TM_MK_reg();
extern void _GluSynapse_andras_reg();
extern void _Ih_reg();
extern void _Im_reg();
extern void _K_Pst_reg();
extern void _K_Tst_reg();
extern void _NMDA_reg();
extern void _NaTa_t_reg();
extern void _NaTs2_t_reg();
extern void _Nap_Et2_reg();
extern void _ProbAMPA_reg();
extern void _ProbAMPANMDA_reg();
extern void _ProbAMPANMDA2_ratio_reg();
extern void _ProbAMPANMDA_EMS_reg();
extern void _ProbGABAA_reg();
extern void _ProbGABAAB_EMS_reg();
extern void _ProbGABAA_EMS_reg();
extern void _ProbNMDA_reg();
extern void _SK_E2_reg();
extern void _SKv3_1_reg();
extern void _double_exp_AMPA_NMDA_reg();
extern void _double_exp_AMPA_NMDA_Behabadi_reg();
extern void _epsp_reg();
extern void _kv_reg();
extern void _na_reg();

void modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," CaDynamics_E2.mod");
fprintf(stderr," Ca_HVA.mod");
fprintf(stderr," Ca_LVAst.mod");
fprintf(stderr," GluSynapse.mod");
fprintf(stderr," GluSynapse_TM_MK.mod");
fprintf(stderr," GluSynapse_andras.mod");
fprintf(stderr," Ih.mod");
fprintf(stderr," Im.mod");
fprintf(stderr," K_Pst.mod");
fprintf(stderr," K_Tst.mod");
fprintf(stderr," NMDA.mod");
fprintf(stderr," NaTa_t.mod");
fprintf(stderr," NaTs2_t.mod");
fprintf(stderr," Nap_Et2.mod");
fprintf(stderr," ProbAMPA.mod");
fprintf(stderr," ProbAMPANMDA.mod");
fprintf(stderr," ProbAMPANMDA2_ratio.mod");
fprintf(stderr," ProbAMPANMDA_EMS.mod");
fprintf(stderr," ProbGABAA.mod");
fprintf(stderr," ProbGABAAB_EMS.mod");
fprintf(stderr," ProbGABAA_EMS.mod");
fprintf(stderr," ProbNMDA.mod");
fprintf(stderr," SK_E2.mod");
fprintf(stderr," SKv3_1.mod");
fprintf(stderr," double_exp_AMPA_NMDA.mod");
fprintf(stderr," double_exp_AMPA_NMDA_Behabadi.mod");
fprintf(stderr," epsp.mod");
fprintf(stderr," kv.mod");
fprintf(stderr," na.mod");
fprintf(stderr, "\n");
    }
_CaDynamics_E2_reg();
_Ca_HVA_reg();
_Ca_LVAst_reg();
_GluSynapse_reg();
_GluSynapse_TM_MK_reg();
_GluSynapse_andras_reg();
_Ih_reg();
_Im_reg();
_K_Pst_reg();
_K_Tst_reg();
_NMDA_reg();
_NaTa_t_reg();
_NaTs2_t_reg();
_Nap_Et2_reg();
_ProbAMPA_reg();
_ProbAMPANMDA_reg();
_ProbAMPANMDA2_ratio_reg();
_ProbAMPANMDA_EMS_reg();
_ProbGABAA_reg();
_ProbGABAAB_EMS_reg();
_ProbGABAA_EMS_reg();
_ProbNMDA_reg();
_SK_E2_reg();
_SKv3_1_reg();
_double_exp_AMPA_NMDA_reg();
_double_exp_AMPA_NMDA_Behabadi_reg();
_epsp_reg();
_kv_reg();
_na_reg();
}
