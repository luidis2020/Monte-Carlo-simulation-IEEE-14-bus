
case_IEEE14_start;
ieee14_load=barras_ipfc(:,8:9);
activate_PEM=3;
activate_normalPF=2;
load_flow_analysis_original;
V_normal=V;
gen_normal=gen;
Pcal_det=Pcal;
Qcal_det=Qcal;
T_det=T;
V_det=V;
J_det=J;
Y_det=Y;