function [T,Y,yinit,param, allNames, allValues] = M1R_Cascade_with_Hodgkin_Huxley_Experiment(argTimeSpan,argYinit,argParam)
% [T,Y,yinit,param] = M1R_Cascade_with_Hodgkin_Huxley_Experiment(argTimeSpan,argYinit,argParam)
%
% input:
%     argTimeSpan is a vector of start and stop times (e.g. timeSpan = [0 10.0])
%     argYinit is a vector of initial conditions for the state variables (optional)
%     argParam is a vector of values for the parameters (optional)
%
% output:
%     T is the vector of times
%     Y is the vector of state variables
%     yinit is the initial conditions that were used
%     param is the parameter vector that was used
%     allNames is the output solution variable names
%     allValues is the output solution variable values corresponding to the names
%
%     example of running this file: [T,Y,yinit,param,allNames,allValues] = myMatlabFunc; <-(your main function name)
%

%
% Default time span
%
timeSpan = [0.0 1.0];

% output variable lengh and names
numVars = 153;
allNames = {'PIP2_M';'h_open';'Ga_GDP_M';'RLGbeta_M';'PH_CFP_C';'n_open';'RLG_GDP_M';'R_M';'PI_M';'PH_CFP_PIP2_M';'m_open';'Gbeta_M';'RGbeta_M';'GaGTP_PLC_M';'DAG_M';'IP3_PH_YFP_C';'GaGTP_M';'IP3_C';'Voltage_M';'RL_M';'PI4P_M';'KCNQ_M';'alpha_n';'alpha_m';'J_PI4K_4Pase';'alpha_h';'UnitFactor_uM_um3_molecules_neg_1';'K_RG_GDP_M_total';'RG_GDP_M';'J_NE_RG';'unitFactor_Na_Flux';'K_G_GDP_M_total';'unitFactor_K_flux';'Kr_PH_YFP_PIP2';'J_IP3Pase';'Kf_G2';'KfG1';'Kf_G1';'K_PLC_M_total';'PLC_M';'J_PLCassoc';'KG2';'KG1';'K_C';'K_m_closed_total';'Na_EX';'K_oxoM_EX_total';'oxoM_EX';'KfL1';'KfL2';'Kf_L2beta';'KL2';'KrL2';'Kr_L2beta';'J_L2beta';'Kr_G2beta';'Kf_PH_YFP_IP3';'Kr_PH_YFP_IP3';'K_PH_YFP_C_total';'PH_YFP_C';'J_PH_YFP_IP3';'K_PH_YFP_PIP2_M_total';'KrG1';'K_EX';'Kv_KCNQ2_3';'K_KCNQ_PIP2_M_total';'KCNQ_PIP2_M';'Kp_KCNQ2_3';'PPo';'C_KCNQ2_3';'I_KCNQ2_3';'unitFactor_KCNQ2_3';'J_KCNQ2_3';'K_GaGDP_PLC_M_total';'PH_YFP_PIP2_M';'Kr_G1beta';'UnitFactor_mV_pF_s_neg_1_pA_neg_1';'G_GDP_M';'J_NE_G';'KFlux_M_C';'unitFactor_Leak_flux';'Kf_PIP2bindKCNQ';'KD_PIP2_KCNQ';'Kr_PIP2bindKCNQ';'J_PIP2bindKCNQ';'K_IP3_PH_CFP_C_total';'IP3_PH_CFP_C';'Na_C';'C_Na_Flux';'I_Na_Flux';'Kr_PH_CFP_PIP2';'Kf_PH_YFP_PIP2';'GaGDP_PLC_M';'allGaPLC';'beta_n';'beta_m';'device_totalCurrClampElectrode.F';'device_totalCurrClampElectrode.I';'beta_h';'Kf_PH_CFP_IP3';'Kr_PH_CFP_IP3';'J_PH_CFP_IP3';'K_h_closed_total';'h_closed';'I_Leak_flux';'C_K_flux';'I_K_flux';'F_M';'J_GTPase_GaP';'Kf_G2beta';'KFlux_M_EX';'J_PH_YFP_PIP2';'m_closed';'J_Na_Gate_Act';'J_NE_GaP';'J_PLCdiss';'Kf_G1beta';'J_Leak_flux';'J_NE_RLG';'J_PIP2hydr';'I_M';'Kr_L2';'Kf_L2';'J_L2';'J_PLC_on_PI4P';'Kf_L1';'Kr_L1';'J_L1';'Kf_PH_CFP_PIP2';'J_PH_CFP_PIP2';'J_reconstitution';'J_Na_Flux';'device_M.Capacitance';'allRL';'J_Na_Gate_Inact';'Kr_G2';'J_G2';'Kr_G1';'J_G1';'K_n_closed_total';'Ion_EX';'allPH_cytoplasm';'J_G2beta';'Ion_C';'n_closed';'J_K_Gate';'allPIP2M';'J_G1beta';'allGabg';'J_DAGPase';'J_GTPase_Ga';'J_K_flux';'allRGbeta';};

if nargin >= 1
	if length(argTimeSpan) > 0
		%
		% TimeSpan overridden by function arguments
		%
		timeSpan = argTimeSpan;
	end
end
%
% Default Initial Conditions
%
yinit = [
	5000.0;		% yinit(1) is the initial condition for 'PIP2_M'
	0.57;		% yinit(2) is the initial condition for 'h_open'
	0.0;		% yinit(3) is the initial condition for 'Ga_GDP_M'
	0.0;		% yinit(4) is the initial condition for 'RLGbeta_M'
	0.0;		% yinit(5) is the initial condition for 'PH_CFP_C'
	0.33;		% yinit(6) is the initial condition for 'n_open'
	0.0;		% yinit(7) is the initial condition for 'RLG_GDP_M'
	1.0;		% yinit(8) is the initial condition for 'R_M'
	140000.0;		% yinit(9) is the initial condition for 'PI_M'
	0.0;		% yinit(10) is the initial condition for 'PH_CFP_PIP2_M'
	0.058;		% yinit(11) is the initial condition for 'm_open'
	0.0;		% yinit(12) is the initial condition for 'Gbeta_M'
	0.0;		% yinit(13) is the initial condition for 'RGbeta_M'
	0.0;		% yinit(14) is the initial condition for 'GaGTP_PLC_M'
	2000.0;		% yinit(15) is the initial condition for 'DAG_M'
	0.0;		% yinit(16) is the initial condition for 'IP3_PH_YFP_C'
	0.0;		% yinit(17) is the initial condition for 'GaGTP_M'
	0.16;		% yinit(18) is the initial condition for 'IP3_C'
	-65.0;		% yinit(19) is the initial condition for 'Voltage_M'
	0.0;		% yinit(20) is the initial condition for 'RL_M'
	4000.0;		% yinit(21) is the initial condition for 'PI4P_M'
	4.0;		% yinit(22) is the initial condition for 'KCNQ_M'
];
if nargin >= 2
	if length(argYinit) > 0
		%
		% initial conditions overridden by function arguments
		%
		yinit = argYinit;
	end
end
%
% Default Parameters
%   constants are only those "Constants" from the Math Description that are just floating point numbers (no identifiers)
%   note: constants of the form "A_init" are really initial conditions and are treated in "yinit"
%
param = [
	1.0;		% param(1) is 'netValence_Na_Gate_Inact'
	0.16;		% param(2) is 'IP3_C_init_uM'
	96480.0;		% param(3) is 'mlabfix_F_'
	1.0;		% param(4) is 'netValence_PIP2bindKCNQ'
	4.0;		% param(5) is 'KCNQ_M_init_molecules_um_2'
	0.0;		% param(6) is 'Kr_PLCassoc'
	0.67;		% param(7) is 'n_closed_init_molecules_um_2'
	1.0;		% param(8) is 'netValence_G2'
	1.0;		% param(9) is 'netValence_G1'
	0.65;		% param(10) is 'Kf_NE_RLG'
	0.006;		% param(11) is 'Kr_PI4K_4Pase'
	-51.4;		% param(12) is 'V_leak'
	1.0;		% param(13) is 'netValence_GTPase_Ga'
	1.0;		% param(14) is 'netValence_PH_YFP_PIP2'
	0.43;		% param(15) is 'h_closed_init_molecules_um_2'
	0.1;		% param(16) is 'K_plc'
	0.0;		% param(17) is 'IP3_PH_YFP_C_init_uM'
	1.0;		% param(18) is 'netValence_Na_Gate_Act'
	1.0;		% param(19) is 'C_M'
	0.0;		% param(20) is 'RLG_GDP_M_init_molecules_um_2'
	1.0;		% param(21) is 'carrierValence_K_flux'
	0.0;		% param(22) is 'RG_GDP_M_init_molecules_um_2'
	2000.0;		% param(23) is 'DAG_M_init_molecules_um_2'
	0.026;		% param(24) is 'Kf_GTPase_Ga'
	0.68;		% param(25) is 'KrG2'
	0.0;		% param(26) is 'GaGTP_M_init_molecules_um_2'
	1.0;		% param(27) is 'Kf_reconstitution'
	40.0;		% param(28) is 'G_GDP_M_init_molecules_um_2'
	0.0;		% param(29) is 'theta'
	0.71;		% param(30) is 'Kf_PLCdiss'
	6563.0;		% param(31) is 'Size_M'
	1.0;		% param(32) is 'netValence_L2beta'
	50000.0;		% param(33) is 'Size_C'
	5000.0;		% param(34) is 'PIP2_M_init_molecules_um_2'
	2.6E-4;		% param(35) is 'Kf_PI4K_4Pase'
	1.866;		% param(36) is 'z_Kv'
	0.0;		% param(37) is 'PH_YFP_PIP2_M_init_molecules_um_2'
	0.0026666;		% param(38) is 'KfG2'
	45.12;		% param(39) is 'Kv0'
	1.0;		% param(40) is 'carrierValence_Leak_flux'
	0.0;		% param(41) is 'RGbeta_M_init_molecules_um_2'
	1.0;		% param(42) is 'carrierValence_KCNQ2_3'
	1.0;		% param(43) is 'netValence_PI4K_4Pase'
	1.5E-5;		% param(44) is 'Kf_NE_G'
	1.0;		% param(45) is 'netValence_PLCdiss'
	4000.0;		% param(46) is 'PI4P_M_init_molecules_um_2'
	0.0025;		% param(47) is 'PLC_basal'
	1.5E-5;		% param(48) is 'Kf_NE_RG'
	0.035;		% param(49) is 'gKCNQ_max'
	-65.0;		% param(50) is 'V_K_closed'
	100.0;		% param(51) is 'Ion_EX_init_uM'
	-35.0;		% param(52) is 'V_Na_Inact_close'
	1.0;		% param(53) is 'netValence_GTPase_GaP'
	1000.0;		% param(54) is 'K_millivolts_per_volt'
	1.0;		% param(55) is 'netValence_PLC_on_PI4P'
	1.0E-9;		% param(56) is 'mlabfix_K_GHK_'
	0.02;		% param(57) is 'Kg_KCNQ2_3'
	1.0;		% param(58) is 'netValence_PH_CFP_PIP2'
	0.05;		% param(59) is 'speed_PIP2_KCNQ'
	0.75;		% param(60) is 'kp0'
	9.648E-5;		% param(61) is 'mlabfix_F_nmol_'
	0.33;		% param(62) is 'n_open_init_molecules_um_2'
	2.0;		% param(63) is 'KD_PH_PIP2'
	0.058;		% param(64) is 'm_open_init_molecules_um_2'
	0.942;		% param(65) is 'm_closed_init_molecules_um_2'
	0.0;		% param(66) is 'IP3_PH_CFP_C_init_uM'
	0.08;		% param(67) is 'K_IP3ase'
	1.0;		% param(68) is 'netValence_PIP2hydr'
	0.14;		% param(69) is 'PLC_efficiency_PIP'
	1.0;		% param(70) is 'R_M_init_molecules_um_2'
	300.0;		% param(71) is 'mlabfix_T_'
	1.0;		% param(72) is 'netValence_NE_GaP'
	140000.0;		% param(73) is 'PI_M_init_molecules_um_2'
	0.0;		% param(74) is 'KCNQ_PIP2_M_init_molecules_um_2'
	0.0;		% param(75) is 'Ga_GDP_M_init_molecules_um_2'
	0.1;		% param(76) is 'KD_PH_IP3'
	0.0;		% param(77) is 'Kr_NE_G'
	0.0;		% param(78) is 'PH_YFP_C_init_uM'
	8314.0;		% param(79) is 'mlabfix_R_'
	-65.0;		% param(80) is 'V_Na_Act_close'
	1.0;		% param(81) is 'Kf_PLCassoc'
	1.0;		% param(82) is 'netValence_NE_RLG'
	0.57;		% param(83) is 'h_open_init_molecules_um_2'
	1.0;		% param(84) is 'Hill_binding'
	0.0;		% param(85) is 'GaGTP_PLC_M_init_molecules_um_2'
	397000.0;		% param(86) is 'K_C_init_uM'
	1.0;		% param(87) is 'netValence_NE_RG'
	0.0;		% param(88) is 'RL_M_init_molecules_um_2'
	0.0;		% param(89) is 'GaGDP_PLC_M_init_molecules_um_2'
	0.02;		% param(90) is 'Kf_DAGPase'
	50000.0;		% param(91) is 'Na_C_init_uM'
	50000.0;		% param(92) is 'Size_EX'
	0.0;		% param(93) is 'PH_CFP_PIP2_M_init_molecules_um_2'
	20000.0;		% param(94) is 'K_EX_init_uM'
	-3.0E-4;		% param(95) is 'gL'
	1.0;		% param(96) is 'netValence_DAGPase'
	-65.0;		% param(97) is 'Voltage_M_init'
	100.0;		% param(98) is 'Ion_C_init_uM'
	1.0;		% param(99) is 'carrierValence_Na_Flux'
	-55.0;		% param(100) is 'V_K_open'
	1.0;		% param(101) is 'speed_PH_PIP2'
	3.141592653589793;		% param(102) is 'mlabfix_PI_'
	-65.0;		% param(103) is 'V_Na_Inact_open'
	8.0;		% param(104) is 'theta_KCNQ2_3'
	0.0;		% param(105) is 'PH_CFP_C_init_uM'
	0.0;		% param(106) is 'Kv'
	-40.0;		% param(107) is 'V_Na_Act_open'
	100.0;		% param(108) is 'alpha'
	0.0;		% param(109) is 'Kp'
	0.0;		% param(110) is 'Kr_PLCdiss'
	0.0;		% param(111) is 'Kg'
	0.0;		% param(112) is 'oxoM_EX_init_uM'
	1.0;		% param(113) is 'netValence_L2'
	1.0;		% param(114) is 'netValence_L1'
	2000.0;		% param(115) is 'KA_PIP2_KCNQ'
	6.02E11;		% param(116) is 'mlabfix_N_pmol_'
	1.0;		% param(117) is 'netValence_reconstitution'
	0.0;		% param(118) is 'RLGbeta_M_init_molecules_um_2'
	1.0;		% param(119) is 'netValence_G2beta'
	0.0;		% param(120) is 'Kr_NE_RG'
	2.0;		% param(121) is 'KL1'
	437000.0;		% param(122) is 'Na_EX_init_uM'
	0.12;		% param(123) is 'gNa_max'
	1.0;		% param(124) is 'netValence_K_Gate'
	(1.0 ./ 602.0);		% param(125) is 'KMOLE'
	1.0;		% param(126) is 'netValence_NE_G'
	0.0;		% param(127) is 'Gbeta_M_init_molecules_um_2'
	0.036;		% param(128) is 'gK_max'
	1.0;		% param(129) is 'netValence_G1beta'
	1.0;		% param(130) is 'netValence_PLCassoc'
	5.555555555555555;		% param(131) is 'KrL1'
	10.0;		% param(132) is 'PLC_M_init_molecules_um_2'
	10.0;		% param(133) is 'speed_PH_IP3'
];
if nargin >= 3
	if length(argParam) > 0
		%
		% parameter values overridden by function arguments
		%
		param = argParam;
	end
end
%
% invoke the integrator
%
[T,Y] = ode15s(@f,timeSpan,yinit,odeset('OutputFcn',@odeplot),param,yinit);

% get the solution
all = zeros(size(T), numVars);
for i = 1:size(T)
	all(i,:) = getRow(T(i), Y(i,:), yinit, param);
end

allValues = all;
end

% -------------------------------------------------------
% get row data
function rowValue = getRow(t,y,y0,p)
	% State Variables
	PIP2_M = y(1);
	h_open = y(2);
	Ga_GDP_M = y(3);
	RLGbeta_M = y(4);
	PH_CFP_C = y(5);
	n_open = y(6);
	RLG_GDP_M = y(7);
	R_M = y(8);
	PI_M = y(9);
	PH_CFP_PIP2_M = y(10);
	m_open = y(11);
	Gbeta_M = y(12);
	RGbeta_M = y(13);
	GaGTP_PLC_M = y(14);
	DAG_M = y(15);
	IP3_PH_YFP_C = y(16);
	GaGTP_M = y(17);
	IP3_C = y(18);
	Voltage_M = y(19);
	RL_M = y(20);
	PI4P_M = y(21);
	KCNQ_M = y(22);
	% Constants
	netValence_Na_Gate_Inact = p(1);
	IP3_C_init_uM = p(2);
	mlabfix_F_ = p(3);
	netValence_PIP2bindKCNQ = p(4);
	KCNQ_M_init_molecules_um_2 = p(5);
	Kr_PLCassoc = p(6);
	n_closed_init_molecules_um_2 = p(7);
	netValence_G2 = p(8);
	netValence_G1 = p(9);
	Kf_NE_RLG = p(10);
	Kr_PI4K_4Pase = p(11);
	V_leak = p(12);
	netValence_GTPase_Ga = p(13);
	netValence_PH_YFP_PIP2 = p(14);
	h_closed_init_molecules_um_2 = p(15);
	K_plc = p(16);
	IP3_PH_YFP_C_init_uM = p(17);
	netValence_Na_Gate_Act = p(18);
	C_M = p(19);
	RLG_GDP_M_init_molecules_um_2 = p(20);
	carrierValence_K_flux = p(21);
	RG_GDP_M_init_molecules_um_2 = p(22);
	DAG_M_init_molecules_um_2 = p(23);
	Kf_GTPase_Ga = p(24);
	KrG2 = p(25);
	GaGTP_M_init_molecules_um_2 = p(26);
	Kf_reconstitution = p(27);
	G_GDP_M_init_molecules_um_2 = p(28);
	theta = p(29);
	Kf_PLCdiss = p(30);
	Size_M = p(31);
	netValence_L2beta = p(32);
	Size_C = p(33);
	PIP2_M_init_molecules_um_2 = p(34);
	Kf_PI4K_4Pase = p(35);
	z_Kv = p(36);
	PH_YFP_PIP2_M_init_molecules_um_2 = p(37);
	KfG2 = p(38);
	Kv0 = p(39);
	carrierValence_Leak_flux = p(40);
	RGbeta_M_init_molecules_um_2 = p(41);
	carrierValence_KCNQ2_3 = p(42);
	netValence_PI4K_4Pase = p(43);
	Kf_NE_G = p(44);
	netValence_PLCdiss = p(45);
	PI4P_M_init_molecules_um_2 = p(46);
	PLC_basal = p(47);
	Kf_NE_RG = p(48);
	gKCNQ_max = p(49);
	V_K_closed = p(50);
	Ion_EX_init_uM = p(51);
	V_Na_Inact_close = p(52);
	netValence_GTPase_GaP = p(53);
	K_millivolts_per_volt = p(54);
	netValence_PLC_on_PI4P = p(55);
	mlabfix_K_GHK_ = p(56);
	Kg_KCNQ2_3 = p(57);
	netValence_PH_CFP_PIP2 = p(58);
	speed_PIP2_KCNQ = p(59);
	kp0 = p(60);
	mlabfix_F_nmol_ = p(61);
	n_open_init_molecules_um_2 = p(62);
	KD_PH_PIP2 = p(63);
	m_open_init_molecules_um_2 = p(64);
	m_closed_init_molecules_um_2 = p(65);
	IP3_PH_CFP_C_init_uM = p(66);
	K_IP3ase = p(67);
	netValence_PIP2hydr = p(68);
	PLC_efficiency_PIP = p(69);
	R_M_init_molecules_um_2 = p(70);
	mlabfix_T_ = p(71);
	netValence_NE_GaP = p(72);
	PI_M_init_molecules_um_2 = p(73);
	KCNQ_PIP2_M_init_molecules_um_2 = p(74);
	Ga_GDP_M_init_molecules_um_2 = p(75);
	KD_PH_IP3 = p(76);
	Kr_NE_G = p(77);
	PH_YFP_C_init_uM = p(78);
	mlabfix_R_ = p(79);
	V_Na_Act_close = p(80);
	Kf_PLCassoc = p(81);
	netValence_NE_RLG = p(82);
	h_open_init_molecules_um_2 = p(83);
	Hill_binding = p(84);
	GaGTP_PLC_M_init_molecules_um_2 = p(85);
	K_C_init_uM = p(86);
	netValence_NE_RG = p(87);
	RL_M_init_molecules_um_2 = p(88);
	GaGDP_PLC_M_init_molecules_um_2 = p(89);
	Kf_DAGPase = p(90);
	Na_C_init_uM = p(91);
	Size_EX = p(92);
	PH_CFP_PIP2_M_init_molecules_um_2 = p(93);
	K_EX_init_uM = p(94);
	gL = p(95);
	netValence_DAGPase = p(96);
	Voltage_M_init = p(97);
	Ion_C_init_uM = p(98);
	carrierValence_Na_Flux = p(99);
	V_K_open = p(100);
	speed_PH_PIP2 = p(101);
	mlabfix_PI_ = p(102);
	V_Na_Inact_open = p(103);
	theta_KCNQ2_3 = p(104);
	PH_CFP_C_init_uM = p(105);
	Kv = p(106);
	V_Na_Act_open = p(107);
	alpha = p(108);
	Kp = p(109);
	Kr_PLCdiss = p(110);
	Kg = p(111);
	oxoM_EX_init_uM = p(112);
	netValence_L2 = p(113);
	netValence_L1 = p(114);
	KA_PIP2_KCNQ = p(115);
	mlabfix_N_pmol_ = p(116);
	netValence_reconstitution = p(117);
	RLGbeta_M_init_molecules_um_2 = p(118);
	netValence_G2beta = p(119);
	Kr_NE_RG = p(120);
	KL1 = p(121);
	Na_EX_init_uM = p(122);
	gNa_max = p(123);
	netValence_K_Gate = p(124);
	KMOLE = p(125);
	netValence_NE_G = p(126);
	Gbeta_M_init_molecules_um_2 = p(127);
	gK_max = p(128);
	netValence_G1beta = p(129);
	netValence_PLCassoc = p(130);
	KrL1 = p(131);
	PLC_M_init_molecules_um_2 = p(132);
	speed_PH_IP3 = p(133);
	% Functions
	alpha_n = (0.01 .* (Voltage_M - V_K_open) ./ (1.0 - exp(( - 0.1 .* (Voltage_M - V_K_open)))));
	alpha_m = (0.1 .* (Voltage_M - V_Na_Act_open) ./ (1.0 - exp(( - 0.1 .* (Voltage_M - V_Na_Act_open)))));
	J_PI4K_4Pase = ((Kf_PI4K_4Pase .* PI_M) - (Kr_PI4K_4Pase .* PI4P_M));
	alpha_h = (0.07 .* exp( - (0.05 .* (Voltage_M - V_Na_Inact_open))));
	UnitFactor_uM_um3_molecules_neg_1 = (1.0 ./ 602.0);
	K_RG_GDP_M_total = ((Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RLGbeta_M_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RL_M_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* R_M_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RGbeta_M_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RG_GDP_M_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RLG_GDP_M_init_molecules_um_2));
	RG_GDP_M = (( - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RLGbeta_M) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RL_M) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* R_M) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RGbeta_M) + K_RG_GDP_M_total - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RLG_GDP_M)) ./ (Size_M .* UnitFactor_uM_um3_molecules_neg_1));
	J_NE_RG = ((Kf_NE_RG .* RG_GDP_M) - ((Kr_NE_RG .* GaGTP_M) .* RGbeta_M));
	unitFactor_Na_Flux = (1.0E9 ./ 1.0);
	K_G_GDP_M_total = ((Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* Gbeta_M_init_molecules_um_2) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RL_M_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* G_GDP_M_init_molecules_um_2) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* R_M_init_molecules_um_2));
	unitFactor_K_flux = (1.0E9 ./ 1.0);
	Kr_PH_YFP_PIP2 = (KD_PH_PIP2 .* speed_PH_PIP2);
	J_IP3Pase = (K_IP3ase .* IP3_C);
	Kf_G2 = KfG2;
	KfG1 = (0.1 .* KfG2);
	Kf_G1 = KfG1;
	K_PLC_M_total = ((Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* Gbeta_M_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RLGbeta_M_init_molecules_um_2) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* GaGTP_M_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RGbeta_M_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* PLC_M_init_molecules_um_2) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* Ga_GDP_M_init_molecules_um_2));
	PLC_M = (( - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* Gbeta_M) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RLGbeta_M) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* GaGTP_M) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RGbeta_M) + K_PLC_M_total + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* Ga_GDP_M)) ./ (Size_M .* UnitFactor_uM_um3_molecules_neg_1));
	J_PLCassoc = (((Kf_PLCassoc .* GaGTP_M) .* PLC_M) - (Kr_PLCassoc .* GaGTP_PLC_M));
	KG2 = (KrG2 ./ KfG2);
	KG1 = (KG2 .* alpha);
	K_C = K_C_init_uM;
	K_m_closed_total = ((Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* m_closed_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* m_open_init_molecules_um_2));
	Na_EX = Na_EX_init_uM;
	K_oxoM_EX_total = ((Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RLGbeta_M_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RL_M_init_molecules_um_2) + (Size_EX .* oxoM_EX_init_uM) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RLG_GDP_M_init_molecules_um_2));
	oxoM_EX = (( - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RLGbeta_M) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RL_M) + K_oxoM_EX_total - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RLG_GDP_M)) ./ Size_EX);
	KfL1 = (KrL1 ./ KL1);
	KfL2 = KfL1;
	Kf_L2beta = KfL2;
	KL2 = (KL1 ./ alpha);
	KrL2 = (KL2 .* KfL2);
	Kr_L2beta = KrL2;
	J_L2beta = (((Kf_L2beta .* RGbeta_M) .* oxoM_EX) - (Kr_L2beta .* RLGbeta_M));
	Kr_G2beta = KrG2;
	Kf_PH_YFP_IP3 = speed_PH_IP3;
	Kr_PH_YFP_IP3 = (KD_PH_IP3 .* speed_PH_IP3);
	K_PH_YFP_C_total = ( - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* PI4P_M_init_molecules_um_2) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* PI_M_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* KCNQ_M_init_molecules_um_2) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* PH_CFP_PIP2_M_init_molecules_um_2) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* PIP2_M_init_molecules_um_2) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* DAG_M_init_molecules_um_2) + (Size_C .* PH_YFP_C_init_uM) + (Size_C .* IP3_PH_YFP_C_init_uM));
	PH_YFP_C = (((Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* PI4P_M) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* PI_M) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* KCNQ_M) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* PH_CFP_PIP2_M) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* PIP2_M) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* DAG_M) + K_PH_YFP_C_total - (Size_C .* IP3_PH_YFP_C)) ./ Size_C);
	J_PH_YFP_IP3 = (((Kf_PH_YFP_IP3 .* IP3_C) .* PH_YFP_C) - (Kr_PH_YFP_IP3 .* IP3_PH_YFP_C));
	K_PH_YFP_PIP2_M_total = ((Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* PI4P_M_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* PI_M_init_molecules_um_2) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* KCNQ_M_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* PH_CFP_PIP2_M_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* PIP2_M_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* PH_YFP_PIP2_M_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* DAG_M_init_molecules_um_2));
	KrG1 = (KfG1 .* KG1);
	K_EX = K_EX_init_uM;
	Kv_KCNQ2_3 = (Kv0 .* exp((z_Kv .* mlabfix_F_ .* Voltage_M ./ (mlabfix_R_ .* mlabfix_T_))));
	K_KCNQ_PIP2_M_total = ((Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* KCNQ_PIP2_M_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* KCNQ_M_init_molecules_um_2));
	KCNQ_PIP2_M = ((K_KCNQ_PIP2_M_total - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* KCNQ_M)) ./ (Size_M .* UnitFactor_uM_um3_molecules_neg_1));
	Kp_KCNQ2_3 = (kp0 .* KCNQ_PIP2_M);
	PPo = (((Kg_KCNQ2_3 + (Kp_KCNQ2_3 .* Kg_KCNQ2_3) + (theta_KCNQ2_3 .* Kv_KCNQ2_3 .* Kp_KCNQ2_3 .* Kg_KCNQ2_3)) ./ (1.0 + Kv_KCNQ2_3 + Kg_KCNQ2_3 + Kp_KCNQ2_3 + (Kv_KCNQ2_3 .* Kg_KCNQ2_3) + (Kv_KCNQ2_3 .* Kp_KCNQ2_3) + (Kp_KCNQ2_3 .* Kg_KCNQ2_3) + (theta_KCNQ2_3 .* Kv_KCNQ2_3 .* Kp_KCNQ2_3 .* Kg_KCNQ2_3))) ./ ((Kg_KCNQ2_3 + (theta_KCNQ2_3 .* Kp_KCNQ2_3 .* Kg_KCNQ2_3)) ./ (1.0 + Kg_KCNQ2_3 + Kp_KCNQ2_3 + (theta_KCNQ2_3 .* Kp_KCNQ2_3 .* Kg_KCNQ2_3))));
	C_KCNQ2_3 = (gKCNQ_max .* PPo);
	I_KCNQ2_3 = (C_KCNQ2_3 .* ((((mlabfix_R_ .* mlabfix_T_) ./ (carrierValence_KCNQ2_3 .* mlabfix_F_)) .* log((K_EX ./ K_C))) - Voltage_M));
	unitFactor_KCNQ2_3 = (1.0E9 ./ 1.0);
	J_KCNQ2_3 = ((I_KCNQ2_3 ./ (carrierValence_KCNQ2_3 .* mlabfix_F_)) .* unitFactor_KCNQ2_3);
	K_GaGDP_PLC_M_total = ((Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* GaGTP_PLC_M_init_molecules_um_2) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* Gbeta_M_init_molecules_um_2) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RLGbeta_M_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* GaGTP_M_init_molecules_um_2) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RGbeta_M_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* GaGDP_PLC_M_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* Ga_GDP_M_init_molecules_um_2));
	PH_YFP_PIP2_M = (( - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* PI4P_M) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* PI_M) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* KCNQ_M) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* PH_CFP_PIP2_M) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* PIP2_M) + K_PH_YFP_PIP2_M_total - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* DAG_M)) ./ (Size_M .* UnitFactor_uM_um3_molecules_neg_1));
	Kr_G1beta = KrG1;
	UnitFactor_mV_pF_s_neg_1_pA_neg_1 = (1000.0 ./ 1.0);
	G_GDP_M = (( - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* Gbeta_M) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RL_M) + K_G_GDP_M_total + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* R_M)) ./ (Size_M .* UnitFactor_uM_um3_molecules_neg_1));
	J_NE_G = ((Kf_NE_G .* G_GDP_M) - ((Kr_NE_G .* Gbeta_M) .* GaGTP_M));
	KFlux_M_C = (Size_M ./ Size_C);
	unitFactor_Leak_flux = (1.0E9 ./ 1.0);
	Kf_PIP2bindKCNQ = speed_PIP2_KCNQ;
	KD_PIP2_KCNQ = (KA_PIP2_KCNQ ^ Hill_binding);
	Kr_PIP2bindKCNQ = (KD_PIP2_KCNQ .* Kf_PIP2bindKCNQ);
	J_PIP2bindKCNQ = ((Kf_PIP2bindKCNQ .* (PIP2_M ^ Hill_binding) .* KCNQ_M) - (Kr_PIP2bindKCNQ .* KCNQ_PIP2_M));
	K_IP3_PH_CFP_C_total = ((Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* PH_CFP_PIP2_M_init_molecules_um_2) + (Size_C .* IP3_PH_CFP_C_init_uM) + (Size_C .* PH_CFP_C_init_uM));
	IP3_PH_CFP_C = (( - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* PH_CFP_PIP2_M) + K_IP3_PH_CFP_C_total - (Size_C .* PH_CFP_C)) ./ Size_C);
	Na_C = Na_C_init_uM;
	C_Na_Flux = (gNa_max .* (m_open ^ 3.0) .* h_open);
	I_Na_Flux = (C_Na_Flux .* ((((mlabfix_R_ .* mlabfix_T_) ./ (carrierValence_Na_Flux .* mlabfix_F_)) .* log((Na_EX ./ Na_C))) - Voltage_M));
	Kr_PH_CFP_PIP2 = (KD_PH_PIP2 .* speed_PH_PIP2);
	Kf_PH_YFP_PIP2 = speed_PH_PIP2;
	GaGDP_PLC_M = (( - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* GaGTP_PLC_M) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* Gbeta_M) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RLGbeta_M) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* GaGTP_M) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RGbeta_M) + K_GaGDP_PLC_M_total - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* Ga_GDP_M)) ./ (Size_M .* UnitFactor_uM_um3_molecules_neg_1));
	allGaPLC = (GaGTP_PLC_M + GaGDP_PLC_M);
	beta_n = (0.125 .* exp(( - 0.0125 .* (Voltage_M - V_K_closed))));
	beta_m = (4.0 .* exp(( - (Voltage_M - V_Na_Act_close) ./ 18.0)));
	device_totalCurrClampElectrode.F = (100.0 .* (((t > 10.0) && (t < 20.0)) || ((t > 2000.0) && (t < 2010.0))));
	device_totalCurrClampElectrode.I = device_totalCurrClampElectrode.F;
	beta_h = (1.0 ./ (1.0 + exp(( - 0.1 .* (Voltage_M - V_Na_Inact_close)))));
	Kf_PH_CFP_IP3 = speed_PH_IP3;
	Kr_PH_CFP_IP3 = (KD_PH_IP3 .* speed_PH_IP3);
	J_PH_CFP_IP3 = (((Kf_PH_CFP_IP3 .* IP3_C) .* PH_CFP_C) - (Kr_PH_CFP_IP3 .* IP3_PH_CFP_C));
	K_h_closed_total = ((Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* h_closed_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* h_open_init_molecules_um_2));
	h_closed = ((K_h_closed_total - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* h_open)) ./ (Size_M .* UnitFactor_uM_um3_molecules_neg_1));
	I_Leak_flux = (gL .* (Voltage_M - V_leak));
	C_K_flux = (gK_max .* (n_open ^ 4.0));
	I_K_flux = (C_K_flux .* ((((mlabfix_R_ .* mlabfix_T_) ./ (carrierValence_K_flux .* mlabfix_F_)) .* log((K_EX ./ K_C))) - Voltage_M));
	F_M = ( - (I_Leak_flux .* Size_M) - (I_Na_Flux .* Size_M) - (I_KCNQ2_3 .* Size_M) - (I_K_flux .* Size_M));
	J_GTPase_GaP = (15.0 .* GaGTP_PLC_M);
	Kf_G2beta = KfG2;
	KFlux_M_EX = (Size_M ./ Size_EX);
	J_PH_YFP_PIP2 = (((Kf_PH_YFP_PIP2 .* PH_YFP_C) .* PIP2_M) - (Kr_PH_YFP_PIP2 .* PH_YFP_PIP2_M));
	m_closed = ((K_m_closed_total - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* m_open)) ./ (Size_M .* UnitFactor_uM_um3_molecules_neg_1));
	J_Na_Gate_Act = ((alpha_m .* m_closed) - (beta_m .* m_open));
	J_NE_GaP = (4.7 .* GaGDP_PLC_M);
	J_PLCdiss = ((Kf_PLCdiss .* GaGDP_PLC_M) - ((Kr_PLCdiss .* PLC_M) .* Ga_GDP_M));
	Kf_G1beta = KfG1;
	J_Leak_flux = ((I_Leak_flux ./ (carrierValence_Leak_flux .* mlabfix_F_)) .* unitFactor_Leak_flux);
	J_NE_RLG = (Kf_NE_RLG .* RLG_GDP_M);
	J_PIP2hydr = (PIP2_M .* (PLC_basal + (K_plc .* GaGTP_PLC_M)));
	I_M = device_totalCurrClampElectrode.F;
	Kr_L2 = KrL2;
	Kf_L2 = KfL2;
	J_L2 = (((Kf_L2 .* oxoM_EX) .* RG_GDP_M) - (Kr_L2 .* RLG_GDP_M));
	J_PLC_on_PI4P = (PI4P_M .* PLC_efficiency_PIP .* (PLC_basal + (GaGTP_PLC_M .* K_plc)));
	Kf_L1 = KfL1;
	Kr_L1 = KrL1;
	J_L1 = (((Kf_L1 .* oxoM_EX) .* R_M) - (Kr_L1 .* RL_M));
	Kf_PH_CFP_PIP2 = speed_PH_PIP2;
	J_PH_CFP_PIP2 = (((Kf_PH_CFP_PIP2 .* PIP2_M) .* PH_CFP_C) - (Kr_PH_CFP_PIP2 .* PH_CFP_PIP2_M));
	J_reconstitution = (Kf_reconstitution .* Gbeta_M .* Ga_GDP_M);
	J_Na_Flux = ((I_Na_Flux ./ (carrierValence_Na_Flux .* mlabfix_F_)) .* unitFactor_Na_Flux);
	device_M.Capacitance = (C_M .* Size_M);
	allRL = (RL_M + RLG_GDP_M + RLGbeta_M);
	J_Na_Gate_Inact = ((alpha_h .* h_closed) - (beta_h .* h_open));
	Kr_G2 = KrG2;
	J_G2 = (((Kf_G2 .* RL_M) .* G_GDP_M) - (Kr_G2 .* RLG_GDP_M));
	Kr_G1 = KrG1;
	J_G1 = (((Kf_G1 .* G_GDP_M) .* R_M) - (Kr_G1 .* RG_GDP_M));
	K_n_closed_total = ((Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* n_closed_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* n_open_init_molecules_um_2));
	Ion_EX = Ion_EX_init_uM;
	allPH_cytoplasm = (PH_CFP_C + PH_YFP_C + IP3_PH_CFP_C + IP3_PH_YFP_C);
	J_G2beta = (((Kf_G2beta .* RL_M) .* Gbeta_M) - (Kr_G2beta .* RLGbeta_M));
	Ion_C = Ion_C_init_uM;
	n_closed = ((K_n_closed_total - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* n_open)) ./ (Size_M .* UnitFactor_uM_um3_molecules_neg_1));
	J_K_Gate = ((alpha_n .* n_closed) - (beta_n .* n_open));
	allPIP2M = (PIP2_M + KCNQ_PIP2_M + PH_CFP_PIP2_M + PH_YFP_PIP2_M);
	J_G1beta = (((Kf_G1beta .* Gbeta_M) .* R_M) - (Kr_G1beta .* RGbeta_M));
	allGabg = (G_GDP_M + RG_GDP_M + RLG_GDP_M);
	J_DAGPase = (Kf_DAGPase .* DAG_M);
	J_GTPase_Ga = (Kf_GTPase_Ga .* GaGTP_M);
	J_K_flux = ((I_K_flux ./ (carrierValence_K_flux .* mlabfix_F_)) .* unitFactor_K_flux);
	allRGbeta = (RG_GDP_M + RLG_GDP_M + RLGbeta_M + RGbeta_M);

	rowValue = [PIP2_M h_open Ga_GDP_M RLGbeta_M PH_CFP_C n_open RLG_GDP_M R_M PI_M PH_CFP_PIP2_M m_open Gbeta_M RGbeta_M GaGTP_PLC_M DAG_M IP3_PH_YFP_C GaGTP_M IP3_C Voltage_M RL_M PI4P_M KCNQ_M alpha_n alpha_m J_PI4K_4Pase alpha_h UnitFactor_uM_um3_molecules_neg_1 K_RG_GDP_M_total RG_GDP_M J_NE_RG unitFactor_Na_Flux K_G_GDP_M_total unitFactor_K_flux Kr_PH_YFP_PIP2 J_IP3Pase Kf_G2 KfG1 Kf_G1 K_PLC_M_total PLC_M J_PLCassoc KG2 KG1 K_C K_m_closed_total Na_EX K_oxoM_EX_total oxoM_EX KfL1 KfL2 Kf_L2beta KL2 KrL2 Kr_L2beta J_L2beta Kr_G2beta Kf_PH_YFP_IP3 Kr_PH_YFP_IP3 K_PH_YFP_C_total PH_YFP_C J_PH_YFP_IP3 K_PH_YFP_PIP2_M_total KrG1 K_EX Kv_KCNQ2_3 K_KCNQ_PIP2_M_total KCNQ_PIP2_M Kp_KCNQ2_3 PPo C_KCNQ2_3 I_KCNQ2_3 unitFactor_KCNQ2_3 J_KCNQ2_3 K_GaGDP_PLC_M_total PH_YFP_PIP2_M Kr_G1beta UnitFactor_mV_pF_s_neg_1_pA_neg_1 G_GDP_M J_NE_G KFlux_M_C unitFactor_Leak_flux Kf_PIP2bindKCNQ KD_PIP2_KCNQ Kr_PIP2bindKCNQ J_PIP2bindKCNQ K_IP3_PH_CFP_C_total IP3_PH_CFP_C Na_C C_Na_Flux I_Na_Flux Kr_PH_CFP_PIP2 Kf_PH_YFP_PIP2 GaGDP_PLC_M allGaPLC beta_n beta_m device_totalCurrClampElectrode.F device_totalCurrClampElectrode.I beta_h Kf_PH_CFP_IP3 Kr_PH_CFP_IP3 J_PH_CFP_IP3 K_h_closed_total h_closed I_Leak_flux C_K_flux I_K_flux F_M J_GTPase_GaP Kf_G2beta KFlux_M_EX J_PH_YFP_PIP2 m_closed J_Na_Gate_Act J_NE_GaP J_PLCdiss Kf_G1beta J_Leak_flux J_NE_RLG J_PIP2hydr I_M Kr_L2 Kf_L2 J_L2 J_PLC_on_PI4P Kf_L1 Kr_L1 J_L1 Kf_PH_CFP_PIP2 J_PH_CFP_PIP2 J_reconstitution J_Na_Flux device_M.Capacitance allRL J_Na_Gate_Inact Kr_G2 J_G2 Kr_G1 J_G1 K_n_closed_total Ion_EX allPH_cytoplasm J_G2beta Ion_C n_closed J_K_Gate allPIP2M J_G1beta allGabg J_DAGPase J_GTPase_Ga J_K_flux allRGbeta ];
end

% -------------------------------------------------------
% ode rate
function dydt = f(t,y,p,y0)
	% State Variables
	PIP2_M = y(1);
	h_open = y(2);
	Ga_GDP_M = y(3);
	RLGbeta_M = y(4);
	PH_CFP_C = y(5);
	n_open = y(6);
	RLG_GDP_M = y(7);
	R_M = y(8);
	PI_M = y(9);
	PH_CFP_PIP2_M = y(10);
	m_open = y(11);
	Gbeta_M = y(12);
	RGbeta_M = y(13);
	GaGTP_PLC_M = y(14);
	DAG_M = y(15);
	IP3_PH_YFP_C = y(16);
	GaGTP_M = y(17);
	IP3_C = y(18);
	Voltage_M = y(19);
	RL_M = y(20);
	PI4P_M = y(21);
	KCNQ_M = y(22);
	% Constants
	netValence_Na_Gate_Inact = p(1);
	IP3_C_init_uM = p(2);
	mlabfix_F_ = p(3);
	netValence_PIP2bindKCNQ = p(4);
	KCNQ_M_init_molecules_um_2 = p(5);
	Kr_PLCassoc = p(6);
	n_closed_init_molecules_um_2 = p(7);
	netValence_G2 = p(8);
	netValence_G1 = p(9);
	Kf_NE_RLG = p(10);
	Kr_PI4K_4Pase = p(11);
	V_leak = p(12);
	netValence_GTPase_Ga = p(13);
	netValence_PH_YFP_PIP2 = p(14);
	h_closed_init_molecules_um_2 = p(15);
	K_plc = p(16);
	IP3_PH_YFP_C_init_uM = p(17);
	netValence_Na_Gate_Act = p(18);
	C_M = p(19);
	RLG_GDP_M_init_molecules_um_2 = p(20);
	carrierValence_K_flux = p(21);
	RG_GDP_M_init_molecules_um_2 = p(22);
	DAG_M_init_molecules_um_2 = p(23);
	Kf_GTPase_Ga = p(24);
	KrG2 = p(25);
	GaGTP_M_init_molecules_um_2 = p(26);
	Kf_reconstitution = p(27);
	G_GDP_M_init_molecules_um_2 = p(28);
	theta = p(29);
	Kf_PLCdiss = p(30);
	Size_M = p(31);
	netValence_L2beta = p(32);
	Size_C = p(33);
	PIP2_M_init_molecules_um_2 = p(34);
	Kf_PI4K_4Pase = p(35);
	z_Kv = p(36);
	PH_YFP_PIP2_M_init_molecules_um_2 = p(37);
	KfG2 = p(38);
	Kv0 = p(39);
	carrierValence_Leak_flux = p(40);
	RGbeta_M_init_molecules_um_2 = p(41);
	carrierValence_KCNQ2_3 = p(42);
	netValence_PI4K_4Pase = p(43);
	Kf_NE_G = p(44);
	netValence_PLCdiss = p(45);
	PI4P_M_init_molecules_um_2 = p(46);
	PLC_basal = p(47);
	Kf_NE_RG = p(48);
	gKCNQ_max = p(49);
	V_K_closed = p(50);
	Ion_EX_init_uM = p(51);
	V_Na_Inact_close = p(52);
	netValence_GTPase_GaP = p(53);
	K_millivolts_per_volt = p(54);
	netValence_PLC_on_PI4P = p(55);
	mlabfix_K_GHK_ = p(56);
	Kg_KCNQ2_3 = p(57);
	netValence_PH_CFP_PIP2 = p(58);
	speed_PIP2_KCNQ = p(59);
	kp0 = p(60);
	mlabfix_F_nmol_ = p(61);
	n_open_init_molecules_um_2 = p(62);
	KD_PH_PIP2 = p(63);
	m_open_init_molecules_um_2 = p(64);
	m_closed_init_molecules_um_2 = p(65);
	IP3_PH_CFP_C_init_uM = p(66);
	K_IP3ase = p(67);
	netValence_PIP2hydr = p(68);
	PLC_efficiency_PIP = p(69);
	R_M_init_molecules_um_2 = p(70);
	mlabfix_T_ = p(71);
	netValence_NE_GaP = p(72);
	PI_M_init_molecules_um_2 = p(73);
	KCNQ_PIP2_M_init_molecules_um_2 = p(74);
	Ga_GDP_M_init_molecules_um_2 = p(75);
	KD_PH_IP3 = p(76);
	Kr_NE_G = p(77);
	PH_YFP_C_init_uM = p(78);
	mlabfix_R_ = p(79);
	V_Na_Act_close = p(80);
	Kf_PLCassoc = p(81);
	netValence_NE_RLG = p(82);
	h_open_init_molecules_um_2 = p(83);
	Hill_binding = p(84);
	GaGTP_PLC_M_init_molecules_um_2 = p(85);
	K_C_init_uM = p(86);
	netValence_NE_RG = p(87);
	RL_M_init_molecules_um_2 = p(88);
	GaGDP_PLC_M_init_molecules_um_2 = p(89);
	Kf_DAGPase = p(90);
	Na_C_init_uM = p(91);
	Size_EX = p(92);
	PH_CFP_PIP2_M_init_molecules_um_2 = p(93);
	K_EX_init_uM = p(94);
	gL = p(95);
	netValence_DAGPase = p(96);
	Voltage_M_init = p(97);
	Ion_C_init_uM = p(98);
	carrierValence_Na_Flux = p(99);
	V_K_open = p(100);
	speed_PH_PIP2 = p(101);
	mlabfix_PI_ = p(102);
	V_Na_Inact_open = p(103);
	theta_KCNQ2_3 = p(104);
	PH_CFP_C_init_uM = p(105);
	Kv = p(106);
	V_Na_Act_open = p(107);
	alpha = p(108);
	Kp = p(109);
	Kr_PLCdiss = p(110);
	Kg = p(111);
	oxoM_EX_init_uM = p(112);
	netValence_L2 = p(113);
	netValence_L1 = p(114);
	KA_PIP2_KCNQ = p(115);
	mlabfix_N_pmol_ = p(116);
	netValence_reconstitution = p(117);
	RLGbeta_M_init_molecules_um_2 = p(118);
	netValence_G2beta = p(119);
	Kr_NE_RG = p(120);
	KL1 = p(121);
	Na_EX_init_uM = p(122);
	gNa_max = p(123);
	netValence_K_Gate = p(124);
	KMOLE = p(125);
	netValence_NE_G = p(126);
	Gbeta_M_init_molecules_um_2 = p(127);
	gK_max = p(128);
	netValence_G1beta = p(129);
	netValence_PLCassoc = p(130);
	KrL1 = p(131);
	PLC_M_init_molecules_um_2 = p(132);
	speed_PH_IP3 = p(133);
	% Functions
	alpha_n = (0.01 .* (Voltage_M - V_K_open) ./ (1.0 - exp(( - 0.1 .* (Voltage_M - V_K_open)))));
	alpha_m = (0.1 .* (Voltage_M - V_Na_Act_open) ./ (1.0 - exp(( - 0.1 .* (Voltage_M - V_Na_Act_open)))));
	J_PI4K_4Pase = ((Kf_PI4K_4Pase .* PI_M) - (Kr_PI4K_4Pase .* PI4P_M));
	alpha_h = (0.07 .* exp( - (0.05 .* (Voltage_M - V_Na_Inact_open))));
	UnitFactor_uM_um3_molecules_neg_1 = (1.0 ./ 602.0);
	K_RG_GDP_M_total = ((Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RLGbeta_M_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RL_M_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* R_M_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RGbeta_M_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RG_GDP_M_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RLG_GDP_M_init_molecules_um_2));
	RG_GDP_M = (( - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RLGbeta_M) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RL_M) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* R_M) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RGbeta_M) + K_RG_GDP_M_total - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RLG_GDP_M)) ./ (Size_M .* UnitFactor_uM_um3_molecules_neg_1));
	J_NE_RG = ((Kf_NE_RG .* RG_GDP_M) - ((Kr_NE_RG .* GaGTP_M) .* RGbeta_M));
	unitFactor_Na_Flux = (1.0E9 ./ 1.0);
	K_G_GDP_M_total = ((Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* Gbeta_M_init_molecules_um_2) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RL_M_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* G_GDP_M_init_molecules_um_2) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* R_M_init_molecules_um_2));
	unitFactor_K_flux = (1.0E9 ./ 1.0);
	Kr_PH_YFP_PIP2 = (KD_PH_PIP2 .* speed_PH_PIP2);
	J_IP3Pase = (K_IP3ase .* IP3_C);
	Kf_G2 = KfG2;
	KfG1 = (0.1 .* KfG2);
	Kf_G1 = KfG1;
	K_PLC_M_total = ((Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* Gbeta_M_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RLGbeta_M_init_molecules_um_2) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* GaGTP_M_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RGbeta_M_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* PLC_M_init_molecules_um_2) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* Ga_GDP_M_init_molecules_um_2));
	PLC_M = (( - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* Gbeta_M) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RLGbeta_M) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* GaGTP_M) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RGbeta_M) + K_PLC_M_total + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* Ga_GDP_M)) ./ (Size_M .* UnitFactor_uM_um3_molecules_neg_1));
	J_PLCassoc = (((Kf_PLCassoc .* GaGTP_M) .* PLC_M) - (Kr_PLCassoc .* GaGTP_PLC_M));
	KG2 = (KrG2 ./ KfG2);
	KG1 = (KG2 .* alpha);
	K_C = K_C_init_uM;
	K_m_closed_total = ((Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* m_closed_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* m_open_init_molecules_um_2));
	Na_EX = Na_EX_init_uM;
	K_oxoM_EX_total = ((Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RLGbeta_M_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RL_M_init_molecules_um_2) + (Size_EX .* oxoM_EX_init_uM) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RLG_GDP_M_init_molecules_um_2));
	oxoM_EX = (( - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RLGbeta_M) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RL_M) + K_oxoM_EX_total - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RLG_GDP_M)) ./ Size_EX);
	KfL1 = (KrL1 ./ KL1);
	KfL2 = KfL1;
	Kf_L2beta = KfL2;
	KL2 = (KL1 ./ alpha);
	KrL2 = (KL2 .* KfL2);
	Kr_L2beta = KrL2;
	J_L2beta = (((Kf_L2beta .* RGbeta_M) .* oxoM_EX) - (Kr_L2beta .* RLGbeta_M));
	Kr_G2beta = KrG2;
	Kf_PH_YFP_IP3 = speed_PH_IP3;
	Kr_PH_YFP_IP3 = (KD_PH_IP3 .* speed_PH_IP3);
	K_PH_YFP_C_total = ( - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* PI4P_M_init_molecules_um_2) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* PI_M_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* KCNQ_M_init_molecules_um_2) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* PH_CFP_PIP2_M_init_molecules_um_2) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* PIP2_M_init_molecules_um_2) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* DAG_M_init_molecules_um_2) + (Size_C .* PH_YFP_C_init_uM) + (Size_C .* IP3_PH_YFP_C_init_uM));
	PH_YFP_C = (((Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* PI4P_M) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* PI_M) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* KCNQ_M) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* PH_CFP_PIP2_M) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* PIP2_M) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* DAG_M) + K_PH_YFP_C_total - (Size_C .* IP3_PH_YFP_C)) ./ Size_C);
	J_PH_YFP_IP3 = (((Kf_PH_YFP_IP3 .* IP3_C) .* PH_YFP_C) - (Kr_PH_YFP_IP3 .* IP3_PH_YFP_C));
	K_PH_YFP_PIP2_M_total = ((Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* PI4P_M_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* PI_M_init_molecules_um_2) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* KCNQ_M_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* PH_CFP_PIP2_M_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* PIP2_M_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* PH_YFP_PIP2_M_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* DAG_M_init_molecules_um_2));
	KrG1 = (KfG1 .* KG1);
	K_EX = K_EX_init_uM;
	Kv_KCNQ2_3 = (Kv0 .* exp((z_Kv .* mlabfix_F_ .* Voltage_M ./ (mlabfix_R_ .* mlabfix_T_))));
	K_KCNQ_PIP2_M_total = ((Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* KCNQ_PIP2_M_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* KCNQ_M_init_molecules_um_2));
	KCNQ_PIP2_M = ((K_KCNQ_PIP2_M_total - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* KCNQ_M)) ./ (Size_M .* UnitFactor_uM_um3_molecules_neg_1));
	Kp_KCNQ2_3 = (kp0 .* KCNQ_PIP2_M);
	PPo = (((Kg_KCNQ2_3 + (Kp_KCNQ2_3 .* Kg_KCNQ2_3) + (theta_KCNQ2_3 .* Kv_KCNQ2_3 .* Kp_KCNQ2_3 .* Kg_KCNQ2_3)) ./ (1.0 + Kv_KCNQ2_3 + Kg_KCNQ2_3 + Kp_KCNQ2_3 + (Kv_KCNQ2_3 .* Kg_KCNQ2_3) + (Kv_KCNQ2_3 .* Kp_KCNQ2_3) + (Kp_KCNQ2_3 .* Kg_KCNQ2_3) + (theta_KCNQ2_3 .* Kv_KCNQ2_3 .* Kp_KCNQ2_3 .* Kg_KCNQ2_3))) ./ ((Kg_KCNQ2_3 + (theta_KCNQ2_3 .* Kp_KCNQ2_3 .* Kg_KCNQ2_3)) ./ (1.0 + Kg_KCNQ2_3 + Kp_KCNQ2_3 + (theta_KCNQ2_3 .* Kp_KCNQ2_3 .* Kg_KCNQ2_3))));
	C_KCNQ2_3 = (gKCNQ_max .* PPo);
	I_KCNQ2_3 = (C_KCNQ2_3 .* ((((mlabfix_R_ .* mlabfix_T_) ./ (carrierValence_KCNQ2_3 .* mlabfix_F_)) .* log((K_EX ./ K_C))) - Voltage_M));
	unitFactor_KCNQ2_3 = (1.0E9 ./ 1.0);
	J_KCNQ2_3 = ((I_KCNQ2_3 ./ (carrierValence_KCNQ2_3 .* mlabfix_F_)) .* unitFactor_KCNQ2_3);
	K_GaGDP_PLC_M_total = ((Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* GaGTP_PLC_M_init_molecules_um_2) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* Gbeta_M_init_molecules_um_2) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RLGbeta_M_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* GaGTP_M_init_molecules_um_2) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RGbeta_M_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* GaGDP_PLC_M_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* Ga_GDP_M_init_molecules_um_2));
	PH_YFP_PIP2_M = (( - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* PI4P_M) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* PI_M) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* KCNQ_M) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* PH_CFP_PIP2_M) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* PIP2_M) + K_PH_YFP_PIP2_M_total - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* DAG_M)) ./ (Size_M .* UnitFactor_uM_um3_molecules_neg_1));
	Kr_G1beta = KrG1;
	UnitFactor_mV_pF_s_neg_1_pA_neg_1 = (1000.0 ./ 1.0);
	G_GDP_M = (( - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* Gbeta_M) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RL_M) + K_G_GDP_M_total + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* R_M)) ./ (Size_M .* UnitFactor_uM_um3_molecules_neg_1));
	J_NE_G = ((Kf_NE_G .* G_GDP_M) - ((Kr_NE_G .* Gbeta_M) .* GaGTP_M));
	KFlux_M_C = (Size_M ./ Size_C);
	unitFactor_Leak_flux = (1.0E9 ./ 1.0);
	Kf_PIP2bindKCNQ = speed_PIP2_KCNQ;
	KD_PIP2_KCNQ = (KA_PIP2_KCNQ ^ Hill_binding);
	Kr_PIP2bindKCNQ = (KD_PIP2_KCNQ .* Kf_PIP2bindKCNQ);
	J_PIP2bindKCNQ = ((Kf_PIP2bindKCNQ .* (PIP2_M ^ Hill_binding) .* KCNQ_M) - (Kr_PIP2bindKCNQ .* KCNQ_PIP2_M));
	K_IP3_PH_CFP_C_total = ((Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* PH_CFP_PIP2_M_init_molecules_um_2) + (Size_C .* IP3_PH_CFP_C_init_uM) + (Size_C .* PH_CFP_C_init_uM));
	IP3_PH_CFP_C = (( - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* PH_CFP_PIP2_M) + K_IP3_PH_CFP_C_total - (Size_C .* PH_CFP_C)) ./ Size_C);
	Na_C = Na_C_init_uM;
	C_Na_Flux = (gNa_max .* (m_open ^ 3.0) .* h_open);
	I_Na_Flux = (C_Na_Flux .* ((((mlabfix_R_ .* mlabfix_T_) ./ (carrierValence_Na_Flux .* mlabfix_F_)) .* log((Na_EX ./ Na_C))) - Voltage_M));
	Kr_PH_CFP_PIP2 = (KD_PH_PIP2 .* speed_PH_PIP2);
	Kf_PH_YFP_PIP2 = speed_PH_PIP2;
	GaGDP_PLC_M = (( - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* GaGTP_PLC_M) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* Gbeta_M) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RLGbeta_M) - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* GaGTP_M) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* RGbeta_M) + K_GaGDP_PLC_M_total - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* Ga_GDP_M)) ./ (Size_M .* UnitFactor_uM_um3_molecules_neg_1));
	allGaPLC = (GaGTP_PLC_M + GaGDP_PLC_M);
	beta_n = (0.125 .* exp(( - 0.0125 .* (Voltage_M - V_K_closed))));
	beta_m = (4.0 .* exp(( - (Voltage_M - V_Na_Act_close) ./ 18.0)));
	device_totalCurrClampElectrode.F = (100.0 .* (((t > 10.0) && (t < 20.0)) || ((t > 2000.0) && (t < 2010.0))));
	device_totalCurrClampElectrode.I = device_totalCurrClampElectrode.F;
	beta_h = (1.0 ./ (1.0 + exp(( - 0.1 .* (Voltage_M - V_Na_Inact_close)))));
	Kf_PH_CFP_IP3 = speed_PH_IP3;
	Kr_PH_CFP_IP3 = (KD_PH_IP3 .* speed_PH_IP3);
	J_PH_CFP_IP3 = (((Kf_PH_CFP_IP3 .* IP3_C) .* PH_CFP_C) - (Kr_PH_CFP_IP3 .* IP3_PH_CFP_C));
	K_h_closed_total = ((Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* h_closed_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* h_open_init_molecules_um_2));
	h_closed = ((K_h_closed_total - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* h_open)) ./ (Size_M .* UnitFactor_uM_um3_molecules_neg_1));
	I_Leak_flux = (gL .* (Voltage_M - V_leak));
	C_K_flux = (gK_max .* (n_open ^ 4.0));
	I_K_flux = (C_K_flux .* ((((mlabfix_R_ .* mlabfix_T_) ./ (carrierValence_K_flux .* mlabfix_F_)) .* log((K_EX ./ K_C))) - Voltage_M));
	F_M = ( - (I_Leak_flux .* Size_M) - (I_Na_Flux .* Size_M) - (I_KCNQ2_3 .* Size_M) - (I_K_flux .* Size_M));
	J_GTPase_GaP = (15.0 .* GaGTP_PLC_M);
	Kf_G2beta = KfG2;
	KFlux_M_EX = (Size_M ./ Size_EX);
	J_PH_YFP_PIP2 = (((Kf_PH_YFP_PIP2 .* PH_YFP_C) .* PIP2_M) - (Kr_PH_YFP_PIP2 .* PH_YFP_PIP2_M));
	m_closed = ((K_m_closed_total - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* m_open)) ./ (Size_M .* UnitFactor_uM_um3_molecules_neg_1));
	J_Na_Gate_Act = ((alpha_m .* m_closed) - (beta_m .* m_open));
	J_NE_GaP = (4.7 .* GaGDP_PLC_M);
	J_PLCdiss = ((Kf_PLCdiss .* GaGDP_PLC_M) - ((Kr_PLCdiss .* PLC_M) .* Ga_GDP_M));
	Kf_G1beta = KfG1;
	J_Leak_flux = ((I_Leak_flux ./ (carrierValence_Leak_flux .* mlabfix_F_)) .* unitFactor_Leak_flux);
	J_NE_RLG = (Kf_NE_RLG .* RLG_GDP_M);
	J_PIP2hydr = (PIP2_M .* (PLC_basal + (K_plc .* GaGTP_PLC_M)));
	I_M = device_totalCurrClampElectrode.F;
	Kr_L2 = KrL2;
	Kf_L2 = KfL2;
	J_L2 = (((Kf_L2 .* oxoM_EX) .* RG_GDP_M) - (Kr_L2 .* RLG_GDP_M));
	J_PLC_on_PI4P = (PI4P_M .* PLC_efficiency_PIP .* (PLC_basal + (GaGTP_PLC_M .* K_plc)));
	Kf_L1 = KfL1;
	Kr_L1 = KrL1;
	J_L1 = (((Kf_L1 .* oxoM_EX) .* R_M) - (Kr_L1 .* RL_M));
	Kf_PH_CFP_PIP2 = speed_PH_PIP2;
	J_PH_CFP_PIP2 = (((Kf_PH_CFP_PIP2 .* PIP2_M) .* PH_CFP_C) - (Kr_PH_CFP_PIP2 .* PH_CFP_PIP2_M));
	J_reconstitution = (Kf_reconstitution .* Gbeta_M .* Ga_GDP_M);
	J_Na_Flux = ((I_Na_Flux ./ (carrierValence_Na_Flux .* mlabfix_F_)) .* unitFactor_Na_Flux);
	device_M.Capacitance = (C_M .* Size_M);
	allRL = (RL_M + RLG_GDP_M + RLGbeta_M);
	J_Na_Gate_Inact = ((alpha_h .* h_closed) - (beta_h .* h_open));
	Kr_G2 = KrG2;
	J_G2 = (((Kf_G2 .* RL_M) .* G_GDP_M) - (Kr_G2 .* RLG_GDP_M));
	Kr_G1 = KrG1;
	J_G1 = (((Kf_G1 .* G_GDP_M) .* R_M) - (Kr_G1 .* RG_GDP_M));
	K_n_closed_total = ((Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* n_closed_init_molecules_um_2) + (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* n_open_init_molecules_um_2));
	Ion_EX = Ion_EX_init_uM;
	allPH_cytoplasm = (PH_CFP_C + PH_YFP_C + IP3_PH_CFP_C + IP3_PH_YFP_C);
	J_G2beta = (((Kf_G2beta .* RL_M) .* Gbeta_M) - (Kr_G2beta .* RLGbeta_M));
	Ion_C = Ion_C_init_uM;
	n_closed = ((K_n_closed_total - (Size_M .* UnitFactor_uM_um3_molecules_neg_1 .* n_open)) ./ (Size_M .* UnitFactor_uM_um3_molecules_neg_1));
	J_K_Gate = ((alpha_n .* n_closed) - (beta_n .* n_open));
	allPIP2M = (PIP2_M + KCNQ_PIP2_M + PH_CFP_PIP2_M + PH_YFP_PIP2_M);
	J_G1beta = (((Kf_G1beta .* Gbeta_M) .* R_M) - (Kr_G1beta .* RGbeta_M));
	allGabg = (G_GDP_M + RG_GDP_M + RLG_GDP_M);
	J_DAGPase = (Kf_DAGPase .* DAG_M);
	J_GTPase_Ga = (Kf_GTPase_Ga .* GaGTP_M);
	J_K_flux = ((I_K_flux ./ (carrierValence_K_flux .* mlabfix_F_)) .* unitFactor_K_flux);
	allRGbeta = (RG_GDP_M + RLG_GDP_M + RLGbeta_M + RGbeta_M);
	% Rates
	dydt = [
		( - J_PIP2hydr - J_PH_YFP_PIP2 - J_PIP2bindKCNQ - J_PH_CFP_PIP2);    % rate for PIP2_M
		J_Na_Gate_Inact;    % rate for h_open
		(J_GTPase_Ga + J_PLCdiss - J_reconstitution);    % rate for Ga_GDP_M
		(J_L2beta + J_G2beta + J_NE_RLG);    % rate for RLGbeta_M
		( - (UnitFactor_uM_um3_molecules_neg_1 .* KFlux_M_C .* J_PH_CFP_PIP2) - J_PH_CFP_IP3);    % rate for PH_CFP_C
		J_K_Gate;    % rate for n_open
		(J_G2 - J_NE_RLG + J_L2);    % rate for RLG_GDP_M
		( - J_L1 - J_G1beta - J_G1);    % rate for R_M
		(J_DAGPase - J_PI4K_4Pase);    % rate for PI_M
		J_PH_CFP_PIP2;    % rate for PH_CFP_PIP2_M
		J_Na_Gate_Act;    % rate for m_open
		( - J_reconstitution + J_NE_G - J_G2beta - J_G1beta);    % rate for Gbeta_M
		(J_NE_RG + J_G1beta - J_L2beta);    % rate for RGbeta_M
		(J_PLCassoc - J_GTPase_GaP + J_NE_GaP);    % rate for GaGTP_PLC_M
		(J_PLC_on_PI4P - J_DAGPase + J_PIP2hydr);    % rate for DAG_M
		J_PH_YFP_IP3;    % rate for IP3_PH_YFP_C
		( - J_GTPase_Ga + J_NE_RG - J_PLCassoc + J_NE_RLG + J_NE_G);    % rate for GaGTP_M
		( - J_PH_CFP_IP3 + (UnitFactor_uM_um3_molecules_neg_1 .* KFlux_M_C .* J_PIP2hydr) - J_PH_YFP_IP3 - J_IP3Pase);    % rate for IP3_C
		((UnitFactor_mV_pF_s_neg_1_pA_neg_1 ./ device_M.Capacitance) .* (I_M - ( - (I_Leak_flux .* Size_M) - (I_Na_Flux .* Size_M) - (I_KCNQ2_3 .* Size_M) - (I_K_flux .* Size_M))));    % rate for Voltage_M
		( - J_G2 - J_G2beta + J_L1);    % rate for RL_M
		(J_PI4K_4Pase - J_PLC_on_PI4P);    % rate for PI4P_M
		 - J_PIP2bindKCNQ;    % rate for KCNQ_M
	];
end
