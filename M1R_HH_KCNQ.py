# Peter Duggins
# Cellular Mathematical Biology, AMATH 900, University of Waterloo
# M1R Signaling Cascade with KCNQ Channels in a Hodgkin-Huxley Neuron
# Dec 21, 2015

#Credit to http://www.gribblelab.org/compneuro/3_Modelling_Action_Potentials.html for much of the HH model

from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

def main():

        # set up a dictionary of parameters
        params = {

                'V_0' : -70e-3, # initial conditions
                't_max' : 3e0, 
                't_step' : 1e-2,
                'clamped_state' : 0, #False, 0=voltage, 3=pip2

                'oxoM_EX_0' : 0e1,
                'PIP2_M_0' : 5e3,    
                'KCNQ_M_0' : 0.0,     # initial state
                'KCNQ_PIP2_M_0' : 0.0,   
                'RLG_GDP_M_0' : 0.0,   
                'RG_GDP_M_0' : 0.0,   
                'DAG_M_0' : 2.0e3,
                'Ga_GTP_M_0' : 0.0,   
                'G_GDP_M_0' : 40.0,  
                'RG_beta_M_0' : 0.0,   
                'PI4P_M_0' : 4.0e3,
                'R_M_0' : 1.0,   
                'PI_M_0' : 1.4e5,      
                'Ga_GDP_M_0' : 0.0,   
                'Ga_GTP_PLC_M_0' : 0.0,   
                'RL_M_0' : 0.0,   
                'Ga_GDP_PLC_M_0' : 0.0,   
                'RLG_beta_M_0' : 0.0,    
                'G_beta_M_0' : 0.0,    
                'PLC_M_0' : 10.0,
                'IP3_C_0' : 0.16,
                'h_0' : 0.57,
                'm_0' : 0.058,
                'n_0' : 0.33,  

                'I_ext'  : 1.0e-9,   # injected current
                'I_start' : 0,
                'I_end' : 0.03,

                'leak_E' : -7.0e-2,     # channel paramters
                'leak_G' : 3.0e-09,
                'Na_E'          : 5.0e-2,
                'Na_G'          : 1.0e-6,
                'k_E'           : -9.0e-2,
                'k_G'           : 2.0/3.5*1e-8, #2.0e-7
                'KCNQ_E' : -9.0e-2,
                'KCNQ_G': 2.0e-8,
                'k_Na_act'      : 3.0e+0,
                'A_alpha_m' : 2.0e+5,
                'B_alpha_m' : -4.0e-2,
                'C_alpha_m' : 1.0e-3,
                'A_beta_m'  : 6.0e+4,
                'B_beta_m'  : -4.9e-2,
                'C_beta_m'  : 2.0e-2,
                'l_Na_'    : 1.0e+0,
                'A_alpha_h' : 8.0e+4,
                'B_alpha_h' : -4.0e-2,
                'C_alpha_h' : 1.0e-3,
                'A_beta_h'  : 4.0e+2,
                'B_beta_h'  : -3.6e-2,
                'C_beta_h'  : 2.0e-3,
                'k_K'           : 4.0e+0,
                'A_alpha_n' : 2.0e+4,
                'B_alpha_n' : -3.1e-2,
                'C_alpha_n' : 8.0e-4,
                'A_beta_n'  : 5.0e+3,
                'B_beta_n'  : -2.8e-2,
                'C_beta_n'  : 4.0e-4,
                'z' : 1.866,
                'F' : 96480,
                'R' : 8314,
                'T' : 300,
                'theta' : 10000,
                'kg' : 0.02,
                'kv0' : 4.512,
                'kp0' : 0.75,
                'C_M'           : 3.0e-11,

                'Kr_PLCassoc' : 0.0,         #intracellular kinetic paramters
                'Kf_NE_RLG' : 0.65,  
                'Kr_PI4K_4Pase' : 0.006, 
                'K_plc' : 0.1,   
                'Kf_GTPase_Ga' : 0.026, 
                'KrG2' : 0.68,  
                'Kf_reconstitution' : 1.0,   
                'Kf_PLCdiss' : 0.71,  
                'Kf_PI4K_4Pase' : 2.6E-4,
                'KfG2' : 0.0026666,     
                'Kf_NE_G' : 1.5E-5,
                'PLC_basal' : 0.0025,
                'Kf_NE_RG' : 1.5E-5,
                'speed_PIP2_KCNQ' : 0.05,  
                'KD_PH_PIP2' : 2.0,   
                'K_IP3ase' : 0.08,  
                'PLC_efficiency_PIP' : 0.14,  
                'KD_PH_IP3' : 0.1,    
                'Kr_NE_G' : 0.0,   
                'Kf_PLCassoc' : 1.0,   
                'Hill_binding' : 1.0,   
                'Kf_DAGPase' : 0.02,  
                'speed_PH_PIP2' : 1.0,    
                'alpha' : 100.0,  
                'Kr_PLCdiss' : 0.0,   
                'KA_PIP2_KCNQ' : 2000.0, 
                'Kr_NE_RG' : 0.0,    
                'KL1' : 2.0,    
                'KrL1' : 5.555555555555555,      
                'speed_PH_IP3' : 10.0,
}

        # set initial states and time vector
        state0 = [
                params['V_0'],
                params['KCNQ_M_0'],
                params['KCNQ_PIP2_M_0'],
                params['PIP2_M_0'],
                params['RLG_GDP_M_0'],
                params['RG_GDP_M_0'],
                params['RG_beta_M_0'],
                params['RLG_beta_M_0'],
                params['R_M_0'],
                params['RL_M_0'],
                params['Ga_GTP_M_0'],
                params['G_GDP_M_0'],
                params['Ga_GDP_M_0'],
                params['Ga_GTP_PLC_M_0'],
                params['Ga_GDP_PLC_M_0'],
                params['G_beta_M_0'],
                params['oxoM_EX_0'],
                params['PLC_M_0'],
                params['IP3_C_0'],
                params['DAG_M_0'],
                params['PI4P_M_0'],
                params['PI_M_0'],
                params['m_0'],
                params['h_0'],
                params['n_0']
        ]

        t = np.arange(0, params['t_max'],params['t_step'])

        # run simulation
        # state = odeint(neuron, state0, t, args=(params,), False, False)
        # voltage_values=[-120*1e-3,-110*1e-3,-100*1e-3,-90*1e-3,-80*1e-3,-70*1e-3,-60*1e-3,
                # -50*1e-3,-40*1e-3,-30*1e-3,-20*1e-3,-10*1e-3,0*1e-3]
        voltage_values=np.linspace(-120e-3,80e-3,20)
        print voltage_values
        conductance_voltage_experiment(state0,t,params,voltage_values)

        # plot the results
        # fig=plt.figure(figsize=(8,12))
        # plt.subplot(2,2,1)
        # plt.plot(t, state[:,0])
        # plt.title('membrane potential V (mV)')
        # plt.subplot(2,2,2)
        # plt.plot(t, state[:,-3])
        # plt.title('Na activation (m) and inactivation (h)')
        # plt.subplot(2,2,2)
        # plt.plot(t, state[:,-2])
        # #plt.title('Na inactivation (h)')
        # plt.subplot(2,2,3)
        # plt.plot(t, state[:,-1])
        # plt.title('K channel activation (n)')
        # #quick recalculation of KCNQ_open
        # KCNQ_open_recalc=[]
        # kv0=params['kv0']
        # z=params['z']
        # F=params['F']
        # R=params['R']
        # T=params['T']
        # kp0=params['kp0']
        # kg=params['kg']
        # theta=params['theta']
        # for i in range(len(t)):
        #         kv=kv0*np.exp(z*F*(state[:,0][i]*1e3)/(R*T))
        #         kp=kp0*state[:,2][i]
        #         PP0=(kg+kv*kg+kp*kg+theta*kv*kp*kg)/(1+kg+kp+kv+kv*kg+kp*kg+kv*kp+theta*kv*kp*kg)
        #         PP0_max=(kg+theta*kp*kg)/(1+kg+kp+theta*kp*kg)
        #         KCNQ_open_recalc.append(PP0/PP0_max)
        # plt.subplot(2,2,4)
        # plt.plot(t, KCNQ_open_recalc)
        # plt.title('KCNQ channel activation')
        # plt.show()
        # fig.savefig('states0,22-24,kcqn,oxom=0')

def neuron(state, t, params):

        # State lookup
        V = state[0]
        KCNQ_M = state[1]
        KCNQ_PIP2_M = state[2]
        PIP2_M = state[3]
        RLG_GDP_M = state[4]
        RG_GDP_M = state[5]
        RG_beta_M = state[6]
        RLG_beta_M = state[7]
        R_M = state[8]
        RL_M = state[9]
        Ga_GTP_M = state[10]
        G_GDP_M = state[11]
        Ga_GDP_M = state[12]
        Ga_GTP_PLC_M = state[13]
        Ga_GDP_PLC_M = state[14]
        G_beta_M = state[15]
        oxoM_EX = state[16]
        PLC_M = state[17]
        IP3_C = state[18]
        DAG_M = state[19]
        PI4P_M = state[20]
        PI_M = state[21]
        m = state[22]
        h = state[23]
        n = state[24]

        # parameter lookup
        KL1 = params['KL1']
        KrL1 = params['KrL1']
        Kr_PLCassoc = params['Kr_PLCassoc']
        Kf_NE_RLG = params['Kf_NE_RLG']
        Kr_PI4K_4Pase = params['Kr_PI4K_4Pase']
        K_plc = params['K_plc']
        Kf_GTPase_Ga = params['Kf_GTPase_Ga']
        KrG2 = params['KrG2']
        Kf_reconstitution = params['Kf_reconstitution']
        Kf_PLCdiss = params['Kf_PLCdiss']
        Kf_PI4K_4Pase = params['Kf_PI4K_4Pase']
        KfG2 = params['KfG2']
        Kf_NE_G = params['Kf_NE_G']
        Kf_NE_RG = params['Kf_NE_RG']
        KD_PH_PIP2 = params['KD_PH_PIP2']
        K_IP3ase = params['K_IP3ase']
        KD_PH_IP3 = params['KD_PH_IP3']
        Kr_NE_G = params['Kr_NE_G']
        Kf_PLCassoc = params['Kf_PLCassoc']
        Kf_DAGPase = params['Kf_DAGPase']
        Kr_PLCdiss = params['Kr_PLCdiss']
        KA_PIP2_KCNQ = params['KA_PIP2_KCNQ']
        Kr_NE_RG = params['Kr_NE_RG']
        PLC_basal = params['PLC_basal']
        PLC_efficiency_PIP = params['PLC_efficiency_PIP']
        Hill_binding = params['Hill_binding']
        alpha = params['alpha']
        theta = params['theta']
        speed_PH_PIP2 = params['speed_PH_PIP2']
        speed_PIP2_KCNQ = params['speed_PIP2_KCNQ']
        speed_PH_IP3 = params['speed_PH_IP3']
        C_M = params['C_M']
        kv0=params['kv0']
        z=params['z']
        F=params['F']
        R=params['R']
        T=params['T']
        kp0=params['kp0']
        kg=params['kg']
        theta=params['theta']

        # Intracellular species kinetics
        K_RG_GDP_M_total = RLG_beta_M + RL_M + R_M + RG_beta_M + RG_GDP_M + RLG_GDP_M
        RG_GDP_M = - RLG_beta_M - RL_M - R_M - RG_beta_M + K_RG_GDP_M_total - RLG_GDP_M
        K_G_GDP_M_total = G_beta_M - RL_M + G_GDP_M - R_M
        Kf_G2 = KfG2
        KfG1 = 0.1 * KfG2
        Kf_G1 = KfG1
        K_PLC_M_total = G_beta_M + RLG_beta_M - Ga_GTP_M + RG_beta_M + PLC_M - Ga_GDP_M
        PLC_M = - G_beta_M - RLG_beta_M + Ga_GTP_M - RG_beta_M + K_PLC_M_total + Ga_GDP_M
        KG2 = (KrG2 / KfG2)
        KG1 = (KG2 * alpha)
        # K_oxoM_EX_total = RLG_beta_M + RL_M + oxoM_EX + RLG_GDP_M
        # oxoM_EX = - RLG_beta_M - RL_M + K_oxoM_EX_total - RLG_GDP_M
        KfL1 = KrL1 / KL1
        KfL2 = KfL1
        Kf_L2beta = KfL2
        KL2 = KL1 / alpha
        KrL2 = KL2 * KfL2
        Kr_L2beta = KrL2
        Kr_G2beta = KrG2
        KrG1 = KfG1 * KG1
        K_KCNQ_PIP2_M_total = KCNQ_PIP2_M + KCNQ_M
        KCNQ_PIP2_M = K_KCNQ_PIP2_M_total - KCNQ_M
        K_Ga_GDP_PLC_M_total = Ga_GTP_PLC_M - G_beta_M - RLG_beta_M + Ga_GTP_M - RG_beta_M + Ga_GDP_PLC_M + Ga_GDP_M
        Kr_G1beta = KrG1
        G_GDP_M =  -G_beta_M + RL_M + K_G_GDP_M_total + R_M
        Kf_PIP2bindKCNQ = speed_PIP2_KCNQ
        KD_PIP2_KCNQ = KA_PIP2_KCNQ ** Hill_binding
        Kr_PIP2bindKCNQ = KD_PIP2_KCNQ * Kf_PIP2bindKCNQ
        Ga_GDP_PLC_M = - Ga_GTP_PLC_M + G_beta_M + RLG_beta_M - Ga_GTP_M + RG_beta_M + K_Ga_GDP_PLC_M_total - Ga_GDP_M
        allGa_PLC = Ga_GTP_PLC_M + Ga_GDP_PLC_M
        Kf_G2beta = KfG2
        Kf_G1beta = KfG1
        Kr_L2 = KrL2
        Kf_L2 = KfL2
        Kf_L1 = KfL1
        Kr_L1 = KrL1
        allRL = RL_M + RLG_GDP_M + RLG_beta_M
        Kr_G2 = KrG2
        Kr_G1 = KrG1
        allPIP2M = PIP2_M + KCNQ_PIP2_M
        allGabg = G_GDP_M + RG_GDP_M + RLG_GDP_M
        allRG_beta = RG_GDP_M + RLG_GDP_M + RLG_beta_M + RG_beta_M

        J_PI4K_4Pase = (Kf_PI4K_4Pase * PI_M) - (Kr_PI4K_4Pase * PI4P_M)
        J_NE_RG = (Kf_NE_RG * RG_GDP_M - (Kr_NE_RG * Ga_GTP_M) * RG_beta_M)
        J_IP3Pase = K_IP3ase * IP3_C
        J_PLCassoc = (((Kf_PLCassoc * Ga_GTP_M) * PLC_M) - (Kr_PLCassoc * Ga_GTP_PLC_M))
        J_L2beta = (((Kf_L2beta * RG_beta_M) * oxoM_EX) - (Kr_L2beta * RLG_beta_M))
        J_NE_G = (Kf_NE_G * G_GDP_M) - (Kr_NE_G * G_beta_M * Ga_GTP_M)
        J_PIP2bindKCNQ = Kf_PIP2bindKCNQ * (PIP2_M ** Hill_binding) * KCNQ_M - (Kr_PIP2bindKCNQ * KCNQ_PIP2_M)
        J_NE_GaP = 4.7 * Ga_GDP_PLC_M
        J_GTPase_GaP = 15.0 * Ga_GTP_PLC_M
        J_L2 = Kf_L2 * oxoM_EX * RG_GDP_M - Kr_L2 * RLG_GDP_M
        J_PLC_on_PI4P = PI4P_M * PLC_efficiency_PIP * (PLC_basal + Ga_GTP_PLC_M * K_plc)
        J_PLCdiss = Kf_PLCdiss * Ga_GDP_PLC_M - Kr_PLCdiss * PLC_M * Ga_GDP_M
        J_NE_RLG = Kf_NE_RLG * RLG_GDP_M
        J_PIP2hydr = PIP2_M * (PLC_basal + K_plc * Ga_GTP_PLC_M)
        J_L1 = Kf_L1 * oxoM_EX * R_M - Kr_L1 * RL_M
        J_reconstitution = Kf_reconstitution * G_beta_M * Ga_GDP_M
        J_G2 = Kf_G2 * RL_M * G_GDP_M - Kr_G2 * RLG_GDP_M
        J_G1 = Kf_G1 * G_GDP_M * R_M - Kr_G1 * RG_GDP_M
        J_G2beta = Kf_G2beta * RL_M * G_beta_M - Kr_G2beta * RLG_beta_M
        J_G1beta = Kf_G1beta * G_beta_M * R_M - Kr_G1beta * RG_beta_M
        J_DAGPase = Kf_DAGPase * DAG_M
        J_GTPase_Ga = Kf_GTPase_Ga * Ga_GTP_M

        # Na Current
        alpha_act = params['A_alpha_m'] * (V-params['B_alpha_m']) / (1.0 - np.exp((params['B_alpha_m']-V) / params['C_alpha_m']))
        beta_act = params['A_beta_m'] * (params['B_beta_m']-V) / (1.0 - np.exp((V-params['B_beta_m']) / params['C_beta_m']) )
        alpha_inact = params['A_alpha_h'] * (params['B_alpha_h']-V) / (1.0 - np.exp((V-params['B_alpha_h']) / params['C_alpha_m']))
        beta_inact  = params['A_beta_h'] / (1.0 + (np.exp((params['B_beta_h']-V) / params['C_beta_h'])))
        I_Na =(params['Na_E']-V) * params['Na_G'] * (m**params['k_Na_act']) * h

        # K Current
        alpha_kal = params['A_alpha_n'] * (V-params['B_alpha_n']) / (1.0 - np.exp((params['B_alpha_n']-V) / params['C_alpha_n']))
        beta_kal = params['A_beta_n'] * (params['B_beta_n']-V) / (1.0 - np.exp((V-params['B_beta_n']) / params['C_beta_n']))
        I_K = (params['k_E']-V) * params['k_G'] * n**params['k_K']

        # Leak current
        I_leak = (params['leak_E']-V) * params['leak_G']

        # KCNQ current, activation kinetics assumed to be instantaneous
        kv=kv0*np.exp(z*F*(V*1e3)/(R*T))
        kp = kp0 * KCNQ_PIP2_M
        PP0=((kg+kv*kg+kp*kg+theta*kv*kp*kg)/(1+kg+kp+kv+kv*kg+kp*kg+kv*kp+theta*kv*kp*kg))
        PP0_max=((kg+theta*kp*kg)/(1+kg+kp+theta*kp*kg))
        KCNQ_open=PP0/PP0_max
        I_KCNQ = (params['KCNQ_E'] - V) * params['KCNQ_G'] * KCNQ_open

        # External current
        if t > params['I_start'] and t < params['I_end']:
                I_ext = params['I_ext']
        else:
                I_ext = 0

        # calculate derivatives
        dVdt = (I_leak + I_K + I_Na + I_ext + I_KCNQ) / C_M
        dKCNQ_Mdt = - J_PIP2bindKCNQ
        dKCNQ_PIP2_Mdt = J_PIP2bindKCNQ
        dPIP2_Mdt = - J_PIP2hydr - J_PIP2bindKCNQ
        dRLG_GDP_M = J_G2 - J_NE_RLG + J_L2
        dRG_GDP_Mdt = J_G1 - J_NE_RG - J_L2
        dRG_beta_Mdt = J_NE_RG + J_G1beta - J_L2beta
        dRLG_beta_Mdt = J_L2beta + J_G2beta + J_NE_RLG
        dR_Mdt = - J_L1 - J_G1beta - J_G1
        dRL_Mdt = - J_G2 - J_G2beta + J_L1
        dGa_GTP_Mdt = - J_GTPase_Ga + J_NE_RG - J_PLCassoc + J_NE_RLG + J_NE_G
        dG_GDP_Mdt = J_reconstitution - J_G2 - J_G1 - J_NE_G
        dGa_GDP_Mdt = J_GTPase_Ga + J_PLCdiss - J_reconstitution
        dGa_GTP_PLC_Mdt = J_PLCassoc - J_GTPase_GaP + J_NE_GaP
        dGa_GDP_PLC_Mdt = J_GTPase_GaP - J_NE_GaP - J_PLCdiss
        dG_Beta_Mdt = -J_reconstitution + J_NE_G - J_G2beta - J_G1beta
        dOxo_M_Exdt = 0#- J_L1 - J_L2beta - J_L2
        dPLC_Mdt = J_PLCdiss - J_PLCassoc
        dIP3_Cdt = J_PIP2hydr - J_IP3Pase
        dDAG_Mdt = J_PLC_on_PI4P - J_DAGPase + J_PIP2hydr
        dPI4P_Mdt = J_PI4K_4Pase - J_PLC_on_PI4P
        dPI_Mdt = J_DAGPase - J_PI4K_4Pase
        dmdt = ( alpha_act * (1.0 - m) ) - ( beta_act * m )
        dhdt = ( alpha_inact*(1.0 - h) ) - ( beta_inact*h )
        dndt = ( alpha_kal*(1.0 - n) ) - ( beta_kal*n )

        # new state
        state_new = [
                dVdt,
                dKCNQ_Mdt,
                dKCNQ_PIP2_Mdt,
                dPIP2_Mdt,
                dRLG_GDP_M,
                dRG_GDP_Mdt,
                dRG_beta_Mdt,
                dRLG_beta_Mdt,
                dR_Mdt,
                dRL_Mdt,
                dGa_GTP_Mdt,
                dG_GDP_Mdt,
                dGa_GDP_Mdt,
                dGa_GTP_PLC_Mdt,
                dGa_GDP_PLC_Mdt,
                dG_Beta_Mdt,
                dOxo_M_Exdt,
                dPLC_Mdt,
                dIP3_Cdt,
                dDAG_Mdt,
                dPI4P_Mdt,
                dPI_Mdt,
                dmdt,
                dhdt,
                dndt
        ]
        
        #set the derivative of the state numbered "clamped_state" to zero
        if params['clamped_state'] != 9000:
                # print state_new
                state_new.insert(params['clamped_state'],0) 
                state_new.pop(params['clamped_state'] + 1)
                # print state_new

        return state_new


def conductance_voltage_experiment(state0,t,params,value_list):

        kv0=params['kv0']
        z=params['z']
        F=params['F']
        R=params['R']
        T=params['T']
        kp0=params['kp0']
        kg=params['kg']
        theta=params['theta']
        y_values=[]

        for i in range(len(value_list)):
                params['V_0'] = value_list[i]
                state0 = [
                        params['V_0'],
                        params['KCNQ_M_0'],
                        params['KCNQ_PIP2_M_0'],
                        params['PIP2_M_0'],
                        params['RLG_GDP_M_0'],
                        params['RG_GDP_M_0'],
                        params['RG_beta_M_0'],
                        params['RLG_beta_M_0'],
                        params['R_M_0'],
                        params['RL_M_0'],
                        params['Ga_GTP_M_0'],
                        params['G_GDP_M_0'],
                        params['Ga_GDP_M_0'],
                        params['Ga_GTP_PLC_M_0'],
                        params['Ga_GDP_PLC_M_0'],
                        params['G_beta_M_0'],
                        params['oxoM_EX_0'],
                        params['PLC_M_0'],
                        params['IP3_C_0'],
                        params['DAG_M_0'],
                        params['PI4P_M_0'],
                        params['PI_M_0'],
                        params['m_0'],
                        params['h_0'],
                        params['n_0']
                ]
                state = odeint(neuron, state0, t, args=(params,))
                kv=kv0*np.exp(z*F*(state[:,0][-1]*1e3)/(R*T))
                kp=kp0*state[:,2][-1]
                PP0=((kg+kv*kg+kp*kg+theta*kv*kp*kg)/(1+kg+kp+kv+kv*kg+kp*kg+kv*kp+theta*kv*kp*kg))
                PP0_max=((kg+theta*kp*kg)/(1+kg+kp+theta*kp*kg))
                KCNQ_open_recalc = PP0/PP0_max
                y_values.append(KCNQ_open_recalc) # measure the final conductance value, divided by its maximum value 
                # print kv
                # print kp
                # print KCNQ_open_recalc

        fig=plt.figure(figsize=(8,12))
        plt.subplot(1,1,1)
        plt.plot(value_list,y_values)
        plt.xlabel('Membrane Voltage (V)')
        plt.ylabel('G/G_max of the KCNQ Channel')
        plt.show()
        # plt.title('K channel activation (n)')
        # fig.savefig('Conductance-Voltage Relationship for KCNQ2/3 Channel.png')

def conductance_pip2_experiment(state0,t,params,value_list):

        kv0=params['kv0']
        z=params['z']
        F=params['F']
        R=params['R']
        T=params['T']
        kp0=params['kp0']
        kg=params['kg']
        theta=params['theta']
        y_values=[]

        for i in range(len(value_list)):
                params['PIP2_M_0'] = value_list[i]
                state0 = [
                        params['V_0'],
                        params['KCNQ_M_0'],
                        params['KCNQ_PIP2_M_0'],
                        params['PIP2_M_0'],
                        params['RLG_GDP_M_0'],
                        params['RG_GDP_M_0'],
                        params['RG_beta_M_0'],
                        params['RLG_beta_M_0'],
                        params['R_M_0'],
                        params['RL_M_0'],
                        params['Ga_GTP_M_0'],
                        params['G_GDP_M_0'],
                        params['Ga_GDP_M_0'],
                        params['Ga_GTP_PLC_M_0'],
                        params['Ga_GDP_PLC_M_0'],
                        params['G_beta_M_0'],
                        params['oxoM_EX_0'],
                        params['PLC_M_0'],
                        params['IP3_C_0'],
                        params['DAG_M_0'],
                        params['PI4P_M_0'],
                        params['PI_M_0'],
                        params['m_0'],
                        params['h_0'],
                        params['n_0']
                ]
                state = odeint(neuron, state0, t, args=(params,))
                kv=kv0*np.exp(z*F*(state[:,0][-1]*1e3)/(R*T))
                kp=kp0*state[:,2][-1]
                PP0=(kg+kv*kg+kp*kg+theta*kv*kp*kg)/(1+kg+kp+kv+kv*kg+kp*kg+kv*kp+theta*kv*kp*kg)
                PP0_max=(kg+theta*kp*kg)/(1+kg+kp+theta*kp*kg)
                KCNQ_open_recalc = PP0/PP0_max
                y_values.append(KCNQ_open_recalc) # measure the final conductance value, divided by its maximum value 
                # print kv
                # print kp
                # print KCNQ_open_recalc

        fig=plt.figure(figsize=(8,12))
        plt.subplot(1,1,1)
        plt.plot(value_list,y_values)
        plt.xlabel('Membrane Voltage (V)')
        plt.ylabel('G/G_max of the KCNQ Channel')
        plt.show()
        # plt.title('K channel activation (n)')
        # fig.savefig('Conductance-Voltage Relationship for KCNQ2/3 Channel.png')




main()