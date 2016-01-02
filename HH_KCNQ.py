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

                'I_ext'  : 0,
                'I_start' : 0,
                'I_end' : 0,
                't_max' : 0.2,
                't_step' : 0.001,

                'leak_E' : -7.0e-2,
                'leak_G' : 3.0e-09,
                'Na_E'          : 5.0e-2,
                'Na_G'          : 1.0e-6,
                'k_E'           : -9.0e-2,
                'k_G'           : 2.0e-7,

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
                'C_m'           : 3.0e-11,
                'k_K'           : 4.0e+0,
                'A_alpha_n' : 2.0e+4,
                'B_alpha_n' : -3.1e-2,
                'C_alpha_n' : 8.0e-4,
                'A_beta_n'  : 5.0e+3,
                'B_beta_n'  : -2.8e-2,
                'C_beta_n'  : 4.0e-4
        }

        # set initial states and time vector
        state0 = [-70e-03, 0, 1, 0]
        t = np.arange(0, params['t_max'],params['t_step'])

        # external current
        params['I_ext'] = 1.0e-9
        params['I_start'] = 0.0
        params['I_end'] = 0.03

        # run simulation
        state = odeint(neuron, state0, t, args=(params,))

        # plot the results
        plt.figure(figsize=(8,12))
        plt.subplot(4,1,1)
        plt.plot(t, state[:,0])
        plt.title('membrane potential V (mV)')
        plt.subplot(4,1,2)
        plt.plot(t, state[:,1])
        plt.title('Na activation (m)')
        plt.subplot(4,1,3)
        plt.plot(t, state[:,2])
        plt.title('Na inactivation (h)')
        plt.subplot(4,1,4)
        plt.plot(t, state[:,3])
        plt.title('K channel activation (n)')
        plt.xlabel('time (sec)')
        plt.show()

def neuron(state, t, params):
        V = state[0]
        m = state[1]
        h = state[2]
        n = state[3]

        # Na Current
        alpha_act = params['A_alpha_m'] * (V-params['B_alpha_m']) / (1.0 - np.exp((params['B_alpha_m']-V) / params['C_alpha_m']))
        beta_act = params['A_beta_m'] * (params['B_beta_m']-V) / (1.0 - np.exp((V-params['B_beta_m']) / params['C_beta_m']) )
        alpha_inact = params['A_alpha_h'] * (params['B_alpha_h']-V) / (1.0 - np.exp((V-params['B_alpha_h']) / params['C_alpha_m']))
        beta_inact  = params['A_beta_h'] / (1.0 + (np.exp((params['B_beta_h']-V) / params['C_beta_h'])))
        dmdt = ( alpha_act * (1.0 - m) ) - ( beta_act * m )
        dhdt = ( alpha_inact*(1.0 - h) ) - ( beta_inact*h )
        I_Na =(params['Na_E']-V) * params['Na_G'] * (m**params['k_Na_act']) * h

        # K Current
        alpha_kal = params['A_alpha_n'] * (V-params['B_alpha_n']) / (1.0 - np.exp((params['B_alpha_n']-V) / params['C_alpha_n']))
        beta_kal = params['A_beta_n'] * (params['B_beta_n']-V) / (1.0 - np.exp((V-params['B_beta_n']) / params['C_beta_n']))
        dndt = ( alpha_kal*(1.0 - n) ) - ( beta_kal*n )
        I_K = (params['k_E']-V) * params['k_G'] * n**params['k_K']

        # Leak current
        I_leak = (params['leak_E']-V) * params['leak_G']

        # External current
        if t > params['I_start'] and t < params['I_end']:
                I_ext = params['I_ext']
        else:
                I_ext = 0

        # calculate derivative of E
        dVdt = (I_leak + I_K + I_Na + I_ext) / params['C_m']
        state_new = [dVdt, dmdt, dhdt, dndt]

        return state_new

main()
