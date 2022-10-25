"################### Synthetic data from case study ###################"

import numpy as np
from scipy.integrate import odeint

np.random.seed(47625)

# function that returns dz/dt
def kinetic_model(z,t):
    k_1 = 100 #k1+ * k2+ * k3+
    k_2 = 3 #k1- * k2- * k3-
    k_3 = 4 #k1+ * (k2+ + k3+ + k2-)
    dTdt = (-1)*((k_1*z[0]*z[1])/((1+k_2*z[0]+((k_3**0.5)*(z[1]**0.5)))**3))
    dHdt = (-3)*((k_1*z[0]*z[1])/((1+k_2*z[0]+((k_3**0.5)*(z[1]**0.5)))**3))
    dMdt = ((k_1*z[0]*z[1])/((1+k_2*z[0]+((k_3**0.5)*(z[1]**0.5)))**3))
    dzdt = [dTdt,dHdt,dMdt]
    return dzdt

# time points
t = np.linspace(0,50,10)
number_exp = 3

# CT_beg = 10
# CT_end = 12
# CH_beg  = 12
# CH_end = 19
# CM_beg = 2
# CM_end = 6
# initial_conditions_T = np.linspace(CT_beg,CT_end,number_exp)
# initial_conditions_H = np.linspace(CH_beg,CH_end,number_exp)
# initial_conditions_M = np.linspace(CM_beg,CM_end,number_exp)
# initial_conditions = np.vstack([initial_conditions_T,initial_conditions_H,initial_conditions_M])
# for i in range(initial_conditions.shape[0]):
#     np.random.shuffle(initial_conditions[i])

initial_conditions = np.array([[10,10,12],[19,19,19],[6,2,4]])

# Get the rate
k_1 = 100 #k1+ * k2+ * k3+
k_2 = 3 #k1- * k2- * k3-
k_3 = 4 #k1+ * (k2+ + k3+ + k2-)
rate = np.empty([number_exp,len(t)])
z = np.empty([number_exp,len(t),3])
noisy_data = np.empty([number_exp,len(t),3])
noise_T = np.empty([number_exp,len(t)])
noise_H = np.empty([number_exp,len(t)])
noise_M = np.empty([number_exp,len(t)])
STD_T = 0.46878853
STD_H = 0.5313656
STD_M = 0.23121147

# different experiments at different initial conditions for A (maybe change B also?)
for i in range(initial_conditions.shape[1]):

    # initial condition
    # step_size = (CA_end - CA_beg) / (number_exp - 1)
    # res = CA_beg / step_size
    z0 = initial_conditions[:,i]
    # p = int(j/step_size - res + 1e-5)
    noise_T[i] = np.random.normal(0,STD_T,len(t))
    noise_H[i] = np.random.normal(0,STD_H,len(t))
    noise_M[i] = np.random.normal(0,STD_M,len(t))
    # noise_T[i] = np.zeros(len(t))
    # noise_H[i] = np.zeros(len(t))
    # noise_M[i] = np.zeros(len(t))

    # solve ODE
    z[i] = odeint(kinetic_model,z0,t)
    noisy_data[i][:,0] = z[i][:,0] + noise_T[i]
    noisy_data[i][:,1] = z[i][:,1] + noise_H[i]
    noisy_data[i][:,2] = z[i][:,2] + noise_M[i]

for a in range(number_exp):
    for b in range(len(t)):
        for c in range(3):
            if noisy_data[a][b][c] < 0:
                noisy_data[a][b][c] = 0

for i in range(initial_conditions.shape[1]):
    np.savetxt(f"data_{i}.csv", noisy_data[i], delimiter=",")
