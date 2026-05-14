import numpy as np
import matplotlib.pyplot as plt

#  Loop de Integração (Heun/Predictor-Corrector)               
def calculate_derivatives_cmlos(r_t, v_t, a_t, r_m, v_m, beta_t, x_filter_state):
    theta_t = np.arctan2(r_t[1], r_t[0])
    theta_m = np.arctan2(r_m[1], r_m[0])
    mod_rm = np.linalg.norm(r_m)
    unit_rm = r_m / mod_rm
    mod_rt = np.linalg.norm(r_t)
    unit_rt = r_t / mod_rt

    v_t = VT * np.array([-np.cos(beta_t), np.sin(beta_t)])
    
    
    # Erro de entrada do filtro (y no diagrama)
    error_linear = mod_rm * (theta_t - theta_m) 
    tau_lag = 1.0 / 20.0
    tau_lead = 1.0 / 2.0
    gain = 10.0
    
    # Derivada do estado do filtro
    dx_filter = (error_linear - x_filter_state) / tau_lag
    
    # Saída do filtro (n_c) considerando o ganho e o avanço (lead)
    # n_c = gain * (error_linear + tau_lead * dx_filter)
    Gs = gain * (x_filter_state + tau_lead * dx_filter)

    dot_rt = np.dot(v_t,unit_rt)
    dot_rm = np.dot(v_m, unit_rm) 

    dot_theta_t = (r_t[0]*v_t[1] - r_t[1]*v_t[0]) / (mod_rt**2)
    ddot_theta_t = (a_t[1]*np.cos(theta_t) - a_t[0]*np.sin(theta_t) - 2*dot_rt*dot_theta_t)/mod_rt

    x_nc = Gs + mod_rm*ddot_theta_t + 2*dot_rm*dot_theta_t

    r_tm = r_t - r_m
    x_lam = np.arctan2(r_tm[1], r_tm[0])
    a_m = x_nc * np.array([-np.sin(x_lam), np.cos(x_lam)])
    betatd = XNT / VT

    return v_t, v_m, a_m, betatd, x_nc, dx_filter


# Inicialização de constantes e variáveis (Baseado em Screenshot 2026-04-25 201010.png)
VM = 3000. # v missile 
VT = 1000. # v target 
XNT = 0.
RM1IC, RM2IC = 0., 1. # initial pos missile
RT1IC, RT2IC = 40000., 10000. # initial pos target
HEDEG = 0. # error angle deg
HE = HEDEG / 57.3 # error angle rad
XNP = 10. # proportional increment
TS = 0.1 # timestep

r_m = np.array([RM1IC, RM2IC])
r_t = np.array([RT1IC, RT2IC])
beta_t = 0. # target orientation
v_t = VT*np.array([-np.cos(beta_t), np.sin(beta_t)]) # 
T = 0.
S = 0.

# Primeiros calculos
r_tm = r_t - r_m
mod_rtm = np.linalg.norm(r_tm)
theta_t = np.arctan2(r_t[1], r_t[0])
v_m = VM*np.array([np.cos(theta_t+HE), np.sin(theta_t+HE)]) # TODO entender
v_tm = v_t - v_m
v_c = -(r_tm[0]*v_tm[0] + r_tm[1]*v_tm[1])/mod_rtm # TODO entender
a_t = np.zeros_like([2,1])

# Inicializar filtro
x_filter = 0.0

rt_history = []
rm_history = []
xnc_history = []
time_history = []


# Abre arquivo para escrita
with open('DATFIL.txt', 'w') as f: # TODO verificar se é necessario
    while v_c >= 0: # Substitui o IF(VC<0.) GOTO 999
        # Determinação do passo H
        if mod_rtm < 1000:
            H = 0.0002 
        else: 
            H = 0.01
        
        # Salvando estados antigos para integração numérica
        beta_told = beta_t         # WARNING: Arrays are referencing objects, so r_told = r_t changes the value of r_told later as we change r_t
        r_told = r_t.copy()        # using .copy() makes me not create a referencing object, so r_told will not change after operations with r_t
        r_mold = r_m.copy()
        v_mold = v_m.copy()
        
        v_tp, v_mp, a_mp, beta_p, xnc_p, dxf_p = calculate_derivatives_cmlos(r_t, v_t, a_t, r_m, v_m, beta_t, x_filter)

        # Salve o estado antigo para o Corretor
        x_filter_old = x_filter

        # Atualize os estados (Preditor)
        x_filter += H * dxf_p

        beta_t += H * beta_p
        r_t += H * v_tp
        r_m += H * v_mp
        v_m += H * a_mp
        T += H

        # Passo do Corretor (Equivalente ao STEP=2 / Bloco 55)
        vt_c, vm_c, am_c, beta_c, xnc_c, dxf_c = calculate_derivatives_cmlos(r_t, v_t, a_t, r_m, v_m, beta_t, x_filter)
        
        # Média de Heun para o filtro
        x_filter = 0.5 * (x_filter_old + x_filter + H * dxf_c)

        beta_t = 0.5 * (beta_told + beta_t + H * beta_c)
        r_t = 0.5 * (r_told + r_t + H * vt_c)
        r_m = 0.5 * (r_mold + r_m + H * vm_c)
        v_m = 0.5 * (v_mold + v_m + H * am_c)
        
        rt_history.append(r_t.copy())
        rm_history.append(r_m.copy())
        xnc_history.append(xnc_c)
        time_history.append(T)

        S += H
        if S >= (TS - 0.0001): # TODO understand
            S = 0.
            # XNC/32.2 converte para G's (Screenshot 2026-04-25 201026.png)
            output = f"{T:.3f} | R_T: {r_t/1000} | R_M: {r_m/1000} | Gs: {xnc_c/32.2:.4f}" 
            print(output)
            f.write(output + "\n")

        # Atualização de controle
        r_tm = r_t - r_m
        mod_rtm = np.linalg.norm(r_tm)
        v_tm = vt_c - v_m # Velocidade relativa corrigida
        v_c = -(np.dot(r_tm, v_tm)) / mod_rtm   # TODO entender

print(f"Final: T={T:.3f}, RTM={mod_rtm:.3f}")

rt_history = np.array(rt_history)
rm_history = np.array(rm_history)
xnc_history = np.array(xnc_history)

plt.subplot(1,2,1)
plt.plot(rt_history[:,0]/1000, rt_history[:,1]/1000)
plt.plot(rm_history[:,0]/1000, rm_history[:,1]/1000)
plt.xlabel("X (1000 ft)")
plt.ylabel("Y (1000 ft)")
plt.title("BRM guidance 2D")

plt.subplot(1,2,2)
plt.plot(time_history, xnc_history/32.2)
plt.xlabel("time (s)")
plt.ylabel("XNC (G)")
plt.title("BRM guidance 2D")

plt.show()