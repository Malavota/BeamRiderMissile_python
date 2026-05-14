import numpy as np
import matplotlib.pyplot as plt

def calculate_derivatives_with_leadlag_3d(r_t, r_m, v_m, beta_t, x_filter_theta, x_filter_phi):
    # 1. Geometria do Alvo (Entrada do sistema)
    theta_t = np.arctan2(r_t[1], r_t[0])
    phi_t = np.arctan2(r_t[2], np.sqrt(r_t[0]**2 + r_t[1]**2))
    
    # 2. Geometria do Míssil (Feedback)
    theta_m = np.arctan2(r_m[1], r_m[0])
    phi_m = np.arctan2(r_m[2], np.sqrt(r_m[0]**2 + r_m[1]**2))
    
    mod_rm = np.linalg.norm(r_m)
    r_tm = r_t - r_m
    mod_rtm = np.linalg.norm(r_tm)
    
    # 3. Lei de Guia (Dois canais: Azimute e Elevação)
    # n_c_theta atua no plano XY, n_c_phi atua no plano vertical
    error_theta = mod_rm * (theta_t - theta_m)
    error_phi = mod_rm * (phi_t - phi_m)
    
    # 4. Vetor de Aceleração 3D
    # Direção da Linha de Visada (LOS)
    los_vector = r_tm / mod_rtm # nao usado

    tau_lead = 0.5
    tau_lag = 0.05

    dx_filter_theta = (error_theta - x_filter_theta)/tau_lag
    dx_filter_phi = (error_phi - x_filter_phi)/tau_lag

    x_nc_theta = XNP * (x_filter_theta + tau_lead*dx_filter_theta)
    x_nc_phi = XNP * (x_filter_phi + tau_lead*dx_filter_phi)
    
    # Construção de um sistema de coordenadas local perpendicular à LOS
    # Simplificação: Aceleração em eixos transversais puros
    a_m = np.array([
        -x_nc_theta * np.sin(theta_m),
         x_nc_theta * np.cos(theta_m),
         x_nc_phi # Simplificado para o eixo Z
    ])

    # 5. Movimento do Alvo (Mantendo 2D no plano para simplicidade ou expandindo)
    v_t = VT * np.array([-np.cos(beta_t), np.sin(beta_t), 0.0])
    betatd = XNT / VT

    return v_t, v_m, a_m, betatd, np.sqrt(x_nc_theta**2 + x_nc_phi**2), dx_filter_theta, dx_filter_phi

# --- Configuração Inicial 3D ---
XNP, VM, VT, TS = 10.0, 3000.0, 1000.0, 0.1
XNT = 0.0
r_m = np.array([0.0, 1.0, 0.0])
r_t = np.array([40000.0, 10000.0, 5000.0]) # Adicionado Z = 5000ft
beta_t, T, S = 0.0, 0.0, 0.0

# Velocidade inicial apontando para o alvo
thett_init = np.arctan2(r_t[1], r_t[0])
phit_init = np.arctan2(r_t[2], np.sqrt(r_t[0]**2 + r_t[1]**2))
v_m = VM * np.array([
    np.cos(phit_init) * np.cos(thett_init),
    np.cos(phit_init) * np.sin(thett_init),
    np.sin(phit_init)
])
# Inicializar filtro
x_filter_theta = 0.0
x_filter_phi = 0.0

# Listas de histórico
rm_hist, rt_hist, time_hist, g_hist = [], [], [], []

# Loop de Integração
r_tm = r_t - r_m
v_c = 1.0 # Valor inicial positivo

while v_c >= 0:
    mod_rtm = np.linalg.norm(r_tm)
    H = 0.0002 if mod_rtm < 1000 else 0.01
    
    # Backups para Heun
    r_told, r_mold, v_mold, beta_told = r_t.copy(), r_m.copy(), v_m.copy(), beta_t

    # Preditor
    vt_p, vm_p, am_p, btd_p, g_p, dxf_theta_p, dxf_phi_p = calculate_derivatives_with_leadlag_3d(r_t, r_m, v_m, beta_t, x_filter_theta, x_filter_phi)
    r_t += H * vt_p
    r_m += H * vm_p
    v_m += H * am_p
    beta_t += H * btd_p
    
    # x_filter_old = x_filter
    x_filter_theta_old = x_filter_theta
    x_filter_phi_old = x_filter_phi

    # x_filter += H * dxf_p
    x_filter_theta += H * dxf_theta_p
    x_filter_phi += H * dxf_phi_p

    # Corretor
    vt_c, vm_c, am_c, btd_c, g_c, dxf_theta_c, dxf_phi_c = calculate_derivatives_with_leadlag_3d(r_t, r_m, v_m, beta_t, x_filter_theta, x_filter_phi)
 
    # x_filter = 0.5 * (x_filter_old + H * dxf_c + x_filter)
    x_filter_theta = 0.5 * (x_filter_theta_old + H * dxf_theta_c + x_filter_theta)
    x_filter_theta = 0.5 * (x_filter_phi_old + H * dxf_phi_c + x_filter_phi)

    r_t = 0.5 * (r_told + r_t + H * vt_c)
    r_m = 0.5 * (r_mold + r_m + H * vm_c)
    v_m = 0.5 * (v_mold + v_m + H * am_c)
    beta_t = 0.5 * (beta_told + beta_t + H * btd_c)
    
    T += H
    S += H
    
    # Atualização de controle e v_c
    r_tm = r_t - r_m
    v_tm = vt_c - v_m
    v_c = -(np.dot(r_tm, v_tm)) / np.linalg.norm(r_tm)

    if S >= (TS - 0.0001):
        S = 0.
        rm_hist.append(r_m.copy())
        rt_hist.append(r_t.copy())
        time_hist.append(T)
        g_hist.append(g_c / 32.2)

# --- Plotagem 3D ---
rm_hist, rt_hist = np.array(rm_hist), np.array(rt_hist)
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')
ax.plot(rt_hist[:,0], rt_hist[:,1], rt_hist[:,2], label='Alvo', color='red')
ax.plot(rm_hist[:,0], rm_hist[:,1], rm_hist[:,2], label='Míssil', color='blue', linestyle='--')
ax.set_xlabel('X (ft)')
ax.set_ylabel('Y (ft)')
ax.set_zlabel('Z (ft)')
plt.legend()
plt.show()