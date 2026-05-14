import math
import numpy as np

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

RT1, RT2 = RT1IC, RT2IC
RM1, RM2 = RM1IC, RM2IC
BETAT = 0.
VT1 = -VT * math.cos(BETAT)
VT2 = VT * math.sin(BETAT)
XNC = 0.
T = 0.
S = 0.


# Cálculos iniciais
RTM1 = RT1 - RM1
RTM2 = RT2 - RM2
RTM = math.sqrt(RTM1**2 + RTM2**2)
THETT = math.atan2(RT2, RT1)
VM1 = VM * math.cos(THETT + HE)
VM2 = VM * math.sin(THETT + HE)
VTM1 = VT1 - VM1
VTM2 = VT2 - VM2
VC = -(RTM1 * VTM1 + RTM2 * VTM2) / RTM


# Abre arquivo para escrita
with open('DATFIL.txt', 'w') as f: # TODO verificar se é necessario
    while VC >= 0: # Substitui o IF(VC<0.) GOTO 999
        # Determinação do passo H
        if RTM < 1000:
            H = 0.0002 
        else: 
            H = 0.01
        
        # Salvando estados antigos para integração numérica
        BETATOLD, RT1OLD, RT2OLD = BETAT, RT1, RT2
        RM1OLD, RM2OLD = RM1, RM2
        VM1OLD, VM2OLD = VM1, VM2

        # --- Loop de Integração (Heun/Predictor-Corrector) ---
        # No original usa-se STEP para controlar GOTO 200
        
        def calculate_derivatives(rt1, rt2, rm1, rm2, betat_val):
            # Lógica extraída do bloco 200
            th_t = math.atan2(rt2, rt1)
            th_m = math.atan2(rm2, rm1)
            rtm1_loc = rt1 - rm1
            rtm2_loc = rt2 - rm2
            rtm_loc = math.sqrt(rtm1_loc**2 + rtm2_loc**2)
            rm_loc = math.sqrt(rm1**2 + rm2**2)

            # Velocidades e Acelerações
            vt1_loc = -VT * math.cos(betat_val)
            vt2_loc = VT * math.sin(betat_val)
            vtm1_loc = vt1_loc - VM1 # Nota: VM1/VM2 são atualizados no passo anterior # TODO entender
            vtm2_loc = vt2_loc - VM2
            
            x_lam = math.atan2(rtm2_loc, rtm1_loc)
            x_nc = XNP * rm_loc * (th_t - th_m) # TODO the problem could be the use pf rtm_loc instead of rm_loc --> rtm >> rm --> thus huge accel in the beggining
            am1_loc = -x_nc * math.sin(x_lam)
            am2_loc = x_nc * math.cos(x_lam)
            betatd_loc = XNT / VT
            # x_lam = np.arctan2(r_tm[1], r_tm[0])

            return vt1_loc, vt2_loc, VM1, VM2, am1_loc, am2_loc, betatd_loc, rtm_loc, x_nc

        # Passo do Preditor (Equivalente ao STEP=1 e GOTO 200)
        vt1, vt2, vm1, vm2, am1, am2, betatd, rtm_curr, xnc_curr = calculate_derivatives(RT1, RT2, RM1, RM2, BETAT)
        
        BETAT += H * betatd
        RT1 += H * vt1
        RT2 += H * vt2
        RM1 += H * vm1
        RM2 += H * vm2
        VM1 += H * am1
        VM2 += H * am2
        T += H

        # Passo do Corretor (Equivalente ao STEP=2 / Bloco 55)
        vt1_c, vt2_c, vm1_c, vm2_c, am1_c, am2_c, betatd_c, _, _ = calculate_derivatives(RT1, RT2, RM1, RM2, BETAT)
        
        BETAT = 0.5 * (BETATOLD + BETAT + H * betatd_c)
        RT1 = 0.5 * (RT1OLD + RT1 + H * vt1_c)
        RT2 = 0.5 * (RT2OLD + RT2 + H * vt2_c)
        RM1 = 0.5 * (RM1OLD + RM1 + H * vm1_c)
        RM2 = 0.5 * (RM2OLD + RM2 + H * vm2_c)
        VM1 = 0.5 * (VM1OLD + VM1 + H * am1_c)
        VM2 = 0.5 * (VM2OLD + VM2 + H * am2_c)
        
        S += H
        if S >= (TS - 0.0001):
            S = 0.
            # Output formatado
            rt1k, rt2k = RT1/1000., RT2/1000.
            rm1k, rm2k = RM1/1000., RM2/1000.
            output = f"{T:10.3f}{rt1k:10.3f}{rt2k:10.3f}{rm1k:10.3f}{rm2k:10.3f}{xnc_curr/32.2:10.3f}\n"
            print(output.strip())
            f.write(output)

        # Atualiza VC e RTM para a próxima iteração
        RTM1, RTM2 = RT1 - RM1, RT2 - RM2
        RTM = math.sqrt(RTM1**2 + RTM2**2)
        VTM1, VTM2 = vt1_c - VM1, vt2_c - VM2
        VC = -(RTM1 * VTM1 + RTM2 * VTM2) / RTM

    # Finalização (Bloco 999)
    print(f"Final: T={T}, RTM={RTM}")