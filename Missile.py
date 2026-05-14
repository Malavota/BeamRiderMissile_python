from Projectile import Projectile
from Target import Target
import numpy as np

class Missile():
    def __init__(self, x, y, z, velocidade, theta):
        self.pos = np.array([float(x), float(y)])
        self.v_mag = velocidade
        self.theta = theta
        self.lost = False
    
    def beam_rider_guidance(self, laser_pointer, dt, K=10):
        if self.lost: return

        # 1. Vetor diretor do laser
        theta_laser = laser_pointer.theta
        u_laser = laser_pointer.get_beam_unit_vector()
        
        # 2. Projeção da posição do míssil no eixo do laser
        # Distância percorrida ao longo do feixe
        d_proj = np.dot(self.pos, u_laser)
        ponto_centro_feixe = d_proj * u_laser
        
        # 3. Cálculo do erro (distância perpendicular ao feixe)
        erro_vetor = self.pos - ponto_centro_feixe
        distancia_erro = np.linalg.norm(erro_vetor)

        # 4. Checa se o míssil saiu do feixe (perda de sinal)
        if distancia_erro > laser_pointer.max_beam_radius:
            self.lost = True
            print("Míssil perdido: fora do feixe de laser!")
            return

        # 5. Lógica de correção (proporcional ao quadrado do erro relativo)
        R_m = np.sqrt(self.pos[0]**2 + self.pos[1]**2)# missile radius
        y = R_m*np.radians(self.theta-theta_laser)# missile-beam distance
        n_c = K*y# missile acceleration command

        # 6. Atualização da posição cinemática
        self.theta += n_c*(dt**2)/2 # ?


        
        pass
        