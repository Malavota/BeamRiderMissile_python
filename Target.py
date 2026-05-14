import numpy as np
from Projectile import Projectile

class Target(Projectile):
    def __init__(self, start_pos, target_pos, velocidade):
        self.pos = np.array(start_pos, dtype=float)
        self.goal_pos  = np.array(target_pos, dtype=float)
        self.v_mag = velocidade
        # Calculates unitary direction towards goal
        direcao = self.goal_pos - self.pos
        self.u_dir = direcao / np.linalg.norm(direcao)
        # Drone and mission status

    def update(self, dt):
        # check if mission is complete or drone died already
        if not self.alive or self.mission_complete:
            return
        
        # check if goal was reached, otherwise update position
        dist_ao_alvo = np.linalg.norm(self.target - self.pos)
        if dist_ao_alvo < self.exp_rad:
            self.destroy_target()
        else:
            self.pos += self.u_dir * self.v_mag * dt

