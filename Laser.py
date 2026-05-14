import numpy as np

class LaserPointer:
    def __init__(self, theta=0,max_dtheta=20.0):
        # Orientação do feixe (em graus)
        self.theta = theta  # Azimute (plano horizontal)
        self.max_dtheta = max_dtheta # # Velocidades angulares máximas (graus por segundo)
        self.beam_radius = 1.0  # Raio máximo do feixe antes de perder o míssil (metros)

    def track_target(self, target_pos, dt):
        """
        Ajusta a orientação do laser para seguir o drone.
        Simplificado: o laser aponta instantaneamente ou com lag proporcional.
        """
        x, y = target_pos
        r = np.sqrt(x**2 + y**2)
        
        # Cálculo dos ângulos alvo para interceptar o drone
        target_theta = np.degrees(np.arctan2(y, x)) # np.arctan() does not distinguishes quartants

        # Cálculo da rotação necessária
        d_theta = target_theta - self.theta

        if d_theta <= self.max_dtheta*dt:
            self.thets = d_theta
        else:
            self.theta = self.max_dtheta*dt

    def get_beam_unit_vector(self):
        """Retorna o vetor unitário na direção do laser."""
        t = np.radians(self.theta)
        return np.array([
            np.cos(t),
            np.sin(t)
        ])

