class Projectile():
    def __init__(self):
        self.mission_complete = False
        self.alive = True
        self.expl_rad = 2.0
    
    def destroy_target(self):
        self.alive = False
        self.mission_complete = True