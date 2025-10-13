
    def chord(self, y):
        self.x_LE = math.tan(math.radians(self.leading_sweep)) * y
        self.x_TE = math.tan(math.radians(self.trailing_sweep)) * y + self.c_root
        return self.x_TE - self.x_LE