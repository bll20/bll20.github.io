'''
Class for creating of instances of a double pendulum.

works in angular coordinates, which are converted to cartesian for output for animating
has methods for iterating the double pendulum in for each induvidual frame, or for returning an array of positions for many frames ahead of time
'''


import numpy as np



g = 0.01

class double_pendulum:
    def __init__(self, params):
        ''' needs to define 2 masses, 2 lengths, and 2 starting positions
        also needs to create arrays for the positions'''

        self.frame = 0

        self.m1 = params[0][0]
        self.m2 = params[1][0]
        self.r1 = params[0][1]
        self.r2 = params[1][1]

        self.a1 = np.array([])
        self.a2 = np.array([])
        self.a1_v = np.array([])
        self.a2_v = np.array([])
        self.a1_a = np.array([])
        self.a2_a = np.array([])

        self.params = params


    def get_masses(self):
        return self.m1, self.m2
    def get_lengths(self):
        return self.r1, self.r2


    def calculate_angular_acceleration(self):
        num1a = -g * (2*self.m1 + self.m2) * np.sin(self.a1[self.frame-1])
        num1b = -self.m2 * g * np.sin(self.a1[self.frame-1] - 2*self.a2[self.frame-1])
        num1c = -2 * np.sin(self.a1[self.frame-1] - self.a2[self.frame-1]) * self.m2 * (self.a2_v[self.frame-1]**2 * self.r2 + self.a1_v[self.frame-1]**2 * self.r1 * np.cos(self.a1[self.frame-1] - self.a2[self.frame-1]))
        den1 = self.r1 * (2*self.m1 + self.m2 - self.m2*np.cos(2*self.a1[self.frame-1] - 2*self.a1[self.frame-1]))
        a1_new_accn = (num1a + num1b + num1c)/(den1)


        num2a = 2 * np.sin(self.a1[self.frame-1] - self.a2[self.frame-1])
        num2ba = self.a1_v[self.frame-1]**2 * self.r1 * (self.m1 + self.m2)
        num2bb = g * (self.m1 + self.m2) * np.cos(self.a1[self.frame-1])
        num2bc = self.a2_v[self.frame-1]**2 * self.r2 * self.m2 * np.cos(self.a1[self.frame-1] - self.a2[self.frame-1])
        den2 = self.r2 * (2*self.m1 + self.m2 - self.m2 * np.cos(2*self.a1[self.frame-1] - 2*self.a2[self.frame-1]))
        a2_new_accn =  num2a*(num2ba + num2bb + num2bc)/(den2)


        return a1_new_accn, a2_new_accn


    def update_pendulum(self):
        ''' frame 0 is the special case where the pendulum is initiated with its inital conditions '''
        if self.frame == 0:
            self.a1 = np.append(self.a1, self.params[0][2])
            self.a2 = np.append(self.a2, self.params[1][2])
            self.a1_v = np.append(self.a1_v, 0)
            self.a2_v = np.append(self.a2_v, 0)
            self.a1_a = np.append(self.a1_a, 0)
            self.a2_a = np.append(self.a2_a, 0)
        else:
            self.a1 = np.append(self.a1, self.a1[self.frame-1] + self.a1_v[self.frame-1])
            self.a2 = np.append(self.a2, self.a2[self.frame-1] + self.a2_v[self.frame-1])
            self.a1_v = np.append(self.a1_v, self.a1_v[self.frame-1] + self.a1_a[self.frame-1])
            self.a2_v = np.append(self.a2_v, self.a2_v[self.frame-1] + self.a2_a[self.frame-1])

            a1_a_new, a2_a_new = self.calculate_angular_acceleration()
            self.a1_a = np.append(self.a1_a, a1_a_new)
            self.a2_a = np.append(self.a2_a, a2_a_new)

        self.frame += 1



    def run(self):
        total_frames = 10000
        for _ in range(1,total_frames):
            self.update_pendulum()
        return self.a1, self.a2


    def iterate(self):   # performs the same functionality as update_pendulum, but returns only the new position
        if self.frame == 0:
            new_pos_1 = self.params[0][2]
            new_pos_2 = self.params[1][2]
            self.a1 = np.append(self.a1, new_pos_1)
            self.a2 = np.append(self.a2, new_pos_2)
            self.a1_v = np.append(self.a1_v, 0)
            self.a2_v = np.append(self.a2_v, 0)
            self.a1_a = np.append(self.a1_a, 0)
            self.a2_a = np.append(self.a2_a, 0)
        else:
            new_pos_1 = self.a1[self.frame-1] + self.a1_v[self.frame-1]
            new_pos_2 = self.a2[self.frame-1] + self.a2_v[self.frame-1]
            self.a1 = np.append(self.a1, new_pos_1)
            self.a2 = np.append(self.a2, new_pos_2)
            self.a1_v = np.append(self.a1_v, self.a1_v[self.frame-1] + self.a1_a[self.frame-1])
            self.a2_v = np.append(self.a2_v, self.a2_v[self.frame-1] + self.a2_a[self.frame-1])

            a1_a_new, a2_a_new = self.calculate_angular_acceleration()
            self.a1_a = np.append(self.a1_a, a1_a_new)
            self.a2_a = np.append(self.a2_a, a2_a_new)

        self.frame += 1

        return new_pos_1, new_pos_2