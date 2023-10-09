'''
Simple animation of a double pendulum using python multiprocessing
import double pendulum class from double_pendulum.py
creates an arbitrary number of double pendulum animations with small variations in starting parameters
performance is improved through use of multiprocessing
'''



import matplotlib.pyplot as plt
import numpy as np
import multiprocessing as mp
from double_pendulum_class import double_pendulum
import random

def polar_to_cartesian(a1,a2, r1,r2):
    x1 = r1 * np.sin(a1)
    y1 = r2 * np.cos(a1)
    x2 = x1 + r1 * np.sin(a2)
    y2 = y1 + r2 * np.cos(a2)

    return x1,y1, x2,y2




fig,ax = plt.subplots()



def create_animation():

    perturbation = random.randint(0,10)*np.pi/100
    params = [[10,10,np.pi/2],[15,10,np.pi/2+perturbation]]
    pendulum = double_pendulum(params)

    ax.set_box_aspect(1)
    size = 1.1 * (params[0][1] + params[1][1])


    while True:

        a1_new, a2_new = pendulum.iterate()
        x1_new,y1_new, x2_new,y2_new = polar_to_cartesian(a1_new,a2_new, params[0][1],params[1][1])


        ax.clear()
        ax.plot([0,x1_new],[0,y1_new],'-',color='k')   # line 1
        ax.plot([x1_new, x2_new], [y1_new, y2_new],'-',color='k')  # line 2
        ax.plot(x1_new, y1_new, '.', markersize = 10) # mass 1
        ax.plot(x2_new, y2_new, '.', markersize = 15) # mass 2
        # ax.plot(x2[:f], y2[:f],color='red',alpha=0.5)


        ax.set_xlim([-size, size])
        ax.set_ylim([size, -size/2])



        fig.canvas.blit(ax.bbox)
        plt.pause(0.01)


if __name__ == '__main__':

    processes = []
    for _ in range(15):
        processes.append(mp.Process(target=create_animation, daemon=True))

    print(processes)

    for process in processes:
        process.start()

    try:
        while True:
            _=''
    except KeyboardInterrupt:      # main process does not end until keyboard interrupt (control C)
        pass
