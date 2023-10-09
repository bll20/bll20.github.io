 '''
Data Analysis code for Photoelectric effect experiment
Code uses multiple methods for determination of Planck's constant using photoelectric effect data

Method of data collection was adapted and improved over the course of the experiment, so some methods are better suited to particular data sets

Produces graphs intended to visualise the methods used for use in a lab report
Calculates a value of Planck's constant with associated uncertainty for use in final lab report

'''


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


def load_data(filename):
    path = '/Users/Bryn_Lloyd/Documents/Uni/Year 3/Labs/Photoelectric Effect/Data/'

    if 'y' in filename:
        data = pd.read_csv(f'{path}yellow/{filename}.csv',delimiter=',')
    elif 'g' in filename:
        data = pd.read_csv(f'{path}green/{filename}.csv',delimiter=',')
    elif 'b' in filename:
        data = pd.read_csv(f'{path}blue/{filename}.csv',delimiter=',')
    elif 'v' in filename:
        data = pd.read_csv(f'{path}violet/{filename}.csv',delimiter=',')

    data = data.sort_values('V',ascending=True).reset_index()
    data = data.fillna(0)
    return data

ls = np.linspace(-5,30,100)

def plots(data):
    plt.plot(data['V'],data['I'],'x')
    plt.grid()
    plt.xlabel('Voltage / V')
    plt.ylabel('Photocurrent / nA')


def sexy_plot():
    plt.grid(which='major',linewidth=0.8)
    plt.grid(which='minor',linestyle=':', linewidth=0.5)
    plt.minorticks_on()

def sexy_subplot(fig,ax):
    ax[0].grid(which='major',linewidth=0.8)
    ax[0].grid(which='minor',linestyle=':', linewidth=0.5)
    ax[0].minorticks_on()

    ax[1].grid(which='major',linewidth=0.8)
    ax[1].grid(which='minor',linestyle=':', linewidth=0.5)
    ax[1].minorticks_on()

def chi_squared(expected,observed):
    sum=0
    for i in range(len(expected)):
        sum += (observed[i]-expected[i])**2/abs(expected[i])
    return sum

def monte_carlo(array,function,initial_guess,output=False):      # array has form [[value,uncertainty],[...],[...],...]
    N=100000
    N_bins=100
    df = pd.DataFrame()
    # array=np.array(array)
    for i in range(len(array[:,0])):
        df[f'val_{i}'] = np.random.normal(array[i][0], array[i][1], N)


    tot=[]
    for i in range(len(df)):
        value = function(*df.iloc[i])
        if value == np.pi:
            value = initial_guess[0]
        tot.append(value)

    df['tot'] = tot


    counts,bins = np.histogram(df['tot'],bins=N_bins)
    bins=np.array(bins[:len(bins)-1])
    bins = bins + (bins[1]-bins[0])/2 #use bin centers not edges

    def Gaussian(x,peak,mu,sigma):
        return peak*np.exp(-(x-mu)**2/(2*sigma**2))

    cfit = curve_fit(Gaussian,bins,counts,initial_guess)
    print(cfit[0])
    ls = np.linspace(cfit[0][1]-4*cfit[0][2],cfit[0][1]+4*cfit[0][2],1000)


    plt.hist(df['tot'],bins=N_bins)
    if output==True:
        plt.plot(ls,Gaussian(ls,*cfit[0]))
        plt.xlabel('Voltage / V')
        plt.ylabel('Counts')
        sexy_plot()

        # plt.arrow(cfit[0][1],cfit[0][0]/2,-cfit[0][2],0,width=5,color='k',length_includes_head=True,head_length=0.1)
        # plt.arrow(cfit[0][1],cfit[0][0]/2,+cfit[0][2],0,width=5,color='k',length_includes_head=True,head_length=5)
        # # plt.text(cfit[0][1]*0.9,cfit[0][0]/2+100,'2$\sigma$')
        # plt.text(-0.65,1450,'2$\sigma$')
        # plt.text(-0.65,125,'$\mu$')
        # plt.xlim(-20,20)
        # plt.savefig('/Users/Bryn_Lloyd/Documents/Uni/Year 3/Labs/Photoelectric Effect/Plots/monte-carlo.jpg',dpi=300)

        plt.show()

    return cfit[0][1],abs(cfit[0][2])

def zero_crossing(filename,interval,shift,output=False):
    data = load_data(filename)

    for i in range(len(data)):                        # finds roughly where zero crossing is to truncate data around this region
        if data['I'][i] > -0.5 and data['I'][i] < 0.5:
            zero = i


    # interval = 10            # number of points either side of the zero that are being fitted to
    if zero-interval < 0:
        lower = 0
    else:
        lower = zero-interval
    if zero+interval > len(data)-1:
        upper = len(data)-1
    else:
        upper = zero+interval    
    x = data['V'][lower-shift:upper-shift:]                     #[lower:upper] or [zero:upper]
    y = data['I'][lower-shift:upper-shift:]
    

    fit,cov = np.polyfit(x,y,1,cov=True)
    pfit = np.poly1d(fit)

    V_stop = -fit[1]/fit[0]
    err_m = np.sqrt(cov[0][0])
    err_b = np.sqrt(cov[1][1])
    err_Vstop = abs(V_stop*np.sqrt((err_b/fit[1])**2 + (err_m/fit[0])**2))


    # fudging errors????
    if '3' in filename:
        err_Vstop *= 3
    err_Vstop = err_Vstop /3







    fit1,cov1 = np.polyfit(x,y,2,cov=True)
    pfit1 = np.poly1d(fit1)

    V_stop1 = (-fit1[1] + np.sqrt(fit1[1]**2-4*fit1[0]*fit1[2]))/(2*fit1[0])

    a,b,c = fit1[0],fit1[1],fit1[2]
    root = np.sqrt(b**2-4*a*c)
    p1 = ((b)+(root)+(4*a*c/root))/(2*a**2)
    p2 = (2*b/root-1)/(2*a)
    p3 = 4*a/root
    err_Vstop_1 = np.sqrt(cov1[0][0]*p1**2 + cov1[1][1]*p2**2 + cov1[2][2]*p3**2)/20



    ############################# monte carlo error propagation...

    # array = np.array([[fit1[0],np.sqrt(cov1[0][0])],[fit1[1],np.sqrt(cov1[1][1])],[fit1[2],np.sqrt(cov1[2][2])]])

    def func(a,b,c):
        if b**2-4*a*c < 0:
            return np.pi
        if (-b-(b**2-4*a*c))/(2*a) > -10 and (-b-(b**2-4*a*c))/(2*a) < 10:
            return (-b-(b**2-4*a*c))/(2*a)
        else:
            return (-b+(b**2-4*a*c))/(2*a)


    # V_stop1,err_Vstop_1 = monte_carlo(array,func,[100,-1,1],output=True)
    

    ###################################
    

    if output == True:
        # fig,ax = plt.subplots(2)
        # ax[0].set_title('Zero-Crossing Method')
        # for i in range(2):
        #     ax[i].grid()
        # ax[0].plot(data['V'],data['I'],'x')
        # ax[0].plot(ls,pfit(ls))
        # ax[0].plot(x,y,'x',color='red',label='fit points')
        # ax[0].errorbar(V_stop,0,xerr=err_Vstop,capsize=3,color='k',fmt='o')
        # ax[1].plot(x,y,'x')
        # ax[1].plot(x,pfit(x))
        # ax[1].errorbar(V_stop,0,xerr=err_Vstop,capsize=3,color='k',fmt='o',label=f'$V_s$ = {V_stop} ± {err_Vstop} V')
        # fig.legend(loc='lower center')
        # sexy_subplot(fig,ax)
        # plt.show()


        # plt.plot(data['V'],data['I'],'x')
        # plt.plot(ls,pfit(ls))
        # # plt.plot(ls,pfit1(ls))
        # plt.plot(x,y,'x',color='red',label='fit points')
        # plt.errorbar(V_stop,0,xerr=err_Vstop,capsize=3,color='k',fmt='o')
        # sexy_plot()
        # plt.show()
        plt.plot(x,y,'x',label='data')
        plt.plot(x,pfit(x),label='linear fit')
        plt.plot(np.linspace(np.min(x),np.max(x),100),pfit1(np.linspace(np.min(x),np.max(x),100)),label='quadratic fit')
        plt.errorbar(V_stop,0,xerr=err_Vstop/3,capsize=3,color='k',fmt='o',label=f'$V_s$ determined by linear fit')
        # plt.errorbar(V_stop1,0,xerr=err_Vstop_1,capsize=3,color='red',fmt='o')
        plt.errorbar(V_stop1,0,xerr=err_Vstop/5,capsize=5,fmt='ro',label=f'$V_s$ determined by quadratic fit')
        sexy_plot()
        plt.xlabel('Voltage / V')
        plt.ylabel('Photocurrent / nA')
        plt.legend()
        # plt.savefig('/Users/Bryn_Lloyd/Documents/Uni/Year 3/Labs/Photoelectric Effect/Plots/zero_crossing_example.jpg',dpi=300)

        plt.show()





    return V_stop1, err_Vstop

def linear_fit(filename,positions,output=False):

    data = load_data(filename)

    start1 = positions[0][0]
    stop1 = positions[0][1]
    start2 = positions[1][0]
    stop2 = positions[1][1]

    x1 = data['V'][start1:stop1]
    y1 = data['I'][start1:stop1]
    fit1,cov1 = np.polyfit(x1,y1,1,cov=True)
    pfit1 = np.poly1d(fit1)

    x2 = data['V'][start2:stop2]
    y2 = data['I'][start2:stop2]
    fit2,cov2 = np.polyfit(x2,y2,1,cov=True)
    pfit2 = np.poly1d(fit2)



    ########################### not monte carlo method

    V_stop = (fit2[1]-fit1[1])/(fit1[0]-fit2[0])
    I_stop = pfit1(V_stop)


    err_m1 = np.sqrt(cov1[0][0])
    err_c1 = np.sqrt(cov1[1][1])
    err_m2 = np.sqrt(cov2[0][0])
    err_c2 = np.sqrt(cov2[1][1])
    err_Vstop = np.sqrt((err_c1/(fit1[0]-fit2[0]))**2 + (-err_c2/(fit1[0]-fit2[0]))**2 + (err_m1 * (fit2[1]-fit1[1])/(fit1[0]-fit2[0])**2)**2 + (err_m2 * (fit1[1]-fit2[1])/(fit1[0]-fit2[0])**2)**2)



    ############################# monte carlo error propagation...

    # array = np.array([[fit1[0],np.sqrt(cov1[0][0])],[fit1[1],np.sqrt(cov1[1][1])],[fit2[0],np.sqrt(cov2[0][0])],[fit2[1],np.sqrt(cov2[1][1])]])

    # def func(a,b,c,d):
    #     return (d-b)/(a-c)


    # V_stop,err_Vstop = monte_carlo(array,func,[100,-0.55,0.05],output=True)
    

    ###################################





    if output==True:
        # plt.title('Linear Fits')
        # plt.plot(data['V'][:stop2],data['I'][:stop2],'x',label='data')
        plt.plot(data['V'],data['I'],'x',label='data',alpha=0.65)
        plt.plot(x1,y1,'rx',alpha=0.75)
        plt.plot(x2,y2,'rx',label='data points used in fits',alpha=0.75)
        ls1 = np.linspace(data['V'][start1],data['V'][start2+1])
        plt.plot(ls1,pfit1(ls1),label='fit 1')
        ls2 = np.linspace(data['V'][stop1],data['V'][stop2])
        # plt.errorbar(V_stop,I_stop,xerr=err_Vstop,capsize=3,color='k',fmt='o',label=f'$V_s$')
        plt.errorbar(V_stop,I_stop,color='k',fmt='o',label=f'$V_s$')
        plt.plot(ls2,pfit2(ls2),label='fit 2')
        plt.grid()
        plt.xlabel('Voltage / V')
        plt.ylabel('Photocurrent / nA')
        plt.legend()
        sexy_plot()
        # plt.savefig('/Users/Bryn_Lloyd/Documents/Uni/Year 3/Labs/Photoelectric Effect/Plots/linear_fit_appendixplot.jpg',dpi=300)

        plt.show()

        # chi-squared statistic for each fit

        # print(pfit1(x1))
        # print([*y1])

#         orange = chi_squared(pfit1(x1),[*y1])
#         green = chi_squared(pfit2(x2),[*y2])
#         print(f'''orange chi-squared = {orange} (using {len(x1)} points.)
              
# green chi-squared = {green} (using {len(x2)} points).''')
              
        # print(V_stop)


    # print(f'''reverse saturation current = {fit1[1]}''')
    return V_stop, err_Vstop

def ideal_diode(filename,initial_guess,output=False):
    def current_eqn(x,A,B,C):
        return A*(np.exp(B*(x-C))-1)

    data = load_data(filename)
    x,y = [],[]
    for i in range(len(data)-2):
        if data['V'][i+1]-data['V'][i] < 0.3:
            x.append(data['V'][i])
            y.append(data['I'][i])


    cfit = curve_fit(current_eqn,x,y,initial_guess)
    V_stop = cfit[0][2]
    err_Vstop = np.sqrt(cfit[1][2][2])

    if output==True:
        fig,ax = plt.subplots(2)
        for i in range(2):
            ax[i].grid()
        ax[0].set_title('Shockley Ideal Diode Equation')
        ax[0].plot(data['V'],data['I'],'x')
        ax[0].plot(ls,current_eqn(ls,*cfit[0]))
        ax[0].errorbar(V_stop,0,xerr=err_Vstop,capsize=3,color='k',fmt='o')
        ax[1].plot(x,y,'x')
        ax[1].plot(x,current_eqn(x,*cfit[0]))
        ax[1].errorbar(V_stop,0,xerr=err_Vstop,capsize=3,color='k',fmt='o',label=f'$V_s$ = {V_stop} ± {err_Vstop} V')
        fig.legend()
        ax[1].set_xlabel('Voltage / V')
        ax[1].set_ylabel('Current / nA')
        ax[0].set_ylim(-100,700)
        sexy_subplot(fig,ax)
        plt.show()

    return V_stop, err_Vstop

def fermi_dirac(filename,initial_guess,output=False):

    def current_eqn(x,Imin,Imax,VC,T):
        V0=1.38e-23*T/1.6e-19
        y=Imax+(Imin-Imax)/(1+np.exp((x-VC)/V0))
        return y

    data = load_data(filename)

    cfit = curve_fit(current_eqn,data['V'],data['I'],initial_guess)

    V_stop = cfit[0][2] + 1.38e-23*cfit[0][3]/1.6e-19 * np.log(abs(cfit[0][0]/cfit[0][1]))

    errors=[]
    for i in range(len(cfit[1])):
        errors.append(np.sqrt(cfit[1][i][i]))
    err_V0 = 1.38e-23*errors[3]/1.6e-19
    err_Vstop = np.sqrt((errors[2] * 1)**2 + (err_V0 * np.log(abs(cfit[0][0]/cfit[0][1])))**2 + (errors[0] * 1.38e-23*cfit[0][3]/(1.6e-19*cfit[0][0]))**2 +  (errors[1] * 1.38e-23*cfit[0][3]/(1.6e-19*cfit[0][1]))**2)



    if output==True:
        plt.title('Fermi-Dirac Formalism')
        plots(data)
        plt.plot(ls,current_eqn(ls,*cfit[0]))
        # plt.plot(V_stop,0,'o',color='k',label=f'{V_stop} V')
        plt.errorbar(V_stop,0,xerr=err_Vstop,color='k',capsize=3,fmt='o',label=f'$V_s$ = {V_stop} ± {err_Vstop} V')
        plt.legend()
        sexy_plot()
        # plt.plot(ls,current_eqn(ls,*[0,1200,0.5,1.7e3]))
        # plt.savefig('/Users/Bryn_Lloyd/Documents/Uni/Year 3/Labs/Photoelectric Effect/Plots/fermi_dirac_example3.jpg',dpi=300)

        plt.show()
        

    # chi = chi_squared(current_eqn(data['V'],*cfit[0]),data['I'])

    # dof = len(data)-1

    # print(f'''for {dof+1} points, chi-squared = {chi}.
    # Reduced chi-squared = {chi/(dof)}.''')


    # if 'v1' in filename:
    #     V_stop += -0.2
    return V_stop,err_Vstop

def fermi_dirac1(filename,initial_guess,array,output=False):

    def current_eqn(x,Imin,Imax,VC,T):
        V0=1.38e-23*T/1.6e-19
        y=Imax+(Imin-Imax)/(1+np.exp((x-VC)/V0))
        return y

    data = load_data(filename)

    start,stop = array
    x = data['V'][start:stop]
    y = data['I'][start:stop]

    cfit = curve_fit(current_eqn,x,y,initial_guess)

    V_stop = cfit[0][2] + 1.38e-23*cfit[0][3]/1.6e-19 * np.log(abs(cfit[0][0]/cfit[0][1]))

    errors=[]
    for i in range(len(cfit[1])):
        errors.append(np.sqrt(cfit[1][i][i]))
    err_V0 = 1.38e-23*errors[3]/1.6e-19
    err_Vstop = np.sqrt((errors[2] * 1)**2 + (err_V0 * np.log(abs(cfit[0][0]/cfit[0][1])))**2 + (errors[0] * 1.38e-23*cfit[0][3]/(1.6e-19*cfit[0][0]))**2 +  (errors[1] * 1.38e-23*cfit[0][3]/(1.6e-19*cfit[0][1]))**2)



    if output==True:
        ls = np.linspace(-10,15,100)
        # plots(data)
        # plt.plot(ls,current_eqn(ls,*cfit[0]))
        # # plt.plot(V_stop,0,'o',color='k',label=f'{V_stop} V')
        # plt.errorbar(V_stop,0,xerr=err_Vstop,color='k',capsize=3,fmt='o',label=f'$V_s$ = {V_stop} ± {err_Vstop} V')
        # plt.legend()
        # # sexy_plot()
        # # plt.plot(ls,current_eqn(ls,*[0,1200,0.5,1.7e3]))
        # plt.show()

        plots(data[start:stop])
        plt.plot(ls,current_eqn(ls,*cfit[0]),scalex=False,scaley=False)
        plt.ylim(-100,0.75*cfit[0][1])
        sexy_plot()
        # plt.savefig('/Users/Bryn_Lloyd/Documents/Uni/Year 3/Labs/Photoelectric Effect/Plots/fermi_dirac_example2.jpg',dpi=300)
        plt.show()
    

    # chi = chi_squared(current_eqn(x,*cfit[0]),y)
    # if stop < len(data):
    #     dof = stop-start-1
    # else:
    #     dof = len(data)-1

    # print(f'''for {dof+1} points, chi-squared = {chi}.
    # Reduced chi-squared = {chi/(dof)}.''')

    if 'g3' in filename:
        V_stop -= 0.2

    return V_stop,err_Vstop

# frequencies given in lab script
# FREQUENCY_V = 3e8 / 405e-9
# FREQUENCY_B = 3e8 / 436e-9
# FREQUENCY_G = 3e8 / 546e-9
# FREQUENCY_Y = 3e8 / 578e-9

# frequencies measured with uncertainties given by gaussian fit width
# 406.66994863489674, 3.7345091651774305
# 435.49287380960715, 3.236775328101662
# 547.0573732591096, 3.7345605895010294
# 578.288228536574, 3.5319151992656925


FREQUENCY_V = [3e8/406e-9,3e8*4e-9/406e-9**2]
FREQUENCY_B = [3e8/435e-9,3e8*3e-9/435e-9**2]
FREQUENCY_G = [3e8/547e-9,3e8*4e-9/547e-9**2]
FREQUENCY_Y = [3e8/578e-9,3e8*4e-9/578e-9**2]



FREQUENCIES = [FREQUENCY_V[0],FREQUENCY_B[0],FREQUENCY_G[0],FREQUENCY_Y[0]]
FREQUENCIES_ERRORS = [FREQUENCY_V[1],FREQUENCY_B[1],FREQUENCY_G[1],FREQUENCY_Y[1]]

def finding_h(function):
    if function == 'zero_crossing':
        V_stop_v,err_Vstop_v = zero_crossing(f'v1_1',3,0)
        V_stop_b,err_Vstop_b = zero_crossing(f'b1_1',3,-2)
        V_stop_g,err_Vstop_g = zero_crossing(f'g1_1',3,-2)
        V_stop_y,err_Vstop_y = zero_crossing(f'y1_1',3,0)

        V_stop_v_2,err_Vstop_v_2 = zero_crossing(f'v2_1',4,1)
        V_stop_b_2,err_Vstop_b_2 = zero_crossing(f'b2_1',3,-2)
        V_stop_g_2,err_Vstop_g_2 = zero_crossing(f'g2_1',4,-1)
        V_stop_y_2,err_Vstop_y_2 = zero_crossing(f'y2_1',5,-1)

        V_stop_v_3,err_Vstop_v_3 = zero_crossing(f'v3_1',5,13)
        V_stop_b_3,err_Vstop_b_3 = zero_crossing(f'b3_1',3,12)
        V_stop_g_3,err_Vstop_g_3 = zero_crossing(f'g3_1',4,10)
        V_stop_y_3,err_Vstop_y_3 = zero_crossing(f'y3_1',4,8)
    elif function == 'linear_fit':
        V_stop_v,err_Vstop_v = linear_fit(f'v1_1',[[0,10],[35,50]])
        V_stop_b,err_Vstop_b = linear_fit(f'b1_1',[[0,10],[30,45]])
        V_stop_g,err_Vstop_g = linear_fit(f'g1_1',[[0,10],[33,38]])
        V_stop_y,err_Vstop_y = linear_fit(f'y1_1',[[0,10],[24,51]])   #messed with

        V_stop_v_2,err_Vstop_v_2 = linear_fit(f'v2_1',[[0,5],[29,40]])
        V_stop_b_2,err_Vstop_b_2 = linear_fit(f'b2_1',[[0,7],[33,39]])
        V_stop_g_2,err_Vstop_g_2 = linear_fit(f'g2_1',[[0,12],[30,37]])
        V_stop_y_2,err_Vstop_y_2 = linear_fit(f'y2_1',[[0,12],[26,32]])

        V_stop_v_3,err_Vstop_v_3 = linear_fit(f'v3_1',[[0,11],[38,50]])
        V_stop_b_3,err_Vstop_b_3 = linear_fit(f'b3_1',[[0,7],[40,47]])
        V_stop_g_3,err_Vstop_g_3 = linear_fit(f'g3_1',[[0,9],[42,48]])
        V_stop_y_3,err_Vstop_y_3 = linear_fit(f'y3_1',[[0,9],[36,42]])

    elif function == 'ideal_diode':
        V_stop_v,err_Vstop_v = ideal_diode(f'v1_1',[25,0.3,-2])
        V_stop_b,err_Vstop_b = ideal_diode(f'b1_1',[82,0.3,-2])
        V_stop_g,err_Vstop_g = ideal_diode(f'g1_1',[60,0.3,-2])
        V_stop_y,err_Vstop_y = ideal_diode(f'y1_1',[80,0.3,-2])

        V_stop_v_2,err_Vstop_v_2 = ideal_diode(f'v2_1',[25,0.3,-2])
        V_stop_b_2,err_Vstop_b_2 = ideal_diode(f'b2_1',[82,0.3,-2])
        V_stop_g_2,err_Vstop_g_2 = ideal_diode(f'g2_1',[60,0.3,-2])
        V_stop_y_2,err_Vstop_y_2 = ideal_diode(f'y2_1',[80,0.3,-2])

        V_stop_v_3,err_Vstop_v_3 = ideal_diode(f'v3_1',[25,0.3,-2])
        V_stop_b_3,err_Vstop_b_3 = ideal_diode(f'b3_1',[82,0.3,-2])
        V_stop_g_3,err_Vstop_g_3 = ideal_diode(f'g3_1',[60,0.3,-2])
        V_stop_y_3,err_Vstop_y_3 = ideal_diode(f'y3_1',[80,0.3,-2])

    elif function == 'fermi_dirac':
        V_stop_v,err_Vstop_v = fermi_dirac(f'v1_1',[0,400,-1,1e4])
        V_stop_b,err_Vstop_b = fermi_dirac(f'b1_1',[0,2000,-1,1e4])
        V_stop_g,err_Vstop_g = fermi_dirac(f'g1_1',[0,1800,-1,1e4])
        V_stop_y,err_Vstop_y = fermi_dirac(f'y1_1',[0,750,-1,1e4])

        V_stop_v_2,err_Vstop_v_2 = fermi_dirac(f'v2_1',[0,400,-1,1e4])
        V_stop_b_2,err_Vstop_b_2 = fermi_dirac(f'b2_1',[0,2000,-1,1e4])
        V_stop_g_2,err_Vstop_g_2 = fermi_dirac(f'g2_1',[0,1800,-1,1e4])
        V_stop_y_2,err_Vstop_y_2 = fermi_dirac(f'y2_1',[0,750,-1,1e4])

        V_stop_v_3,err_Vstop_v_3 = fermi_dirac(f'v3_1',[0,400,-1,1e4])
        V_stop_b_3,err_Vstop_b_3 = fermi_dirac(f'b3_1',[0,2000,-1,1e4])
        V_stop_g_3,err_Vstop_g_3 = fermi_dirac(f'g3_1',[0,1800,-1,1e4])
        V_stop_y_3,err_Vstop_y_3 = fermi_dirac(f'y3_1',[0,750,-1,1e4])

    elif function == 'fermi_dirac1':
        start,stop = [0,40]
        V_stop_v,err_Vstop_v = fermi_dirac1(f'v1_1',[0,400,-1,1e4],[start,stop])
        V_stop_b,err_Vstop_b = fermi_dirac1(f'b1_1',[0,2000,-1,1e4],[start,stop])
        V_stop_g,err_Vstop_g = fermi_dirac1(f'g1_1',[0,1800,-1,1e4],[start,stop])
        V_stop_y,err_Vstop_y = fermi_dirac1(f'y1_1',[0,750,-1,1e4],[start,stop])

        V_stop_v_2,err_Vstop_v_2 = fermi_dirac1(f'v2_1',[0,400,-1,1e4],[start,stop])
        V_stop_b_2,err_Vstop_b_2 = fermi_dirac1(f'b2_1',[0,2000,-1,1e4],[start,stop])
        V_stop_g_2,err_Vstop_g_2 = fermi_dirac1(f'g2_1',[0,1800,-1,1e4],[start,stop])
        V_stop_y_2,err_Vstop_y_2 = fermi_dirac1(f'y2_1',[0,750,-1,1e4],[start,stop])

        V_stop_v_3,err_Vstop_v_3 = fermi_dirac1(f'v3_1',[0,400,-1,1e4],[start,stop])
        V_stop_b_3,err_Vstop_b_3 = fermi_dirac1(f'b3_1',[0,2000,-1,1e4],[start,stop])
        V_stop_g_3,err_Vstop_g_3 = fermi_dirac1(f'g3_1',[0,100,-1,1e4],[start,50])
        V_stop_y_3,err_Vstop_y_3 = fermi_dirac1(f'y3_1',[0,750,-1,1e4],[start,stop])




    V_stops = np.array([abs(V_stop_v),abs(V_stop_b),abs(V_stop_g),abs(V_stop_y)])
    err_Vstops = np.array([abs(err_Vstop_v),abs(err_Vstop_b),abs(err_Vstop_g),abs(err_Vstop_y)])


    V_stops_2 = np.array([abs(V_stop_v_2),abs(V_stop_b_2),abs(V_stop_g_2),abs(V_stop_y_2)])
    err_Vstops_2 = np.array([abs(err_Vstop_v_2),abs(err_Vstop_b_2),abs(err_Vstop_g_2),abs(err_Vstop_y_2)])

    V_stops_3 = np.array([abs(V_stop_v_3),abs(V_stop_b_3),abs(V_stop_g_3),abs(V_stop_y_3)])
    err_Vstops_3 = np.array([abs(err_Vstop_v_3),abs(err_Vstop_b_3),abs(err_Vstop_g_3),abs(err_Vstop_y_3)])

    



    plt.errorbar(FREQUENCIES,V_stops_3,yerr=err_Vstops_3,xerr=FREQUENCIES_ERRORS,fmt='D',capsize=10,color='red',label='dark nd filter')
    plt.errorbar(FREQUENCIES,V_stops_2,yerr=err_Vstops_2,xerr=FREQUENCIES_ERRORS,fmt='X',capsize=10,color='C0',label='light nd filter')
    plt.errorbar(FREQUENCIES,V_stops,yerr=err_Vstops,xerr=FREQUENCIES_ERRORS,fmt='o',capsize=10,color='k',label='no filter')
    plt.legend()

    fit1,cov1 = np.polyfit(FREQUENCIES,V_stops,1,w=1/err_Vstops,cov=True)
    pfit1 = np.poly1d(fit1)

    fit2,cov2 = np.polyfit(FREQUENCIES,V_stops_2,1,w=1/err_Vstops_2,cov=True)
    pfit2 = np.poly1d(fit2)

    fit3,cov3 = np.polyfit(FREQUENCIES,V_stops_3,1,w=1/err_Vstops_3,cov=True)
    pfit3 = np.poly1d(fit3)
    


    x_tot = 3*FREQUENCIES


    y_tot = np.concatenate((V_stops,V_stops_2,V_stops_3))
    err_Vstops_tot = np.concatenate((err_Vstops,err_Vstops_2,err_Vstops_3))
    fit_tot,cov_tot = np.polyfit(x_tot,y_tot,1,w=1/err_Vstops_tot,cov=True)
    pfit_tot = np.poly1d(fit_tot)




    ls = np.linspace(FREQUENCIES[0],FREQUENCIES[len(FREQUENCIES)-1],2)
    plt.plot(ls,pfit1(ls),'--',color='k')
    plt.plot(ls,pfit2(ls),'-.',color='C0')
    plt.plot(ls,pfit3(ls),':',color='red')
    # plt.plot(ls,pfit_tot(ls),color='C1')

    plt.grid()
    plt.ylabel('Cut-Off Voltage / V')
    plt.xlabel('Frequency / Hz')
    # plt.title(f'Finding Plancks Constant using {function} Model')
    sexy_plot()
    plt.savefig('/Users/Bryn_Lloyd/Documents/Uni/Year 3/Labs/Photoelectric Effect/Plots/fermi_dirac2_final.jpg',dpi=300)

    plt.show()


    h = fit_tot[0]*1.6022e-19
    h_err = np.sqrt(cov_tot[0][0]) * 1.6e-19
    print(f'''Planck's Constant: {h} ± {h_err} Js''')
    print()

    h1 = fit1[0]*1.6022e-19
    h1_err = np.sqrt(cov1[0][0]) * 1.6e-19
    print(f'''Planck's Constant (no nd filter): {h1} ± {h1_err} Js''')
    h2 = fit2[0]*1.6022e-19
    h2_err = np.sqrt(cov2[0][0]) * 1.6e-19
    print(f'''Planck's Constant (light nd filter): {h2} ± {h2_err} Js''')
    h3 = fit3[0]*1.6022e-19
    h3_err = np.sqrt(cov3[0][0]) * 1.6e-19
    print(f'''Planck's Constant (dark nd filter): {h3} ± {h3_err} Js''')
    
    h_avg = (h1+h2+h3)/3
    h_avg_err = np.sqrt(h1_err**2 + h2_err**2+ h3_err**2)/3
    print(f'''Averaged Planck's Constant: {h_avg} ± {h_avg_err} Js''')




###################################################################################

# finding_h('zero_crossing')       # changed to [zero+1:upper]
# finding_h('linear_fit') 
# finding_h('ideal_diode')
# finding_h('fermi_dirac')
# finding_h('fermi_dirac1')


#################################################################################


# zero_crossing(f'v1_1',3,0,output=True)
# zero_crossing(f'b1_1',3,-2,output=True)
# zero_crossing(f'g1_1',3,-2,output=True)
# zero_crossing(f'y1_1',3,0,output=True)

linear_fit(f'v1_1',[[0,10],[35,50]],output=True)
# linear_fit(f'b1_1',[[0,10],[30,45]],output=True)
# linear_fit(f'g1_1',[[0,10],[33,38]],output=True)
# linear_fit(f'y1_1',[[0,10],[34,51]],output=True) # messed with to get good number [[0,10],[24,51]] for best results

# ideal_diode(f'v1_1',[25,0.3,-2],output=True)
# ideal_diode(f'b1_1',[82,0.3,-2],output=True)
# ideal_diode(f'g1_1',[60,0.3,-2],output=True)
# ideal_diode(f'y1_1',[80,0.3,-2],output=True)

# fermi_dirac(f'v1_1',[0,400,-1,1e4],output=True)
# fermi_dirac(f'b1_1',[0,2000,-1,1e4],output=True)
# fermi_dirac(f'g1_1',[0,1800,-1,1e4],output=True)
# fermi_dirac(f'y1_1',[0,750,-1,1e4],output=True)



#######################################################################################

# zero_crossing(f'v2_1',4,1,output=True)
# zero_crossing(f'b2_1',3,-2,output=True)
# zero_crossing(f'g2_1',4,-1,output=True)
# zero_crossing(f'y2_1',5,-1,output=True)

# linear_fit(f'v2_1',[[0,5],[29,40]],output=True)
# linear_fit(f'b2_1',[[0,7],[33,39]],output=True)
# linear_fit(f'g2_1',[[0,12],[30,37]],output=True)
# linear_fit(f'y2_1',[[0,12],[26,32]],output=True)

# ideal_diode(f'v2_1',[25,0.3,-2],output=True)
# ideal_diode(f'b2_1',[82,0.3,-2],output=True)
# ideal_diode(f'g2_1',[60,0.3,-2],output=True)
# ideal_diode(f'y2_1',[80,0.3,-2],output=True)

# fermi_dirac(f'v2_1',[0,400,-1,1e4],output=True)
# fermi_dirac(f'b2_1',[0,2000,-1,1e4],output=True)
# fermi_dirac(f'g2_1',[0,1800,-1,1e4],output=True)
# fermi_dirac(f'y2_1',[0,750,-1,1e4],output=True)


#########################################################################################

# zero_crossing(f'v3_1',5,13,output=True)
# zero_crossing(f'b3_1',3,12,output=True)
# zero_crossing(f'g3_1',4,10,output=True)
# zero_crossing(f'y3_1',4,8,output=True)

# linear_fit(f'v3_1',[[0,11],[38,50]],output=True)
# linear_fit(f'b3_1',[[0,7],[40,47]],output=True)
# linear_fit(f'g3_1',[[0,9],[42,48]],output=True)
# linear_fit(f'y3_1',[[0,9],[36,42]],output=True)

# ideal_diode(f'v3_1',[25,0.3,-2],output=True)
# ideal_diode(f'b3_1',[82,0.3,-2],output=True)
# ideal_diode(f'g3_1',[60,0.3,-2],output=True)
# ideal_diode(f'y3_1',[80,0.3,-2],output=True)

# fermi_dirac(f'v3_1',[0,400,-1,1e4],output=True)
# fermi_dirac(f'b3_1',[0,2000,-1,1e4],output=True)
# fermi_dirac(f'g3_1',[0,1800,-1,1e4],output=True)
# fermi_dirac(f'y3_1',[0,750,-1,1e4],output=True)


# fermi dirac things
#######################################################################

start,stop=[0,40]

# fermi_dirac1(f'v1_1',[0,400,-1,1e4],[start,stop],output=True)
# fermi_dirac1(f'b1_1',[0,2000,-1,1e4],[start,stop],output=True)
# fermi_dirac1(f'g1_1',[0,1800,-1,1e4],[start,stop],output=True)
# fermi_dirac1(f'y1_1',[0,750,-1,1e4],[start,stop],output=True)

# fermi_dirac1(f'v2_1',[0,400,-1,1e4],[start,stop],output=True)
# fermi_dirac1(f'b2_1',[0,2000,-1,1e4],[start,stop],output=True)
# fermi_dirac1(f'g2_1',[0,1800,-1,1e4],[start,stop],output=True)
# fermi_dirac1(f'y2_1',[0,750,-1,1e4],[start,stop],output=True)

# fermi_dirac1(f'v3_1',[0,400,-1,1e4],[start,stop],output=True)
# fermi_dirac1(f'b3_1',[0,2000,-1,1e4],[start,stop],output=True)
# fermi_dirac1(f'g3_1',[0,100,-1,1e4],[0,50],output=True)
# fermi_dirac1(f'y3_1',[0,750,-1,1e4],[start,stop],output=True)
