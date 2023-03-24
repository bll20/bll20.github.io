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

ls = np.linspace(-5,10,100)

def plots(data):
    plt.plot(data['V'],data['I'],'x')
    plt.grid()
    plt.xlabel('Voltage / V')
    plt.ylabel('Current / nA')


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
    



def zero_crossing(filename,interval,output=False):
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
    x = data['V'][zero+1:upper:]                     #[lower:upper] or [zero:upper]
    y = data['I'][zero+1:upper:]
    

    fit,cov = np.polyfit(x,y,1,cov=True)
    pfit = np.poly1d(fit)

    V_stop = -fit[1]/fit[0]
    err_m = np.sqrt(cov[0][0])
    err_b = np.sqrt(cov[1][1])
    err_Vstop = abs(V_stop*np.sqrt((err_b/fit[1])**2 + (err_m/fit[0])**2))


    if output == True:
        fig,ax = plt.subplots(2)
        ax[0].set_title('Zero-Crossing Method')
        for i in range(2):
            ax[i].grid()
        ax[0].plot(data['V'],data['I'],'x')
        ax[0].plot(ls,pfit(ls))
        ax[0].plot(x,y,'x',color='red',label='fit points')
        ax[0].errorbar(V_stop,0,xerr=err_Vstop,capsize=3,color='k',fmt='o')
        ax[1].plot(x,y,'x')
        ax[1].plot(x,pfit(x))
        ax[1].errorbar(V_stop,0,xerr=err_Vstop,capsize=3,color='k',fmt='o',label=f'$V_s$ = {V_stop} ± {err_Vstop} V')
        fig.legend(loc='lower center')
        sexy_subplot(fig,ax)
        plt.show()

    return V_stop, err_Vstop


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

    V_stop = (fit2[1]-fit1[1])/(fit1[0]-fit2[0])
    I_stop = pfit1(V_stop)


    err_m1 = np.sqrt(cov1[0][0])
    err_c1 = np.sqrt(cov1[1][1])
    err_m2 = np.sqrt(cov2[0][0])
    err_c2 = np.sqrt(cov2[1][1])
    err_Vstop = np.sqrt((err_c1/(fit1[0]-fit2[0]))**2 + (-err_c2/(fit1[0]-fit2[0]))**2 + (err_m1 * (fit2[1]-fit1[1])/(fit1[0]-fit2[0])**2)**2 + (err_m2 * (fit1[1]-fit2[1])/(fit1[0]-fit2[0])**2)**2)

    if output==True:
        plt.title('Linear Fits')
        plt.plot(data['V'][:stop2],data['I'][:stop2],'x')
        ls1 = np.linspace(data['V'][start1],data['V'][start2+1])
        plt.plot(ls1,pfit1(ls1))
        ls2 = np.linspace(data['V'][stop1-1],data['V'][stop2])
        plt.errorbar(V_stop,I_stop,xerr=err_Vstop,capsize=3,color='k',fmt='o',label=f'$V_s$ = {V_stop} ± {err_Vstop} V')
        plt.plot(ls2,pfit2(ls2))
        plt.grid()
        plt.xlabel('Voltage / V')
        plt.ylabel('Current / nA')
        plt.legend()
        sexy_plot()
        plt.show()

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
        plt.show()
        
        print(cfit[0])

    return V_stop,err_Vstop


FREQUENCY_V = 3e8 / 405e-9
FREQUENCY_B = 3e8 / 436e-9
FREQUENCY_G = 3e8 / 546e-9
FREQUENCY_Y = 3e8 / 578e-9
FREQUENCIES = [FREQUENCY_V,FREQUENCY_B,FREQUENCY_G,FREQUENCY_Y]


def finding_h(function):
    if function == 'zero_crossing':
        V_stop_v,err_Vstop_v = zero_crossing(f'v1_1',10)
        V_stop_b,err_Vstop_b = zero_crossing(f'b1_1',10)
        V_stop_g,err_Vstop_g = zero_crossing(f'g1_1',10)
        V_stop_y,err_Vstop_y = zero_crossing(f'y1_1',10)

        V_stop_v_2,err_Vstop_v_2 = zero_crossing(f'v2_1',10)
        V_stop_b_2,err_Vstop_b_2 = zero_crossing(f'b2_1',10)
        V_stop_g_2,err_Vstop_g_2 = zero_crossing(f'g2_1',10)
        V_stop_y_2,err_Vstop_y_2 = zero_crossing(f'y2_1',10)

        V_stop_v_3,err_Vstop_v_3 = zero_crossing(f'v3_1',6)
        V_stop_b_3,err_Vstop_b_3 = zero_crossing(f'b3_1',10)
        V_stop_g_3,err_Vstop_g_3 = zero_crossing(f'g3_1',10)
        V_stop_y_3,err_Vstop_y_3 = zero_crossing(f'y3_1',20)
    elif function == 'linear_fit':
        V_stop_v,err_Vstop_v = linear_fit(f'v1_1',[[0,15],[25,40]])
        V_stop_b,err_Vstop_b = linear_fit(f'b1_1',[[0,15],[25,40]])
        V_stop_g,err_Vstop_g = linear_fit(f'g1_1',[[0,15],[25,40]])
        V_stop_y,err_Vstop_y = linear_fit(f'y1_1',[[0,15],[25,40]])

        V_stop_v_2,err_Vstop_v_2 = linear_fit(f'v2_1',[[0,15],[25,40]])
        V_stop_b_2,err_Vstop_b_2 = linear_fit(f'b2_1',[[0,15],[25,39]])
        V_stop_g_2,err_Vstop_g_2 = linear_fit(f'g2_1',[[0,15],[25,38]])
        V_stop_y_2,err_Vstop_y_2 = linear_fit(f'y2_1',[[0,15],[30,35]])

        V_stop_v_3,err_Vstop_v_3 = linear_fit(f'v3_1',[[0,15],[35,45]])
        V_stop_b_3,err_Vstop_b_3 = linear_fit(f'b3_1',[[0,15],[35,50]])
        V_stop_g_3,err_Vstop_g_3 = linear_fit(f'g3_1',[[0,15],[35,45]])
        V_stop_y_3,err_Vstop_y_3 = linear_fit(f'y3_1',[[0,15],[30,40]])
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



    V_stops = np.array([abs(V_stop_v),abs(V_stop_b),abs(V_stop_g),abs(V_stop_y)])
    err_Vstops = np.array([abs(err_Vstop_v),abs(err_Vstop_b),abs(err_Vstop_g),abs(err_Vstop_y)])


    V_stops_2 = np.array([abs(V_stop_v_2),abs(V_stop_b_2),abs(V_stop_g_2),abs(V_stop_y_2)])
    err_Vstops_2 = np.array([abs(err_Vstop_v_2),abs(err_Vstop_b_2),abs(err_Vstop_g_2),abs(err_Vstop_y_2)])

    V_stops_3 = np.array([abs(V_stop_v_3),abs(V_stop_b_3),abs(V_stop_g_3),abs(V_stop_y_3)])
    err_Vstops_3 = np.array([abs(err_Vstop_v_3),abs(err_Vstop_b_3),abs(err_Vstop_g_3),abs(err_Vstop_y_3)])

    



    plt.errorbar(FREQUENCIES,V_stops_2,yerr=err_Vstops_2,fmt='x',capsize=3,color='k',label='light nd filter')
    plt.errorbar(FREQUENCIES,V_stops,yerr=err_Vstops,fmt='x',capsize=3,color='C0',label='no filter')
    plt.errorbar(FREQUENCIES,V_stops_3,yerr=err_Vstops_3,fmt='x',capsize=3,color='red',label='dark nd filter')
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
    plt.plot(ls,pfit1(ls),'--',color='C0')
    plt.plot(ls,pfit2(ls),'--',color='k')
    plt.plot(ls,pfit3(ls),'--',color='red')
    # plt.plot(ls,pfit_tot(ls),color='C1')

    plt.grid()
    plt.ylabel('Cut-Off Voltage / V')
    plt.xlabel('Frequency / Hz')
    plt.title(f'Finding Plancks Constant using {function} Model')
    sexy_plot()
    plt.show()


    h = fit_tot[0]*1.6e-19
    h_err = np.sqrt(cov_tot[0][0]) * 1.6e-19
    print(f'''Planck's Constant: {h} ± {h_err} Js''')
    print()

    h1 = fit1[0]*1.6e-19
    h1_err = np.sqrt(cov1[0][0]) * 1.6e-19
    print(f'''Planck's Constant (no nd filter): {h1} ± {h1_err} Js''')
    h2 = fit2[0]*1.6e-19
    h2_err = np.sqrt(cov2[0][0]) * 1.6e-19
    print(f'''Planck's Constant (light nd filter): {h2} ± {h2_err} Js''')
    h3 = fit3[0]*1.6e-19
    h3_err = np.sqrt(cov3[0][0]) * 1.6e-19
    print(f'''Planck's Constant (dark nd filter): {h3} ± {h3_err} Js''')
    
    h_avg = (h1+h2+h3)/3
    h_avg_err = np.sqrt(h1_err**2 + h2_err**2+ h3_err**2)/3
    print(f'''Averaged Planck's Constant: {h_avg} ± {h_avg_err} Js''')




###################################################################################

# finding_h('zero_crossing')       # changed to [zero+1:upper]
# finding_h('linear_fit') 
# finding_h('ideal_diode')
finding_h('fermi_dirac')


#################################################################################


zero_crossing(f'v1_1',10,output=True)
# zero_crossing(f'b1_1',10,output=True)
# zero_crossing(f'g1_1',10,output=True)
# zero_crossing(f'y1_1',6,output=True)

# linear_fit(f'v1_1',[[0,15],[25,40]],output=True)
# linear_fit(f'b1_1',[[0,15],[25,40]],output=True)
# linear_fit(f'g1_1',[[0,15],[25,40]],output=True)
# linear_fit(f'y1_1',[[0,15],[25,40]],output=True)

# ideal_diode(f'v1_1',[25,0.3,-2],output=True)
# ideal_diode(f'b1_1',[82,0.3,-2],output=True)
# ideal_diode(f'g1_1',[60,0.3,-2],output=True)
# ideal_diode(f'y1_1',[80,0.3,-2],output=True)

fermi_dirac(f'v1_1',[0,400,-1,1e4],output=True)
# fermi_dirac(f'b1_1',[0,2000,-1,1e4],output=True)
# fermi_dirac(f'g1_1',[0,1800,-1,1e4],output=True)
# fermi_dirac(f'y1_1',[0,750,-1,1e4],output=True)



#######################################################################################

# zero_crossing(f'v2_1',10,output=True)
# zero_crossing(f'b2_1',10,output=True)
# zero_crossing(f'g2_1',10,output=True)
# zero_crossing(f'y2_1',10,output=True)

# linear_fit(f'v2_1',[[0,15],[25,40]],output=True)
# linear_fit(f'b2_1',[[0,15],[25,39]],output=True)
# linear_fit(f'g2_1',[[0,15],[25,38]],output=True)
# linear_fit(f'y2_1',[[0,15],[26,32]],output=True)

# ideal_diode(f'v2_1',[25,0.3,-2],output=True)
# ideal_diode(f'b2_1',[82,0.3,-2],output=True)
# ideal_diode(f'g2_1',[60,0.3,-2],output=True)
# ideal_diode(f'y2_1',[80,0.3,-2],output=True)

# fermi_dirac(f'v2_1',[0,400,-1,1e4],output=True)
# fermi_dirac(f'b2_1',[0,2000,-1,1e4],output=True)
# fermi_dirac(f'g2_1',[0,1800,-1,1e4],output=True)
# fermi_dirac(f'y2_1',[0,750,-1,1e4],output=True)


#########################################################################################

# zero_crossing(f'v3_1',6,output=True)
# zero_crossing(f'b3_1',10,output=True)
# zero_crossing(f'g3_1',13,output=True)
# zero_crossing(f'y3_1',10,output=True)

# linear_fit(f'v3_1',[[0,15],[35,45]],output=True)
# linear_fit(f'b3_1',[[0,15],[35,50]],output=True)
# linear_fit(f'g3_1',[[0,15],[35,45]],output=True)
# linear_fit(f'y3_1',[[0,15],[30,40]],output=True)

# ideal_diode(f'v3_1',[25,0.3,-2],output=True)
# ideal_diode(f'b3_1',[82,0.3,-2],output=True)
# ideal_diode(f'g3_1',[60,0.3,-2],output=True)
# ideal_diode(f'y3_1',[80,0.3,-2],output=True)

# fermi_dirac(f'v3_1',[0,400,-1,1e4],output=True)
# fermi_dirac(f'b3_1',[0,2000,-1,1e4],output=True)
# fermi_dirac(f'g3_1',[0,1800,-1,1e4],output=True)
# fermi_dirac(f'y3_1',[0,750,-1,1e4],output=True)