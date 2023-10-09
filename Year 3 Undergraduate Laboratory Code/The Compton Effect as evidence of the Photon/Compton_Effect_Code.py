'''
Data Analysis code for Compton Scattering Experiment

Code uses calibration data to correct data
code calculates differential cross-section from experimental data for comparison with theoretical results
code produces plots comparing experimental with theoretical data with error propagation included
'''


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

AM241_PEAK_LOC = 59.54
CS137_PEAK_LOC = 661.7
path = '/Users/Bryn_Lloyd/Documents/Uni/Year 3/Labs/Compton Effect/Data/'

def load_data(filename):
    data = pd.read_csv(f'{path}{filename}.csv', delimiter=',')
    # data.rename(columns={'E_1 / keV': 'E'},inplace=True)
    # data.rename(columns={'No Al1': 'No_Al'},inplace=True)
    data.drop(len(data)-1,inplace=True) # final row contains a random spike so just remove it

    if 'Cs' in filename:
        # print('Caesium')
        data['No_Al'] = (data['No Al1'] + data['No Al2'] + data['No Al3'])/3
        data['signal'] = data['Al']-data['No_Al']
        return calibrate_Cs(data)
    else:
        # print('Americium')
        data['signal'] = data['Al']-data['No_Al']
        return calibrate_Am(data)



def calibrate_Cs(data):
    scale = 0.399
    scale_error = 0.004
    shift = -12
    shift_error = 1
    data.rename(columns={'E_1 / keV': 'E_uncalibrated'},inplace=True)
    data['E'] = data['E_uncalibrated']*scale + shift
    data['E_error'] = abs( np.sqrt( (scale_error*data['E'])**2 + (0*scale)**2 + (shift_error)**2))    
    return data

# gradient: 0.39883055630366226 ± 0.004861266534202355
# intercept: -12.327894635809786 ± 1.1945632314196613


def calibrate_Am(data):
    scale = 0.0398
    scale_error = 0.0004
    shift = -2.4
    shift_error = 0.5
    data.rename(columns={'E_1 / keV': 'E_uncalibrated'},inplace=True)
    data['E'] = data['E_uncalibrated']*scale + shift
    data['E_error'] = abs( np.sqrt( (scale_error*data['E'])**2 + (0*scale)**2 + (shift_error)**2))    

    return data


# gradient: 0.039788264237061784 ± 0.0003521453922729296
# intercept: -2.3667579797534692 ± 0.49448421692032735




def remove_negative(data):
    for i in range(len(data)):
        if data['E'][i] < 0:
            data=data.drop(i)
    data=data.reset_index()
    return data


def Gaussian(x,peak,mu,sigma):
    return peak*np.exp(-(x-mu)**2/(2*sigma**2))

def chi_squared(expected,observed):
    sum=0
    for i in range(len(expected)):
        sum += ((observed[i]-expected[i])**2)/expected[i]
    return sum



def get_peak(filename,guess,output=False):
    data = load_data(filename)
    data = remove_negative(data)

    peak_xvals,peak_yvals=[],[]
    start,stop = guess[1]-2*guess[2],guess[1]+2*guess[2]
    for i in range(len(data['signal'])):
        if data['E'][i] > start and data['E'][i] < stop:
            peak_xvals.append(data['E'][i])
            peak_yvals.append(data['signal'][i])

    y_errors=[]
    for x in peak_xvals:
        y_errors.append(efficiency_error(x))
    y_errors=np.array(y_errors)
    initial_guess = [guess[0],guess[1],guess[2]]
    cfit = curve_fit(Gaussian,peak_xvals,peak_yvals,initial_guess,sigma=y_errors,absolute_sigma=True)
    if output == True:
        # plt.title('Combined background data')
        # plt.plot(data['E'],data['No Al1'],alpha=0.2)
        # plt.plot(data['E'],data['No Al2'],alpha=0.2)
        # plt.plot(data['E'],data['No Al3'],alpha=0.2)
        # plt.plot(data['E'],data['No_Al'])
        # plt.grid()
        # plt.xlabel('Energy / keV')
        # plt.ylabel('Counts')
        # plt.show()

        
        # plt.title('Signal')
        plt.plot(data['E'],data['signal'])
        x_peak_ls = np.linspace(cfit[0][1]-3*cfit[0][2],cfit[0][1]+3*cfit[0][2],100)
        plt.plot(x_peak_ls,Gaussian(x_peak_ls,*cfit[0]),color='C1')
        plt.grid()
        plt.xlabel('Energy / keV')
        plt.ylabel('Counts')
        plt.ylim(0,np.max(data['signal']))
        plt.title('Signal')
        plt.savefig('/Users/Bryn_Lloyd/Documents/Uni/Year 3/Labs/Compton Effect/Plots/Plot1',dpi=400)

        plt.show()

        # plt.title('Just peak')
        # plt.plot(peak_xvals,peak_yvals)
        # plt.plot(x_peak_ls,Gaussian(x_peak_ls,*cfit[0]),color='C1')
        # plt.grid()
        # plt.show()


        # n=4     # merges every n bins together - just to make the nice graph
        # plt.title(f'signal with 1/{n} as many bins')   
        # merged_bins=[]
        # for i in range(len(data)//n):
        #     merge=0
        #     for j in range(n):
        #         merge += data['signal'][n*i+j]
        #     merged_bins.append(merge/n)
        # plt.plot(data['E'][:len(data)-n:n],merged_bins) 
        # plt.plot(x_peak_ls,Gaussian(x_peak_ls,*cfit[0]),color='C1')
        # plt.grid()
        # plt.show()   


        window = 80
        background_rolling_average = data['No_Al'].rolling(window).mean()
        rolling_average = data['Al'].rolling(window).mean()
        signal_rolling_average = data['signal'].rolling(window).mean()

        # fig,ax = plt.subplots(3,sharex=True)
        # ax[0].set_title('Rolling Averages')
        # ax[1].set_ylabel('Counts')
        # ax[2].set_xlabel('Energy / keV')
        # for i in range(3):
        #     ax[i].grid()
        # ax[0].plot(data['E'],data['No_Al'],color='k',label='Background')
        # ax[0].plot(data['E'],background_rolling_average,color='C1')
        # ax[1].plot(data['E'],data['Al'],color='g',label='Background + Signal')
        # ax[1].plot(data['E'],rolling_average,color='C1')
        # ax[2].plot(data['E'],data['signal'],label='Signal')
        # ax[2].plot(data['E'],signal_rolling_average,color='C1')
        # fig.legend()
        # plt.show()

        plt.plot(data['E'],background_rolling_average,label='Background')
        plt.plot(data['E'],rolling_average,label='Background + Signal')
        plt.plot(data['E'],signal_rolling_average,label = 'Signal')
        plt.grid()
        plt.xlabel('Energy / keV')
        plt.ylabel('Counts')
        plt.legend()
        plt.title('Rolling average plot')
        plt.savefig('/Users/Bryn_Lloyd/Documents/Uni/Year 3/Labs/Compton Effect/Plots/rollingavg',dpi=400)


        plt.show()
        

        print(f'''Peak of {cfit[0][0]} ± {np.sqrt(cfit[1][0][0])}
        Peak at E = {cfit[0][1]} ± {np.sqrt(cfit[1][1][1])} keV
        With width {cfit[0][2]} ± {np.sqrt(cfit[1][2][2])}''')


    peak_pos = cfit[0][1]
    for i in range(len(data)):
        if data['E'][i] > peak_pos-0.2 and data['E'][i] < peak_pos+0.2:
            x = i
    peak_error = np.sqrt(abs(cfit[1][1][1])) + (data['E_error'][x])
    # print('relative_errors')
    # print(peak_error/peak_pos)
    # print()

    return peak_pos,peak_error,cfit

def compton_E(angle,initial_E):
    m_ec2=9.11e-31*(3e8)**2*1e-3/(1.6e-19)
    # initial_E = 661.7
    E_f=(initial_E)/(1+((initial_E)/(m_ec2))*(1-np.cos(angle)))
    return E_f

def compton_E_me(angle,m_e):
    c2=(3e8)**2*1e-3/(1.6e-19)
    initial_E = 661.7
    actualM_e = 9.11e-31
    E_f=(initial_E)/(1+((initial_E)/((m_e*actualM_e)*c2))*(1-np.cos(angle)))
    return E_f


def compton_scattering(array,x='theta',output=False):
    global theta                        # used in multiple functions so just make it a global variable
    theta = np.linspace(-np.pi,np.pi,100)      
    #  predictd y-values
    if 'Cs' in array[0][0]:
        y_vals = compton_E(theta,CS137_PEAK_LOC)
    else:
        y_vals = compton_E(theta,AM241_PEAK_LOC)

    # actual y-values
    peaks,peak_errors=[],[]
    for i in range(len(array)):
        peaks.append(get_peak(*array[i])[0])
        peak_errors.append(get_peak(*array[i])[1])




    x_vals = rad(array[::,2])
    R1 = 0.0042
    R2 = 0.02016/2
    R3 = 0.01971
    L1 = 0.1
    L2 = 0.169
    global angle_error1
    angle_error1 = np.arctan(R3/(2*L2)) + np.arctan(R1/(2*L1))

    angle_error2 = []
    for i in x_vals:
        a = rad(90+i/2)
        L11 = L1**2 + R2**2 -2*L1*R2*np.cos(a)
        L21 = L2**2 + R2**2 -2*L2*R2*np.cos(a)
        angle_error2.append(np.arcsin((R2*np.sin(a)/(L11))) + np.arcsin((R2*np.sin(a)/(L21))))
        
    if output == True:
        if x == 'theta':
            x_vals = rad(array[::,2])
            plt.plot(theta,y_vals,'--',label='Theoretical Values')
        elif x == 'cos':
            x_vals = np.cos(rad(array[::,2]))
            plt.plot(np.cos(theta),y_vals,label='Theoretical Values')


        # plt.plot(x_vals,peaks,'x',label='Data')
        plt.errorbar(x_vals,peaks,yerr = peak_errors,xerr = angle_error1,fmt='x',capsize=5,label='Compton Peaks')

        compton_fit = curve_fit(compton_E,x_vals,peaks)
        # print(compton_fit)


        if x == 'theta':
            # plt.plot(theta,compton_E(theta,*compton_fit[0]),':',label='Fit to data')
            ticks = np.arange(0,2*np.pi,np.pi/4)
            labels = [0,'$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$','$5\pi/4$','$3\pi/2$','$7\pi/4$']
            plt.xticks(ticks=ticks,labels=labels)
            plt.xlabel('Scattering Angle')
            plt.ylabel('Energy / keV')
        elif x == 'cos':
            # plt.plot(np.cos(theta),compton_E(theta,*compton_fit[0]),':',label='Fit to data')
            plt.xlabel('Cos(angle)')
            plt.ylabel('Energy / keV')

        # can curvefit with electron mass as independent variable using compton_E_me
        # m_e_fit = curve_fit(compton_E_me,x_vals,peaks)
        # print(f'predicted electron mass: {m_e_fit[0][0]} x electron mass')

        plt.grid()
        plt.legend()
        plt.xlim(0,np.pi)
        plt.savefig('/Users/Bryn_Lloyd/Documents/Uni/Year 3/Labs/Compton Effect/Plots/final_compton_graph',dpi=400)
        plt.show()

        exp=[]
        for x in x_vals:
            exp.append(compton_E(x,CS137_PEAK_LOC))
        obs = peaks
        chi_sq = chi_squared(exp,obs)
        print(f'compton effect chi-squared = {chi_sq}')



    return peaks

def rad(deg):
    return deg * np.pi / 180
  


def differential_cross_section(E0,array,output=True):
    def diff_cross_sec(E0,E1,angle):
        ratio=E1/E0
        cross_sec=1/2*(((7.94e-30))*(ratio**2)*(ratio+1/ratio -(np.sin(angle))**2))
        return cross_sec 
    def dcs_err(E_0,E1,E1_err,angle,angle_err):
        R=E1/E0
        R_err = E1_err / E0
        err = np.sqrt( R_err**2*(0.5*7.94e-30*(3*R**2 + 1 - 2*R*np.sin(angle)**2))**2 + angle_err**2*(7.94e-30*R**2*np.sin(angle)*np.cos(angle))**2)
        return err




    energies = compton_scattering(array)
    angles = rad(array[::,2])

    theta=np.linspace(-np.pi,np.pi,100)
    Energies=compton_E(theta,E0)
    cross_sections,cross_section_errors,expected=[],[],[]

    peak_errors=[]
    for i in range(len(array)):
        peak_errors.append(get_peak(*array[i])[1])
    peak_errors = np.array(peak_errors)


    for i in range(len(energies)):
        cross_sections.append(diff_cross_sec(E0,energies[i],angles[i]))
        cross_section_errors.append(dcs_err(E0,energies[i],peak_errors[i],angles[i],angle_error1))
    for i in range(100):
        expected.append(diff_cross_sec(E0,Energies[i],theta[i]))


    print(cross_section_errors)


    if output == True:
        # fig = plt.figure()
        # ax = plt.axes(polar=True)
        # if 'Cs' in array[0]:
        #     ax.plot(angles,cross_sections,'x',label='Cs')
        # else:
        #     ax.plot(angles,cross_sections,'x',label='Am')
        # ax.plot(theta,expected)

        # ticks = np.arange(0,2*np.pi,np.pi/4)
        # labels = [0,'$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$','$5\pi/4$','$3\pi/2$','$7\pi/4$']
        # plt.xticks(ticks=ticks,labels=labels)
        # plt.legend(loc='upper right')

        # ax.errorbar(angles,cross_sections,yerr=cross_section_errors,xerr=angle_error1,fmt='x',capsize=5)

        plt.show()
        # plt.title('Klein Nishina cross-sections')
        # plt.plot(angles,cross_sections,'x')
        plt.errorbar(angles,cross_sections,yerr=cross_section_errors,xerr=angle_error1,fmt='x',capsize=5,ecolor='C1')
        plt.grid()
        plt.xticks(ticks=[0,np.pi/4,np.pi/2,3*np.pi/4,np.pi],labels=[0,'$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'])
        plt.xlabel('Scattering Angle')
        plt.ylabel('$d\sigma / d\Omega$ / $m^2$')


        plt.show()
    
    return cross_sections,cross_section_errors,expected

def thomson(theta):
    q = 1.6e-19
    e0 = 8.85e-12
    m = 9.11e-31
    c = 3e8
    coeff = (q**2 / (4*np.pi*e0*m*c**2))**2
    return coeff * (1 + np.cos(theta)**2)/2
    
# [filename, initial guess array for peak position, scattering angle(degrees)]
Cs_10 = ['Cs_10deg',[15,700,70],10]    # [peak height, peak position, peak width]
Cs_20 = ['Cs_20deg',[15,600,50],20]
Cs_30 = ['Cs_30deg',[15,550,70],30]
Cs_40 = ['Cs_40deg',[15,500,70],40]
Cs_50 = ['Cs_50deg',[15,525,70],50]
Cs_60 = ['Cs_60deg',[15,400,70],60]
Cs_70 = ['Cs_70deg',[5,380,20],70]
Cs_80 = ['Cs_80deg',[5,330,20],80]
Cs_90 = ['Cs_90deg',[15,360,70],90]
Cs_100 = ['Cs_100deg',[15,360,70],100]
Cs_110 = ['Cs_110deg',[7,260,50],110]
Cs_120 = ['Cs_120deg',[7,240,30],120]
Cs_130 = ['Cs_130deg',[6,220,30],130]
Cs_array=np.array([Cs_10,Cs_20,Cs_30,Cs_40,Cs_50,Cs_60,Cs_70,Cs_80,Cs_90,Cs_100,Cs_110,Cs_120,Cs_130],dtype='object')


Am_15 = ['Am_15deg',[30,50,20],15]
Am_30 = ['Am_30deg',[30,50,20],30]
Am_45 = ['Am_45deg',[30,50,20],45]
Am_60 = ['Am_60deg',[30,50,20],60]
Am_90 = ['Am_90deg',[30,50,20],90]
Am_120 = ['Am_120deg',[30,50,20],120]
Am_150 = ['Am_150deg',[30,50,20],150]
Am_array=np.array([Am_15,Am_30,Am_45,Am_60,Am_90,Am_120,Am_150],dtype='object')


def flux(array,aquisition_time,area=np.pi*(19.62/2)**2):

    data = load_data(array[0])
    data = data.rename(columns={'signal': 'N'})

    cfit = get_peak(array[0],array[1])[2][0]
    # plt.plot(data['E'],data['N'])
    # plt.plot(data['E'],Gaussian(data['E'],*cfit))

    start = cfit[1] - 3*cfit[2]
    stop = cfit[1] + 3*cfit[2]

    for i in range(len(data)):
        if data['E'][i] > start-0.2 and data['E'][i] < start+0.2:
            start_pos = i
    for i in range(len(data)):
        if data['E'][i] > stop-0.2 and data['E'][i] < stop+0.2:
            stop_pos = i

    count = 0
    for i in range(start_pos,stop_pos):
        count += data['N'][i]



    flux_error = efficiency_error(cfit[1])/(aquisition_time*area)

    return count/(aquisition_time*area), flux_error






def incident_flux():
    incidentflux_data = pd.read_csv(f'{path}Cs_incidentflux.csv', delimiter=',')
    # incidentflux_data = incidentflux_data.rename(columns={'E_1 / keV': 'E'})
    incidentflux_data = incidentflux_data.drop(len(incidentflux_data)-1) # final row contains a random spike so just remove it
    incidentflux_data = calibrate_Cs(incidentflux_data)

    peak = np.max(incidentflux_data['N'])
    peak_loc = np.float64(incidentflux_data.loc[incidentflux_data['N']==peak]['E'])

    initial_guess = [peak,peak_loc,1]
    cfit = curve_fit(Gaussian,incidentflux_data['E'],incidentflux_data['N'], initial_guess)
    initial_guess = cfit[0]
    start = cfit[0][1]- 3*cfit[0][2]
    stop = cfit[0][1]+ 3*cfit[0][2]
    peak_xvals,peak_yvals=[],[]
    for i in range(len(incidentflux_data['E'])):
        if i > start and i < stop:
            peak_xvals.append(incidentflux_data['E'][i])
            peak_yvals.append(incidentflux_data['N'][i])
    cfit1 = curve_fit(Gaussian,peak_xvals,peak_yvals,initial_guess)



    # plt.plot(incidentflux_data['E'],incidentflux_data['N'])
    # plt.plot(peak_xvals,Gaussian(peak_xvals,*cfit1[0]))
    # plt.show()

    aquisition_time = 200
    area=np.pi*(19.62/2)**2                                                                                               ###factor of 1000 out!!!!!##########
    count = 0
    # start_pos = int(incidentflux_data.loc[incidentflux_data['E']==int(start)]['E'])
    # stop_pos = int(incidentflux_data.loc[incidentflux_data['E']==int(stop)]['E'] +1)


    for i in range(len(incidentflux_data)):
        if incidentflux_data['E'][i] > start-0.2 and incidentflux_data['E'][i] < start+0.2:
            start_pos = i
    for i in range(len(incidentflux_data)):
        if incidentflux_data['E'][i] > stop-0.2 and incidentflux_data['E'][i] < stop+0.2:
            stop_pos = i




    for i in range(start_pos,stop_pos):
        count += incidentflux_data['N'][i]
    F_inc = count/(aquisition_time*area)

    # print(f'incident flux = {F_inc}')

    F_inc_error = 1000 * efficiency_error(cfit[0][1])/(aquisition_time*area)       #times by 100.....????? dodgy
    print('in_errors')
    print(F_inc_error/F_inc)
    print()

    return F_inc,F_inc_error


def efficiency_error(x):
    yn_val = 1.04158576e+02 - 2.02938024e-02*x - 6.87104639e-05*x**2 + 7.66234556e-08*x**3 - 2.93776613e-11*x**4
    error=(1-yn_val/100)*x
    return error






def gcs():
    n = 2710*13/(26*1.67e-27)
    t = 0.02016
    omega=4*np.pi*(0.01971/2)**2/((0.169)**2)

    F_inc,F_inc_error = incident_flux()
    aquisition_time = [200,400,400,400,400,600,600,600,600,600,600,600,600]
    F_out,F_out_error,dcs,dcs_error = [],[],[],[]
    for i,array in enumerate(Cs_array):
        F_out.append(flux(array,aquisition_time[i])[0])
        F_out_error.append(flux(array,aquisition_time[i])[1])
        dcs.append(1/(n*t*omega) * F_out[i]/F_inc)
        dcs_error.append(dcs[i]*np.sqrt((F_inc_error/F_inc)**2 + (F_out_error[i]/F_out[i])**2))

    print('out_errors')
    print(np.array(F_out_error)/np.array(F_out))
    print()

    dcs_error_lower=[]
    for i in range(len(Cs_array)):
        dcs_error_lower.append(dcs[i]*np.sqrt((0.15)**2 + (0.20)**2))

    asymmetric_error = [dcs_error_lower,dcs_error]    


    KN_cross_sections,KN_errors,KN_theoretical = differential_cross_section(CS137_PEAK_LOC,Cs_array,output=False)


    plt.xlabel('Scattering Angle')
    plt.ylabel('$d\sigma / d\Omega$ / $m^2$')
    # plt.plot(rad(Cs_array[::,2]),KN_cross_sections,'x',label='Klein-Nishina cross-sections')
    # plt.errorbar(rad(Cs_array[::,2]),KN_cross_sections,yerr=KN_errors,xerr=angle_error1,fmt='x',label='Klein-Nishina cross-sections',capsize=3)
    ls = np.linspace(0,np.pi,50)
    KN_theoretical = np.array(KN_theoretical)[50::]
    plt.plot(ls,KN_theoretical,'--',label='Klein-Nishina differential cross-section')
    # plt.plot(rad(Cs_array[::,2]),dcs,'x',label='General cross-sections')
    # plt.errorbar(rad(Cs_array[::,2]),dcs,yerr=dcs_error,xerr=angle_error1,fmt='x',label='Experimental differential cross-sections',capsize=3)
    plt.errorbar(rad(Cs_array[::,2]),dcs,yerr=asymmetric_error,xerr=angle_error1,fmt='x',label='Experimental differential cross-sections',capsize=3)
    plt.grid()
    plt.xticks(ticks=[0,np.pi/4,np.pi/2,3*np.pi/4,np.pi],labels=[0,'$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'])
    thomson_vals = []
    ls = np.linspace(0,np.pi,100)
    for i in ls:
        thomson_vals.append(thomson(i))
    plt.plot(ls,thomson_vals,':',color='k',label='Thomson differential cross-section')    
    plt.legend()

    plt.show()


    # plt.polar(rad(Cs_array[::,2]),KN_cross_sections,'x',label='Klein-Nishina cross-sections')
    # plt.polar(rad(Cs_array[::,2]),dcs,'x',label='General cross-sections')
    # ticks = np.arange(0,2*np.pi,np.pi/4)
    # labels = [0,'$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$','$5\pi/4$','$3\pi/2$','$7\pi/4$']
    # plt.xticks(ticks=ticks,labels=labels)
    # plt.show()


    


    # fig = plt.figure()
    # ax1 = plt.subplot(121)
    # ax2 = plt.subplot(122, projection='polar')
    # ax1.plot(rad(Cs_array[::,2]),KN_cross_sections,'x',label='Klein-Nishina cross-sections')
    # ax1.plot(rad(Cs_array[::,2]),dcs,'x',label='General cross-sections')
    # ax1.set_xticks(ticks=[0,np.pi/4,np.pi/2,3*np.pi/4],labels=[0,'$\pi/4$','$\pi/2$','$3\pi/4$'])

    # ax1.grid()


    # ax2.plot(rad(Cs_array[::,2]),KN_cross_sections,'x')
    # ax2.plot(rad(Cs_array[::,2]),dcs,'x')
    # ticks = np.arange(0,2*np.pi,np.pi/4)
    # labels = [0,'$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$','$5\pi/4$','$3\pi/2$','$7\pi/4$']
    # ax2.set_xticks(ticks=ticks,labels=labels)

    # fig.legend()






# get_peak(Cs_40[0],Cs_40[1],output=True)
compton_scattering(Cs_array,output=True)
# differential_cross_section(CS137_PEAK_LOC,Cs_array)
# differential_cross_section(AM241_PEAK_LOC,Am_array)
# gcs()





# plt.savefig('/Users/Bryn_Lloyd/Documents/Uni/Year 3/Labs/Compton Effect/Plots/all_dcs_compared',dpi=400)
