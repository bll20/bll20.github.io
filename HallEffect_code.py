import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# THICKNESS = 3e-6


def relative_quadrature(array):          # array should be of the form [err,val]
    rel_err=[]
    for i in range(len(array)):
        rel_err.append(array[0][i]/array[1][i])
    rel_err = np.array(rel_err)
    return np.sqrt(np.sum(rel_err**2))


def quadratic_error(x,a_err,b_err,c_err):
    return x**4 * a_err**2 + x**2 * b_err**2 + c_err**2



def calculate_resistivity(filename,THICKNESS):        # returns a polyfit, and an error in each coefficient
    path = '/Users/Bryn_Lloyd/Documents/Uni/Year 3/Labs/Hall Effect/Data/Resistivity/'
    data = pd.read_csv(f'{path}{filename}.csv', delimiter=',')
    vdp_data = pd.read_csv(f'{path}van_der_pauw_correction_factor.csv', delimiter=',')


    for i in ['1243','2134','4312','3421','1423','4132','2314','3241']:
        resistance = []
        for j in range(len(data)):
            resistance.append(data['V'+i][j] / data['I'][j])
        data['R'+i] = resistance

    R_A = []
    R_B = []
    for j in range(len(data)):
        R_A.append((data['R1243'][j]+data['R2134'][j]+data['R4312'][j]+data['R3421'][j])/4)
        R_B.append((data['R1423'][j]+data['R4132'][j]+data['R2314'][j]+data['R3241'][j])/4)
    data['R_A'] = R_A
    data['R_B'] = R_B

    vdp_fit = np.polyfit(vdp_data['x'],vdp_data['y'],20)
    vdp_pfit = np.poly1d(vdp_fit)

    R_sh = []
    resistivity = []
    for j in range(len(data)):
        f = vdp_pfit(data['R_A'][j] / data['R_B'][j])
        R_sh.append((np.pi * (data['R_A'][j] + data['R_B'][j]))/(2*np.log(2)) * f)
        resistivity.append(abs(R_sh[j]*THICKNESS*100))
    data['R_sh'] = R_sh
    data['resistivity'] = resistivity


    for i in ['1243','2134','4312','3421','1423','4132','2314','3241']:
        R_err = []
        for j in range(len(data)):
            V_err = 0.06*data['V'+i][j]
            I_err = 0.001*data['I'][j]
            R_err.append(np.sqrt((V_err/data['V'+i][j])**2 + (I_err/data['I'][j])**2))
        data['R'+i+'_err'] = R_err   


    R_A_err = []
    R_B_err = []
    for j in range(len(data)):
        Aerror = [data['R1243_err'][j],data['R2134_err'][j],data['R4312_err'][j],data['R3421_err'][j]]
        Avalue = [data['R1243'][j],data['R2134'][j],data['R4312'][j],data['R3421'][j]] 
        R_A_err.append(relative_quadrature([Aerror,Avalue]))
        Berror = [data['R1423_err'][j],data['R4132_err'][j],data['R2314_err'][j],data['R3241_err'][j]]
        Bvalue = [data['R1423'][j],data['R4132'][j],data['R2314'][j],data['R3241'][j]] 
        R_B_err.append(relative_quadrature([Berror,Bvalue]))
    data['R_A_err'] = R_A_err
    data['R_B_err'] = R_B_err


    R_sh_err = []
    resistivity_err = []
    for j in range(len(data)):
        R_sh_err.append(np.sqrt((R_A_err[j])**2 + (R_B_err[j])**2))
        resistivity_err.append(abs(data['resistivity'][j] * R_sh_err[j]/data['R_sh'][j]))
    data['R_sh_err'] = R_sh_err
    data['resistivity_err'] = resistivity_err


    currents = set(data['I'])
    temp = []
    for i in currents:
        temp.append(i)
    currents=np.sort(temp)    
    fits,covs,=[],[]
    fits_err = []
    for i,current in enumerate(currents):
        x,y,y_err = [],[],[]
        for j in range(len(data)):
            if data['I'][j] == current:
                x.append(data['B'][j])
                y.append(data['resistivity'][j])
                y_err.append(data['resistivity_err'][j])
        fit,cov = np.polyfit(x,y,2,cov=True)
        fits.append(fit)
        covs.append(cov)
        fits_err.append([np.sqrt(cov[0][0]),np.sqrt(cov[1][1]),np.sqrt(cov[2][2])])



    fits = np.array(fits)
    average_fit = np.zeros(len(fits[0]))
    for i in range(len(average_fit)):
        average_fit[i] = np.mean(fits[:,i])


    fits_err = np.array(fits_err)
    average_fit_error = np.zeros(len(fits_err[0]))
    for i in range(len(average_fit_error)):
        temp = np.sum((fits_err[:,i]/fits[:,i])**2)
        average_fit_error[i] = average_fit[i] * np.sqrt(temp)

    return average_fit,average_fit_error


def calculate_hall(filename,THICKNESS):                # returns hall_coefficient, error
    path = '/Users/Bryn_Lloyd/Documents/Uni/Year 3/Labs/Hall Effect/Data/Hall_Voltage/vary_B/'
    data = pd.read_csv(f'{path}{filename}.csv', delimiter=',',skipfooter=1)


    y_1,y_2,y_3,y_4=[],[],[],[]
    for i in range(len(data)):
        y_1 = data['V_1']/data['I_1']
        y_2 = data['V_2']/data['I_2']
        y_3 = data['V_3']/data['I_3']
        y_4 = data['V_4']/data['I_4']
    data['y_1'],data['y_2'],data['y_3'],data['y_4'] = y_1,y_2,y_3,y_4

    fit_1,cov_1 = np.polyfit(data['B_1'],data['y_1'],1,cov=True)
    fit_2,cov_2 = np.polyfit(data['B_2'],data['y_2'],1,cov=True)
    fit_3,cov_3 = np.polyfit(data['B_3'],data['y_3'],1,cov=True)
    fit_4,cov_4 = np.polyfit(data['B_4'],data['y_4'],1,cov=True)

    combined_fit = THICKNESS * (fit_1+fit_2+fit_3+fit_4)/4

    for i in ['_1','_2','_3','_4']:
        y_err = []
        for j in range(len(data)):
            V_err = 0.0035*data['V'+i][j]
            I_err = 0.001*data['I'+i][j]
            y_err.append(abs(np.sqrt((V_err/data['V'+i][j])**2 + (I_err/data['I'+i][j])**2)))
        data['y'+i+'_err'] = y_err

    combined_gradient_err = combined_fit[0] * np.sqrt((np.sqrt(cov_1[0][0])/fit_1[0])**2 + (np.sqrt(cov_2[0][0])/fit_2[0])**2 + (np.sqrt(cov_3[0][0])/fit_3[0])**2 + (np.sqrt(cov_4[0][0])/fit_4[0])**2)
    return 1e6 * combined_fit[0], abs(combined_gradient_err)

def single_carrier(R_H,R_H_err,resistivity,resistivity_err):
    mobility = abs(R_H/resistivity)
    mobility_error = mobility * np.sqrt((R_H_err/R_H)**2+(resistivity_err/resistivity)**2)
    carrier_density = 1/(resistivity*mobility*1.6e-19)

    return mobility,mobility_error,carrier_density

def two_carrier(R_H,R_H_err,Res_0,Res_0_err,mobility_value=0):
    mobility = mobility_value - R_H/Res_0
    mobility_error = mobility * np.sqrt((R_H_err/R_H)**2+(Res_0_err/Res_0)**2)
    carrier_density = (mobility-mobility_value)/(1.6e-19 * R_H*(mobility-mobility_value))    
    return mobility,mobility_error,carrier_density

def two_carrier_2(R_H,R_H_err,Res_0,Res_0_err,b,b_err,mobility_value=0):
    # mobility = (np.sqrt(Res_0)*mobility_value + 2*np.sqrt(b))/np.sqrt(Res_0)

    mobility = mobility_value - 4*b/R_H
    mobility_error = mobility * np.sqrt((b_err/b)**2 + (R_H_err/R_H)**2)

    carrier_density = (mobility-mobility_value)/(1.6e-19 * R_H*(mobility-mobility_value))

    return mobility,mobility_error,carrier_density
    





def final_values(resistivity_file,hall_file,THICKNESS,carrier=1):
    resistivity,resistivity_error = calculate_resistivity(resistivity_file,THICKNESS)
    Res_0 = np.poly1d(resistivity)(0)
    Res_0_err = quadratic_error(Res_0,*(resistivity_error))
    hall_coefficient,hall_coefficient_error = calculate_hall(hall_file,THICKNESS)

    single_mobility, single_mobility_error,single_carrier_density = single_carrier(hall_coefficient,hall_coefficient_error,Res_0,Res_0_err)
    if carrier != 1:
        two_mobility, two_mobility_error,two_carrier_density = two_carrier(hall_coefficient,hall_coefficient_error,Res_0,Res_0_err,mobility_value=750)
        two_mobility2,two_mobility2_error,two_carrier_density2 = two_carrier_2(hall_coefficient,hall_coefficient_error,Res_0,Res_0_err,1e8*resistivity[0],resistivity_error[0],mobility_value=750)

    print(f'''Hall Coefficient = {hall_coefficient} ± {hall_coefficient_error} cm³ / C
    Resistivity at 0 field = {Res_0} ± {Res_0_err} Ωcm
    This gives a single carrier mobility of {single_mobility} ± {single_mobility_error} cm² / Vs
    and a carrier density of {single_carrier_density} cm⁻³''')

    if carrier != 1:
        print(f'''This gives an electron mobility of {two_mobility} ± {two_mobility_error} cm² / Vs
        and a carrier density of {two_carrier_density}
    This assumes no magnetoresistance.
    If we take into account magnetoresistance, we get an electron mobility of {two_mobility2}± {two_mobility2_error} cm² / Vs
    and a carrier density of {two_carrier_density2}''')

        


# print('n-type GaAs:')
# final_values('n-typeGaAs_RinB','n-typeGaAs_vary_B_2',THICKNESS=3e-6)
# print()
# print('p-type GaAs')
# final_values('p-typeGaAs_RinB','p-typeGaAs_vary_B_2',THICKNESS=3e-6)
# print()
print('InSb')
final_values('InSb_R','InSb_vary_B_2',THICKNESS=1e-6,carrier=2)


