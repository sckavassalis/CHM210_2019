from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import math as m
from scipy import integrate

def solve_h2o(relhum=70.0, tempK=273.15, pressure=1):
    """Convert relative humidity to gas phase water concentration"""
    #Converts relative humidity into H2O mixing ratio
    temp = tempK - 273.15

    #Method from Lowe, P.R. and J.M. Ficke, 1974: The computation of saturation vapor pressure. Tech.
    #Paper No. 4-74, Environmental Prediction Research Facility, Naval Postgraduate School,
    #Monterey, CA, 27 pp
    aow = 6.107799961 
    a1w = 0.4436518521
    a2w = .01428945805
    a3w = 2.650648471*10**(-4)
    a4w = 3.031240396*10**(-6)
    a5w = 2.034080948*10**(-8)
    a6w = 6.136820929*10**(-11)

    aoi = 6.109177956 
    a1i = 0.5034698970
    a2i = .01886013408
    a3i = 4.176223716*10**(-4)
    a4i = 5.824720280*10**(-6)
    a5i = 4.838803174*10**(-8)
    a6i = 1.838826904*10**(-10)

    ewater = aow + temp*(a1w+temp*(a2w+temp*(a3w+temp*(a4w+temp*(a5w+temp*a6w)))))
    eice = aoi + temp*(a1i+temp*(a2i+temp*(a3i+temp*(a4i+temp*(a5i+temp*a6i)))))
    eLowe = min(ewater, eice)

    return relhum*eLowe*(10**(-6))*(6.022*10**23)/(8.31447215 *(temp+273.1))

def estimate_PHOx(temp=273.15, pressure=1, relhum=70, o3=60, h2co=4, HOD=12):
    #temp in K, P in atm, [o3] in ppb, [h2co] in ppb, relhum in percent
    #O3 + hv -> O(1D) + O2 				(R1)
    #H2CO + hv ->(O2) CHO + HO2		(R2a)
    #CHO + O2 -> CO + HO2				(R2b)
    #O(1D) + H2O -> 2OH				(R3)
    #O(1D) + O2 -> O(3P) + O2			(R4a)
    #O(1D) +N2 -> O(3P) + N2			(R4b)

    if HOD <= 5: 
        J1 = 0
        J2 = 0
    if HOD > 5 and HOD <=6: 
        J1 = 2.162E-08
        J2 = 2.563E-07
    if HOD > 6 and HOD <=7:
        J1 = 9.258E-07
        J2 = 7.504E-06
    if HOD > 7 and HOD <=8:
        J1 = 5.613E-06
        J2 = 2.207E-05
    if HOD > 8 and HOD <=9:
        J1 = 1.492E-05
        J2 = 3.584E-05
    if HOD > 9 and HOD <=10:
        J1 = 2.545E-05
        J2 = 4.554E-05
    if HOD > 10 and HOD <=11: 
        J1 = 3.363E-05
        J2 = 5.122E-05
    if HOD > 11 and HOD <=12:
        J1 = 3.685E-05
        J2 = 5.316E-05
    if HOD > 12 and HOD <=13:
        J1 = 3.432E-05
        J2 = 5.164E-05   
    if HOD > 13 and HOD <=14:
        J1 = 2.801E-05
        J2 = 4.713E-05
    if HOD > 14 and HOD <=15:
        J1 = 1.713E-05
        J2 = 3.759E-05
    if HOD > 15 and HOD <=16: 
        J1 = 7.131E-06
        J2 = 2.442E-05
    if HOD > 16 and HOD <=17:
        J1 = 1.477E-06
        J2 = 9.910E-06
    if HOD > 17 and HOD <=18:
        J1 = 6.578E-08
        J2 = 7.882E-07
    if HOD > 18:
        J1 = 0
        J2 = 0
    k3 = 1.62E-10*m.exp(0.54/(.0083145*temp)) # 235 - 370 K, 2004DUN/RAV3333-3340
    k4a = 3.2E-11*m.exp(0.56/(.0083145*temp)) #  200 - 350 K, 1997ATK/BAU1329-1499
    k4b = 1.79E-11*m.exp(0.89/(.0083145*temp)) # 100 - 350 K, 1997ATK/BAU1329-1499

    o3mpercc = o3/((10**9/(6.0221413*10**23)) * 82.05736*temp/pressure)  # converts ppb to molecules/cm^3
    h2compercc = h2co/((10**9/(6.0221413*10**23)) * 82.05736*temp/pressure)  # converts ppb to molecules/cm^3

    na =  ((6.0221413*10**23)*pressure/(82.05736*temp)) # number density of air in molecules/cm^3
    h2ompercc = solve_h2o(relhum=relhum, tempK=temp, pressure=pressure)

    PHOxmperccpersec = 2*J1*o3mpercc*k3*h2ompercc/(k3*h2ompercc+.785*na*k4b+.215*na*k4a)+2*J2*h2compercc
    PHOxppbperhour = PHOxmperccpersec*((10**9/(6.0221413*10**23)) * 82.05736*temp/pressure)*3600
    return PHOxmperccpersec

def ozone_prod_contours( temp_K=273.15, RH_percent=70, ozone_ppb=60, H2CO_ppb=4, hour_of_day=12):
    pressure=1 # assumes pressure is 1 atm
    NO2oNO = 3 # ratio of NO2 to NO
    alpha=0.04 # from Cleary
    temp = temp_K
    relhum = RH_percent
    o3 = ozone_ppb
    h2co = H2CO_ppb
    HOD = hour_of_day 
        
    #NO2 + OH + M -> HNO3 + M			(R7), k7
    k7 = 2.6*10**(-30)*(temp/298)**(-2.9) # 1997ATK/BAU1329-1499
    #RO2 + NO -> RO + NO2					(R10), k10effective
    k10eff = 8*10**(-12) #Tyndall (2001), Sander (2003), Geddes (2009)
    #RO2 + R'O2 -> ROOR' + O2				(R11), k11effective
    k11eff = 5*10**(-12) #Tyndall (2001), Sander (2003), Geddes (2009)
    #HO2 + NO -> NO2 + OH					(R8), k8
    k8 = 3.6*10**(-12)*m.exp(2.24/(.0083145*temp)) # 2004ATK/BAU1461-1738
    #HO2 + HO2 + M -> H2O2 + O2 + M		(R9), k9
    k9 = 1.9*10**(-33)*m.exp(8.15/(.0083145*temp)) # 1997ATK/BAU1329-1499
  
    h2oconc = solve_h2o(relhum=relhum, tempK=temp, pressure=pressure)
    o3conc = o3/((10**9/(6.0221413*10**23)) * 82.05736*temp/pressure)
    
    Na =  ((6.0221413*10**23)*pressure/(82.05736*temp)) # number density of air in molecules/cm^3
  
    PHOx = estimate_PHOx(temp=temp, pressure=pressure, relhum=relhum, o3=o3, h2co=h2co, HOD=HOD)

    num_no2=100
    num_vocr=100
    NO2 = np.zeros(num_no2)
    NO2ppb = np.zeros(num_no2)
    logNO2ppb = np.zeros(num_no2)
    VOCR = np.zeros(num_vocr)
    A = np.zeros((num_no2,num_vocr))
    B = np.zeros((num_no2,num_vocr))
    C = np.zeros((num_no2,num_vocr))
    OH = np.zeros((num_no2,num_vocr))
    o3_prod_contours = np.zeros((num_no2,num_vocr))
    for l in range(num_no2):
        NO2ppb[l] = (l+0.1)*.6
        NO2[l] = NO2ppb[l]*Na/10**9
        logNO2ppb[l] = m.log10(NO2ppb[l])
        for n in range(num_vocr):
            VOCR[n] = n*.06 + 0.1
            A[n][l] = k7*NO2[l]*Na+alpha*VOCR[n]
            B[n][l] = (k7*NO2[l]*Na+alpha*VOCR[n])**2 + 24*PHOx*k11eff*(VOCR[n]/(k10eff*NO2[l]/NO2oNO))**2
            C[n][l] = 12*k11eff*(VOCR[n]/(k10eff*NO2[l]/NO2oNO))**2
            OH[n][l] = ((-A[n][l]+m.sqrt(B[n][l]))/C[n][l])
            o3_prod_contours[n][l] = 2*VOCR[n]*OH[n][l]*((10**9/(6.0221413*10**23))*82.05736*temp/pressure)*3600 
    #return NO2ppb, VOCR, o3_prod_contours
    
    cp = plt.contourf(logNO2ppb, VOCR, o3_prod_contours)
    cbar = plt.colorbar(cp)

    plt.title('Ozone Isopleth Graph')
    plt.xlabel('log($NO_2$ ppb/ppb) x ppb')
    plt.ylabel('VOC reactivity ($s^{-1}$)')
    plt.clabel(cp, colors = 'k', fmt = '%2.1f', fontsize=12)
    cbar.set_label('$O_3$ production (ppb/hr)', rotation=270, labelpad=15)
    plt.show()