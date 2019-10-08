from matplotlib import pyplot as plt
import numpy as np
import math as m

def pEpHplot(temp_C=25, P_O2=0.21, P_H2=1, logConc_Fe=-6):
    # theoretical limits of stability
    E_O2_standard = 1.229 # V
# gas constant (in useful units)
    R = 8.314/1000 # kJ/mol K
    # temperature (25°C in K)
    T = temp_C + 273.15 # K
    # Faraday constant (in useful units
    F = 96.5 # kJ mol-1
    # assuming pH = 7
    pH_range = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14]
    # partial pressure of O2
    PO2 = P_O2
    pE_O2_stability = [E_O2_standard/(2.303*R*T/(F)) -pH_range[i] + 1/4*m.log10(PO2) for i in range(0,len(pH_range))]

    # theoretical limits of stability
    E_H2O_standard = 0 # V
    PH2 = P_H2 # atm
    pE_H2O_stability = [E_H2O_standard/(2.303*R*T/(F)) -pH_range[i] + 1/2*m.log10(PH2) for i in range(0,len(pH_range))]

    # the standard cell potential difference (this is an experimental value one has to look up)
    E_Fe3_Fe2_standard = 0.78 # V
    # ratio of Fe3+ to Fe2+
    Fe3_Fe2 = 1 
    pE_Fe3_Fe2_array = [E_Fe3_Fe2_standard/(2.303*R*T/(F)) + 1/8*m.log10(Fe3_Fe2) for i in range(0,len(pH_range))]

    CFe = 10**(logConc_Fe)

    E_Fe3_Fe2_standard = 0.78 # V
    Ksp_FeOH2 = 4.8*10**(-17)
    Ksp_FeOH3 = 1.0*10**(-38)

    E_FeOH3_Fe3 = 0.289
    pH_FeOH3_Fe3 = 1/3*(E_FeOH3_Fe3/(2.303*R*T/(F)) - m.log10(CFe))

    E_FeOH3_Fe2 = 1.06
    pE_FeOH3_Fe2 = [E_FeOH3_Fe2/(2.303*R*T/(F)) - 3*pH_range[i] - m.log10(CFe) for i in range(0,len(pH_range))]

    E_FeOH2_Fe2 = 0.733
    pH_FeOH2_Fe2 = 1/2*(E_FeOH2_Fe2/(2.303*R*T/(F)) - m.log10(CFe))

    E_FeOH3_FeOH2 = 0.327
    pE_FeOH3_FeOH2 = [E_FeOH3_FeOH2/(2.303*R*T/(F)) - pH_range[i] for i in range(0,len(pH_range))]
    
    fig5 = plt.figure(figsize=(6, 8))
    pEpH = fig5.add_subplot(111)

    pH_FeOH3_Fe3_pH = [pH_FeOH3_Fe3, pH_FeOH3_Fe3]
    pH_FeOH3_Fe3_pE = [E_FeOH3_Fe2/(2.303*R*T/(F)) - 3*pH_FeOH3_Fe3 - m.log10(CFe), E_O2_standard/(2.303*R*T/(F)) -pH_FeOH3_Fe3 + 1/4*m.log10(PO2)]

    pH_FeOH2_Fe2_pH = [pH_FeOH2_Fe2, pH_FeOH2_Fe2]
    pH_FeOH2_Fe2_pE = [E_H2O_standard/(2.303*R*T/(F)) -pH_FeOH2_Fe2 + 1/2*m.log10(PH2), E_FeOH3_FeOH2/(2.303*R*T/(F)) - pH_FeOH2_Fe2]


    pE_Fe3_Fe2_pH = [0,pH_FeOH3_Fe3]
    pE_Fe3_Fe2_new = [pE_Fe3_Fe2_array[0],E_Fe3_Fe2_standard/(2.303*R*T/(F)) + 1/8*m.log10(Fe3_Fe2)]

    pE_FeOH3_Fe2_pH = [pH_FeOH3_Fe3, pH_FeOH2_Fe2]
    pE_FeOH3_Fe2_new = [E_FeOH3_Fe2/(2.303*R*T/(F)) - 3*pH_FeOH3_Fe3 - m.log10(CFe), E_FeOH3_Fe2/(2.303*R*T/(F)) - 3*pH_FeOH2_Fe2 - m.log10(CFe)]

    pE_FeOH3_FeOH2_pH = [pH_FeOH2_Fe2, 14]
    pE_FeOH3_FeOH2_new = [E_FeOH3_FeOH2/(2.303*R*T/(F)) - pH_FeOH2_Fe2, E_FeOH3_FeOH2/(2.303*R*T/(F)) - 14]

    pEpH.plot(pH_range,pE_O2_stability,  'g--', label='Upper limit of stability of water')
    pEpH.text(7, 18, 'Upper limit of stability of water')
    pEpH.plot(pH_range, pE_H2O_stability, 'g--', label='Lower limit of stability of water')
    pEpH.text(1, -12, 'Lower limit of stability of water')
    pEpH.plot(pE_Fe3_Fe2_pH, pE_Fe3_Fe2_new, 'g-', label='$Fe^{3+}$:$Fe^{2+}$')
    pEpH.text(1, 5, '$Fe^{2+}$')
    pEpH.plot(pH_FeOH3_Fe3_pH, pH_FeOH3_Fe3_pE, 'g-', label='$Fe^{3+}$:$FeOH_3$')
    pEpH.text(1, 17, '$Fe^{3+}$')
    pEpH.plot(pE_FeOH3_Fe2_pH, pE_FeOH3_Fe2_new, 'g-', label='$Fe^{3+}$:$Fe^{2+}$')
    pEpH.text(12, 0, '$Fe(OH)_3$')
    pEpH.plot(pH_FeOH2_Fe2_pH, pH_FeOH2_Fe2_pE, 'g-', label='$Fe^{2+}$:$FeOH_2$')
    pEpH.plot(pE_FeOH3_FeOH2_pH, pE_FeOH3_FeOH2_new, 'g-', label='$Fe^{3+}$:$Fe^{2+}$')
    pEpH.text(12, -10, '$Fe(OH)_2$')

    pEpH.set_ylabel('pE')
    pEpH.set_xlabel('pH')
    pEpH.set_ylim((-14, 20))
    pEpH.set_xlim((0, 14))
    
    plt.show()
