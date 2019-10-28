from matplotlib import pyplot as plt
import numpy as np
import math as m

def CarbonDiagram1(temp_C=25):
    #defining the pH range
    pH_range = [0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5,12,12.5,13,13.5,14]
    
    T = temp_C + 273.15
    
    pK1 =  3404.71/T + 0.032786*T - 14.8435 #(Harned and Davis, 1943) 
    pK2 =  2902.39/T  + 0.02379*T - 6.4980 #(Harned and Scholes, 1941) 

    D = [(10**(-pH_range[i]))**2 + (10**(-pK1))*(10**(-pH_range[i])) + (10**(-pK1))*10**(-pK2)  for i in range(0,len(pH_range))]
    
    H2CO3 = [(10**(-pH_range[i]))**2/D[i] for i in range(0,len(pH_range))]
    HCO3 = [(10**(-pK1))*(10**(-pH_range[i]))/D[i] for i in range(0,len(pH_range))]
    CO3 = [(10**(-pK1))*(10**(-pK2))/D[i] for i in range(0,len(pH_range))]
    
    fig = plt.figure(figsize=(6, 8))
    carbon = fig.add_subplot(111)
    
    carbon.plot(pH_range,H2CO3,  'g-', label='$H_2CO_3$')
    carbon.plot(pH_range, HCO3, 'r-', label='$HCO_3^-$')
    carbon.plot(pH_range, CO3, 'b-', label='$CO_3^{2-}$')
    
    carbon.set_ylabel('Fraction of carbon as specified species')
    carbon.set_xlabel('pH')
    carbon.set_ylim((-.01, 1.01))
    carbon.set_xlim((0, 14))
    
    plt.legend()
    # Put a legend to the right of the current axis
    carbon.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.show()
    
def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

def Plot_A4(HOD, Temp_C, CO2, pH_f):
    fig, host = plt.subplots()
    fig.subplots_adjust(right=0.75)

    par1 = host.twinx()
    par2 = host.twinx()

    # Offset the right spine of par2.  The ticks and label have already been
    # placed on the right by twinx above.
    par2.spines["right"].set_position(("axes", 1.2))
    # Having been created by twinx, par2 has its frame off, so the line of its
    # detached spine is invisible.  First, activate the frame but make the patch
    # and spines invisible.
    make_patch_spines_invisible(par2)
    # Second, show the right spine.
    par2.spines["right"].set_visible(True)

    p1, = host.plot(HOD, Temp_C, "b-", label="Temperature (C)")
    p2, = par1.plot(HOD, CO2, "r-", label="$CO_2$ (ppm)")
    p3, = par2.plot(HOD, pH_f, "g-", label="pH")



    host.set_xlabel("Hour of Day")
    host.set_ylabel("Temperature (C)")
    par1.set_ylabel("$CO_2$ (ppm)")
    par2.set_ylabel("pH")

    host.yaxis.label.set_color(p1.get_color())
    par1.yaxis.label.set_color(p2.get_color())
    par2.yaxis.label.set_color(p3.get_color())

    tkw = dict(size=4, width=1.5)
    host.tick_params(axis='y', colors=p1.get_color(), **tkw)
    par1.tick_params(axis='y', colors=p2.get_color(), **tkw)
    par2.tick_params(axis='y', colors=p3.get_color(), **tkw)
    host.tick_params(axis='x', **tkw)

    lines = [p1, p2, p3]


    plt.show()
    
def Plot_A4_2(HOD, Temp_C, CO2):
    fig, host = plt.subplots()
    fig.subplots_adjust(right=0.75)

    par1 = host.twinx()

    p1, = host.plot(HOD, Temp_C, "b-", label="Temperature (C)")
    p2, = par1.plot(HOD, CO2, "r-", label="$CO_2$ (ppm)")
  

    host.set_xlabel("Hour of Day")
    host.set_ylabel("Temperature (C)")
    par1.set_ylabel("$CO_2$ (ppm)")

    host.yaxis.label.set_color(p1.get_color())
    par1.yaxis.label.set_color(p2.get_color())

    tkw = dict(size=4, width=1.5)
    host.tick_params(axis='y', colors=p1.get_color(), **tkw)
    par1.tick_params(axis='y', colors=p2.get_color(), **tkw)
    host.tick_params(axis='x', **tkw)

    lines = [p1, p2]


    plt.show()