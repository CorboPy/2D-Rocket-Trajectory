# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 21:42:24 2023

@author: alexc

notes: to ensure same accuracy, the greater the length of time, the greater N must be (e.g. timex2 requires Nx2 points)
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.ticker as ticker

#f1,f2,f3,f4 planet only (no moon):
def dvxdt(t,x,y,vx,vy,M):
    """Computes dvx/dt (acceleration in x direction)
    

    Parameters
    ----------
    t : Float
        Time elapsed (s)
    x : Float
        Displacement in x axis (m)
    y : Float
        Displacement in y axis (m)
    vx : Float
        Velocity in direction of x axis (m/s)
    vy : Float
        Velocity in direction of y axis (m/s)
    M : Float
        Mass of planet (kg)

    Returns
    -------
    Acceleration in m/s^2 in x direction (dvx/dt) at time t given the rocket's displacement from origin and instantaneous velocity.

    """
    return(((-6.6743*10**(-11))*M*x)/((x**2+y**2)**(3/2)))
def dvydt(t,x,y,vx,vy,M):
    """Computes dvy/dt (acceleration in y direction)
    

    Parameters
    ----------
    t : Float
        Time elapsed (s)
    x : Float
        Displacement in x axis (m)
    y : Float
        Displacement in y axis (m)
    vx : Float
        Velocity in direction of x axis (m/s)
    vy : Float
        Velocity in direction of y axis (m/s)
    M : Float
        Mass of planet (kg)

    Returns
    -------
    Acceleration in m/s^2 in y direction (dvy/dt) at time t given the rocket's displacement from origin and instantaneous velocity.

    """
    return(((-6.6743*10**(-11))*M*y)/((x**2+y**2)**(3/2)))
def dxdt(vx):
    """Identity function used in the calculation of kx's
    

    Parameters
    ----------
    vx : Float
        Velocity in direction of x axis (m/s)

    Returns
    -------
    Velocity in direction of x axis (m/s)

    """
    return(vx)
def dydt(vy):
    """Identity function used in the calculation of ky's
    

    Parameters
    ----------
    vy : Float
        Velocity in direction of y axis (m/s)

    Returns
    -------
    Velocity in direction of y axis (m/s)

    """
    return(vy)

#f1,f2,f3,f4 earth+moon:
def dvxdt_moon(t,x,y,vx,vy,Mp,Mm):
    """Computes dvx/dt (acceleration in x direction) under gravitational force of two bodies (earth + moon)
    

    Parameters
    ----------
    t : Float
        Time elapsed (s)
    x : Float
        Displacement in x axis (m)
    y : Float
        Displacement in y axis (m)
    vx : Float
        Velocity in direction of x axis (m/s)
    vy : Float
        Velocity in direction of y axis (m/s)
    Mp : Float
        Mass of planet (kg)
    Mm : Float
        Mass of moon (kg)

    Returns
    -------
    Acceleration in m/s^2 in x direction (dvx/dt) at time t given the rocket's displacement from origin and instantaneous velocity.

    """
    ax_E=((-6.6743*10**(-11))*Mp*x)/((x**2+y**2)**(3/2))
    ax_M=((-6.6743*10**(-11))*Mm*x)/((x**2+(y-earthmoon)**2)**(3/2))
    a=ax_E + ax_M
    return(a)
def dvydt_moon(t,x,y,vx,vy,Mp,Mm):
    """Computes dvy/dt (acceleration in y direction) under gravitational force of two bodies (earth + moon)
    

    Parameters
    ----------
    t : Float
        Time elapsed (s)
    x : Float
        Displacement in x axis (m)
    y : Float
        Displacement in y axis (m)
    vx : Float
        Velocity in direction of x axis (m/s)
    vy : Float
        Velocity in direction of y axis (m/s)
    Mp : Float
        Mass of planet (kg)
    Mm : Float
        Mass of moon (kg)

    Returns
    -------
    Acceleration in m/s^2 in y direction (dvy/dt) at time t given the rocket's displacement from origin and instantaneous velocity.

    """
    ax_E=((-6.6743*10**(-11))*Mp*(y))/((x**2+y**2)**(3/2))
    ax_M=((-6.6743*10**(-11))*Mm*(y-earthmoon))/((x**2+(y-earthmoon)**2)**(3/2))
    a=ax_E + ax_M
    return(a)
def dxdt_moon(vx):
    """Identity function used in the calculation of kx's
    

    Parameters
    ----------
    vx : Float
        Velocity in direction of x axis (m/s)

    Returns
    -------
    Velocity in direction of x axis (m/s)

    """
    return(vx)
def dydt_moon(vy):
    """Identity function used in the calculation of ky's
    

    Parameters
    ----------
    vy : Float
        Velocity in direction of y axis (m/s)

    Returns
    -------
    Velocity in direction of y axis (m/s)

    """
    return(vy)

#time-stepping functions:
def xi(xprev,h,k1x,k2x,k3x,k4x):
    """Time-stepping function for displacement in x
    

    Parameters
    ----------
    xprev : Float
        Previous x coordinate
    h : Float
        Time-step
    k1x : Float
        First slope k1 for x
    k2x : Float
        Second slope k2 for x
    k3x : Float
        Third slope k3 for x
    k4x : Float
        Fourth slope k4 for x

    Returns
    -------
    Subsequent x coordinate.

    """
    newx = xprev+((h/6)*(k1x+(2*k2x)+(2*k3x)+k4x))
    return(newx)
def yi(yprev,h,k1y,k2y,k3y,k4y):
    """Time-stepping function for displacement in y
    

    Parameters
    ----------
    yprev : Float
        Previous y coordinate
    h : Float
        Time-step
    k1y : Float
        First slope k1 for y
    k2y : Float
        Second slope k2 for y
    k3y : Float
        Third slope k3 for y
    k4y : Float
        Fourth slope k4 for y

    Returns
    -------
    Subsequent y coordinate.

    """
    newy = yprev+((h/6)*(k1y+(2*k2y)+(2*k3y)+k4y))
    return(newy)
def vxi(xvprev,h,k1vx,k2vx,k3vx,k4vx):
    """Time-stepping function for velocity in x direction
    

    Parameters
    ----------
    xvprev : Float
        Previous velocity in x direction
    h : Float
        Time-step
    k1vx : Float
        First slope k1 for vx
    k2vx : Float
        Second slope k2 for vx
    k3vx : Float
        Third slope k3 for vx
    k4vx : Float
        Fourth slope k4 for vx

    Returns
    -------
    Subsequent velocity in x direction.

    """
    newxv = xvprev+((h/6)*(k1vx+(2*k2vx)+(2*k3vx)+k4vx))
    return(newxv)
def vyi(yvprev,h,k1vy,k2vy,k3vy,k4vy):
    """Time-stepping function for velocity in y direction
    

    Parameters
    ----------
    yvprev : Float
        Previous velocity in x direction
    h : Float
        Time-step
    k1vy : Float
        First slope k1 for vy
    k2vy : Float
        Second slope k2 for vy
    k3vy : Float
        Third slope k3 for vy
    k4vy : Float
        Fourth slope k4 for vy

    Returns
    -------
    Subsequent velocity in y direction.

    """
    newyv = yvprev+((h/6)*(k1vy+(2*k2vy)+(2*k3vy)+k4vy))
    return(newyv)

#pot. enregy funcitons (for part a and part b respectively)
def pot_a(x,y,mr,Mp):
    """Used for the calculation of gravitational potential energy in part a (planet only)
    

    Parameters
    ----------
    x : Float
        Displacement in x axis (m)
    y : Float
        Displacement in y axis (m)
    mr : Float
        Mass of rocket (kg)
    Mp : Float
        Mass of planet (kg)

    Returns
    -------
    Gravitational potential energy of rocket around planet (in joules) at coordinate (x,y)

    """
    pot = ((-6.6743*10**(-11))*Mp*(mr))/((x**2+y**2)**(1/2))
    return(pot)
def pot_b(x,y,mr,Mp,Mm):
    """Used for the calculation of gravitational potential energy in part b (planet + moon)
    

    Parameters
    ----------
    x : Float
        Displacement in x axis (m)
    y : Float
        Displacement in y axis (m)
    mr : Float
        Mass of rocket (kg)
    Mp : Float
        Mass of planet (kg)
    Mm : Float
        Mass of moon (kg)

    Returns
    -------
    Gravitational potential energy of rocket around planet and moon (in joules) at coordinate (x,y)

    """
    pot_E=((-6.6743*10**(-11))*Mp*(mr))/((x**2+y**2)**(1/2))
    pot_M=((-6.6743*10**(-11))*Mm*(mr))/((x**2+(y-earthmoon)**2)**(1/2))
    return(pot_E+pot_M)

#kin. energy funciton
def kin(vx,vy,mr):
    """Used for the calculation of kinetic energy of rocket
    

    Parameters
    ----------
    vx : Float
        Velocity in direction of x axis (m/s)
    vy : Float
        Velocity in direction of y axis (m/s)
    mr : Float
        Mass of rocket (kg)

    Returns
    -------
    Kinetic energy of rocket (in joules) at velocity vector (vx,vy)
    
    """

    kin_x=0.5*mr*vx**2
    kin_y=0.5*mr*vy**2
    return(kin_x+kin_y)

MyInput = '0'
while MyInput != 'q':
    MyInput = input('\nPlease input one of the following: \n \n - "a" for orbits around the Earth, Jupiter and Mars. \n - "b" for Earth - Moon trajectories. \n - "q" to quit. \n')
    if MyInput == 'a' or MyInput == 'A':
        print('You have chosen part (a): Orbital trajectories.\n')
        
        #planet parameters (simulation will start at the same height from surface of planet, regardless of planet)
        
        EarthRadius=6357.5*10**3
        EarthMass=5.972 * 10**24
        lowearthorbit = 7000*10**3
        
        MarsRadius= 3.3895 *10**6
        MarsMass= 6.39 * 10**23 
        MarsOrbit=MarsRadius++642.5*10**3
        
        JupiterRadius=69.911*10**6
        JupiterMass=1.899*10**27
        JupiterOrbit=JupiterRadius+642.5*10**3   
        
        while True:
            AInput = input('If you would like to edit the parameters, please type "yes". Otherwise, to use default parameters, please type "no": ')
            if AInput == "yes" or AInput == "no":
                if AInput == 'yes':
                    while True:     #testing input_T for a positive integer
                        Input_T = input('\nPlease enter a positive integer value for the time of flight (in hours): ')
                        try:
                            T = int(Input_T)
                        except ValueError:
                            print( 'Invalid number, please try again.')
                            continue
                        if 0 < T:
                            break
                        else:
                            print( 'Invalid number, please try again.')
                            
                    while True:     #testing input_N for a positive integer
                        Input_N = input('\nPlease enter a positive integer value (e.g. 10000) for N (the number of points to consider): ')
                        try:
                            numpoints = int(Input_N)
                        except ValueError:
                            print( 'Invalid number, please try again.')
                            continue
                        if 2 <= numpoints:
                            break
                        else:
                            print( 'Invalid number, please try again.')
                    noofsecs = T*3600
                
                    while True:     #testing velocity input for a float
                        velocity_input = input("\nPlease enter a number for the initial tangential velocity of the rocket (in meters per second): ")
                        try:
                            initialvelocity = float(velocity_input)
                            break
                        except ValueError:
                            print("Invalid input. Please enter a valid number.")
                
                    while True:     #user chooses planet
                        PInput = input('\nPlease choose your planet: \n\nType "e" for Earth, "j" for Jupiter, or "m" for Mars: ')
                        if PInput == "e" or PInput == "j" or PInput == 'm':
                            
                            if PInput == 'e':
                                planetname='Earth'
                                P_Radius = EarthRadius
                                m = EarthMass
                                rgb=(0.4, 0.7, 1.0)
                                start=lowearthorbit
                            
                            elif PInput == 'j':
                                planetname='JUPITER'
                                P_Radius = JupiterRadius
                                m = JupiterMass
                                rgb=rgb=(0.9,0.347,0)
                                start=JupiterOrbit
                                
                            else:
                                planetname='MARS'
                                P_Radius = MarsRadius
                                m = MarsMass
                                rgb=(1,0.1,0)
                                start=MarsOrbit
                                
                            break 
                        else:
                            print("\nInvalid input. Please try again.\n")
                
                elif AInput == 'no':    #default parameters
                    numpoints = 100000
                    noofsecs = 86400
                    
                    while True:     #user chooses planet:
                        PInput = input('\nPlease choose your planet: \n\nType "e" for Earth, "j" for Jupiter, or "m" for Mars: ')
                        if PInput == "e" or PInput == "j" or PInput == 'm':
                            
                            if PInput == 'e':
                                planetname='Earth'
                                P_Radius = EarthRadius
                                m = EarthMass
                                rgb=(0.4, 0.7, 1.0)
                                start=lowearthorbit
                                initialvelocity=9100
                            
                            elif PInput == 'j':
                                planetname='Jupiter'
                                P_Radius = JupiterRadius
                                m = JupiterMass
                                rgb=rgb=(0.9,0.347,0)
                                start=JupiterOrbit
                                initialvelocity=56000
                                
                            else:
                                planetname='Mars'
                                P_Radius = MarsRadius
                                m = MarsMass
                                rgb=(1,0.1,0)
                                start=MarsOrbit
                                initialvelocity=3400
                                
                            break 
                        else:
                            print("\nInvalid input. Please try again.\n")
                    
                break 
            else:
                print("\nInvalid input. Please try again.\n")
        
        print("\nCalculating, please wait...\n")
        
        h = tstep = noofsecs/numpoints #step size
        
        #initialising arrays
        xarray = np.zeros(numpoints)
        yarray = np.zeros(numpoints)
        vxarray = np.zeros(numpoints)
        vyarray = np.zeros(numpoints)
        tarray = np.linspace(0,noofsecs,numpoints)
        
        #initial conditions
        xarray[0] = 0
        yarray[0] = -start
        vxarray[0] = initialvelocity
        vyarray[0] = 0
        
        #progress variables
        I1 = round((numpoints-1)/4)
        I2 = round((numpoints-1)/2)
        I3 = round((3*(numpoints-1))/4)
        I4 = round((numpoints-1))
        
        #radius from planet
        r = np.sqrt((xarray[0]**2)+(yarray[0])**2)
        
        i=0
        index=0
        while r > P_Radius and i < numpoints-1: #conditions for loop to end - r < P_Radius (AKA crashing into planet)
            i = i+1
            
            k1 = [ dxdt(vxarray[i-1]), #k1x
                  dydt(vyarray[i-1]), #k1y
                  dvxdt(tarray[i-1],xarray[i-1],yarray[i-1],vxarray[i-1],vyarray[i-1],m),   #k1vx
                  dvydt(tarray[i-1],xarray[i-1],yarray[i-1],vxarray[i-1],vyarray[i-1],m) ]   #k1vy

            
            k2 = [ dxdt((vxarray[i-1])+(h*k1[2])/2), #k2x
                  dydt((vyarray[i-1])+(h*k1[3])/2), #k2y
                  dvxdt(tarray[i-1],(xarray[i-1])+((h*k1[0])/2),(yarray[i-1])+((h*k1[0])/2),vxarray[i-1],vyarray[i-1],m), #k2vx
                  dvydt(tarray[i-1],(xarray[i-1])+((h*k1[1])/2),(yarray[i-1])+((h*k1[1])/2),vxarray[i-1],vyarray[i-1],m) ] #k2vy

            
            k3 = [ dxdt((vxarray[i-1])+(h*k2[2])/2), #k3x
                  dydt((vyarray[i-1])+(h*k2[3])/2), #k3y
                  dvxdt(tarray[i-1],(xarray[i-1])+(h*k2[0])/2,(yarray[i-1])+(h*k2[0])/2,vxarray[i-1],vyarray[i-1],m), #k3vx
                  dvydt(tarray[i-1],(xarray[i-1])+(h*k2[1])/2,(yarray[i-1])+(h*k2[1])/2,vxarray[i-1],vyarray[i-1],m) ] #k3vy

            
            k4 = [ dxdt((vxarray[i-1])+(h*k3[2])), #k4x
                  dydt((vyarray[i-1])+(h*k3[3])), #k4y
                  dvxdt(tarray[i-1],(xarray[i-1])+(h*k3[0]),(yarray[i-1])+(h*k3[0]),vxarray[i-1],vyarray[i-1],m), #k4vx
                  dvydt(tarray[i-1],(xarray[i-1])+(h*k3[1]),(yarray[i-1])+(h*k3[1]),vxarray[i-1],vyarray[i-1],m) ] #k4vy

            xarray[i] = xi(xarray[i-1],h,k1[0],k2[0],k3[0],k4[0]) #x displacement
            yarray[i] = yi(yarray[i-1],h,k1[1],k2[1],k3[1],k4[1]) #y displacement
            vxarray[i] = vxi(vxarray[i-1],h,k1[2],k2[2],k3[2],k4[2]) #velocity in direction of x
            vyarray[i] = vyi(vyarray[i-1],h,k1[3],k2[3],k3[3],k4[3]) #velocity in direction of y
            
            #updating distance from planet (condition to break is in while r > P_Radius statement):
            r = np.sqrt((xarray[i-1]**2)+(yarray[i-1])**2)
            
            #progress print statements
            if i==I1:
                print("25% complete...")
            elif i==I2:
                print("50% complete...")
            elif i==I3:
                print("75% complete...")
            elif i==I4:
                print("100% complete...")
            else:
                continue
        
        #setting up the arrow for plot later (depends on where i gets to in the loop in the case of rocket crashing)
        arrow = patches.FancyArrowPatch((xarray[i-2], yarray[i-2]), (xarray[i-1], yarray[i-1]), arrowstyle='simple', mutation_scale=10, color='r')
        
        #energy conservation
        potinitial=pot_a(xarray[0],yarray[0],549054,m)
        potfinal=pot_a(xarray[i],yarray[i],549054,m)  
        kininitial=kin(vxarray[0],vyarray[0],549054)
        kinfinal=kin(vxarray[i],vyarray[i],549054)
        initialenergy=potinitial+kininitial
        finalenergy=potfinal+kinfinal
        energydiff=finalenergy-initialenergy
        print("\nPercentage energy difference (percentage uncertainty) =",abs(100*((finalenergy-initialenergy)/initialenergy)),"%")

        
        #'rocket has crashed' info message
        timeofcol=(tarray[i]/3600)
        if xarray[numpoints-1]==0 or yarray[numpoints-1]==0:    #since rocket will never reach (0,0) unless the loop has prematurely been broken as a result of rocket crashing - in this case there will still be zeros left in the array
            print("\n!! THE ROCKET HAS CRASHED INTO",planetname,"!!")
            print("Time of collision = %s" % timeofcol,"hours into flight.")
        
        #if rocket has crashed, this will assign the remaining zeros in the y and x arrays to the final position vector of the rocket before crashing
        limitx = xarray[i]
        limity = yarray[i]
        while i != numpoints-1:
            i=i+1
            xarray[i]= limitx
            yarray[i]= limity    

        #finding max distance rocket travels in order to set the values of x and y ticks
        max_x = np.abs(xarray.max())
        max_y = np.abs(yarray.max())
        max_distance = np.sqrt((max_x)**2+(max_y)**2)
        if max_distance < (P_Radius+2*10**6):
            tickslimit = P_Radius+3*10**6
        else:
            xmax=max_distance
            xmin=max_distance
            ymax=max_distance
            ymin=max_distance
            tickslimit = max_distance

        
        
        fig, ax = plt.subplots(dpi=200)
        #setting x and y ticks limits
        ax.set_xlim([-tickslimit, tickslimit])
        ax.set_ylim([-tickslimit, tickslimit])
        
        #adding planet at origin:
        circle = patches.Circle((0, 0), radius=(P_Radius), facecolor=rgb, linewidth=1, edgecolor='black')
        ax.add_patch(circle)
        sizefont=15*(start+0.3*(start))/tickslimit
        ax.text(0,0,planetname,horizontalalignment='center',color='w',fontsize=sizefont,fontname="Lucida Sans Unicode")
        
        #labels:
        ax.set_xlabel('x (m)',fontname="Calibri")
        ax.set_ylabel('y (m)',fontname="Calibri")
        plt.title("Rocket Trajectory From Orbit Radius of %s km" % int(start/1000),weight='bold',fontname="Calibri")
        ax.xaxis.set_tick_params(labelsize=7)
        ax.yaxis.set_tick_params(labelsize=7)
        ax.xaxis.set_major_formatter(ticker.StrMethodFormatter('{x:.0e}'))
        ax.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:.0e}'))
        txt="(for initial tangential velocity = %s m/s, time of flight = %s hours, N = %s)" % (initialvelocity, round(timeofcol,3), numpoints)
        plt.figtext(0.5, -0.05, txt, wrap=True, horizontalalignment='center', fontsize=7, style='italic')
        
        #plotting:
        ax.plot(xarray,yarray,color='red',linewidth=1,linestyle='dotted')
        ax.set_aspect('equal')

        ax.add_patch(arrow)        
        
        plt.show()
    
    
        
    elif MyInput == 'b' or MyInput == 'B':
        print('You have chosen part (b): Earth - Moon trajectories.\n')
        
        EarthRadius=6357.5*10**3
        EarthMass=5.972 * 10**24
        lowearthorbit = 7000*10**3
        
        earthmoon=384400*10**3
        M_radius= 1737.4*10**3
        M_mass = 7.34767309 * 10**22
        
        while True:
            AInput = input('If you would like to edit the parameters, please type "yes". Otherwise, to use predetermined parameters for a slingshot around the moon, please type "no": ')
            if AInput == "yes" or AInput == "no":
                if AInput == 'yes':
                    while True:     #testing input_T for a positive integer
                        Input_T = input('\nPlease enter a positive integer value for the time of flight in hours (e.g. ~ 288 for Earth-Moon-Earth): ')
                        try:
                            T = int(Input_T)
                        except ValueError:
                            print( 'Invalid number, please try again.')
                            continue
                        if 0 < T:
                            break
                        else:
                            print( 'Invalid number, please try again.')
                            
                    while True:     #testing input_N for a positive integer (and above 2 as simulation is simply useless if N=1)
                        Input_N = input('\nPlease enter a positive integer value (e.g. 1000000) for N (the number of points to consider): ')
                        try:
                            numpoints = int(Input_N)
                        except ValueError:
                            print( 'Invalid number, please try again.')
                            continue
                        if 2 <= numpoints:
                            break
                        else:
                            print( 'Invalid number, please try again.')
                    noofsecs = T*3600
                
                    while True:     #testing velocity input for a float
                        velocity_input = input("\nPlease enter a number for the initial tangential velocity of the rocket in meters per second (sweet spot is ~ 10562.35 m/s): ")
                        try:
                            initialvelocity = float(velocity_input)
                            break
                        except ValueError:
                            print("Invalid input. Please enter a valid number.")
                    
                    moonname='THE MOON'
                    planetname='EARTH'
                    P_Radius = EarthRadius
                    m = EarthMass
                    rgb=(0.4, 0.7, 1.0)
                    start=lowearthorbit
                    
                elif AInput == 'no':    #default parameters
                    numpoints = 1000000
                    noofsecs = 4*259200
                    moonname='THE MOON'
                    planetname='EARTH'
                    P_Radius = EarthRadius
                    m = EarthMass
                    rgb=(0.4, 0.7, 1.0)
                    start=lowearthorbit 
                    initialvelocity=10562.35
                    
                    
                break 
            else:
                print("\nInvalid input. Please try again.\n")
        
        print("\nCalculating, please wait...\n")
        
        mm=M_mass
    
        #P_Radius=6357.5*10**3
        #m = 5.972 * 10**24
        
        h = tstep = noofsecs/numpoints #step size
        
        #initialising arrays
        xarray = np.zeros(numpoints)
        yarray = np.zeros(numpoints)
        vxarray = np.zeros(numpoints)
        vyarray = np.zeros(numpoints)
        tarray = np.linspace(0,noofsecs,numpoints)
        
        #initial conditions
        xarray[0] = 0
        yarray[0] = -start
        vxarray[0] = initialvelocity
        vyarray[0] = 0
        
        #progress variables
        I1 = round((numpoints-1)/4)
        I2 = round((numpoints-1)/2)
        I3 = round((3*(numpoints-1))/4)
        I4 = round((numpoints-1))
        
        #radius from planet
        r = np.sqrt((xarray[0]**2)+(yarray[0])**2)
        
        body='no crash'
        i=0
        while i < numpoints-1: #condition for loop to complete (fill arrays)
            
            i = i+1
            
            k1 = [ dxdt_moon(vxarray[i-1]), #k1x
                  dydt_moon(vyarray[i-1]), #k1y
                  dvxdt_moon(tarray[i-1],xarray[i-1],yarray[i-1],vxarray[i-1],vyarray[i-1],m,mm),   #k1vx
                  dvydt_moon(tarray[i-1],xarray[i-1],yarray[i-1],vxarray[i-1],vyarray[i-1],m,mm) ]   #k1vy

            
            k2 = [ dxdt_moon((vxarray[i-1])+(h*k1[2])/2), #k2x
                  dydt_moon((vyarray[i-1])+(h*k1[3])/2), #k2y
                  dvxdt_moon(tarray[i-1],(xarray[i-1])+((h*k1[0])/2),(yarray[i-1])+((h*k1[0])/2),vxarray[i-1],vyarray[i-1],m,mm), #k2vx
                  dvydt_moon(tarray[i-1],(xarray[i-1])+((h*k1[1])/2),(yarray[i-1])+((h*k1[1])/2),vxarray[i-1],vyarray[i-1],m,mm) ] #k2vy

            
            k3 = [ dxdt_moon((vxarray[i-1])+(h*k2[2])/2), #k3x
                  dydt_moon((vyarray[i-1])+(h*k2[3])/2), #k3y
                  dvxdt_moon(tarray[i-1],(xarray[i-1])+(h*k2[0])/2,(yarray[i-1])+(h*k2[0])/2,vxarray[i-1],vyarray[i-1],m,mm), #k3vx
                  dvydt_moon(tarray[i-1],(xarray[i-1])+(h*k2[1])/2,(yarray[i-1])+(h*k2[1])/2,vxarray[i-1],vyarray[i-1],m,mm) ] #k3vy

            
            k4 = [ dxdt_moon((vxarray[i-1])+(h*k3[2])), #k4x
                  dydt_moon((vyarray[i-1])+(h*k3[3])), #k4y
                  dvxdt_moon(tarray[i-1],(xarray[i-1])+(h*k3[0]),(yarray[i-1])+(h*k3[0]),vxarray[i-1],vyarray[i-1],m,mm), #k4vx
                  dvydt_moon(tarray[i-1],(xarray[i-1])+(h*k3[1]),(yarray[i-1])+(h*k3[1]),vxarray[i-1],vyarray[i-1],m,mm) ] #k4vy

            xarray[i] = xi(xarray[i-1],h,k1[0],k2[0],k3[0],k4[0]) #x displacement
            yarray[i] = yi(yarray[i-1],h,k1[1],k2[1],k3[1],k4[1]) #y displacement
            vxarray[i] = vxi(vxarray[i-1],h,k1[2],k2[2],k3[2],k4[2]) #v in x direction
            vyarray[i] = vyi(vyarray[i-1],h,k1[3],k2[3],k3[3],k4[3]) #v in y direction
            
            #updating distance of rocket from earth and moon respectively:
            r_E = np.sqrt((xarray[i-1]**2)+(yarray[i-1])**2)
            r_M = np.sqrt((xarray[i-1])**2+(yarray[i-1]-earthmoon)**2)
            
            #tests to see if rocket has collided with earth or moon
            if r_E < P_Radius:
                body=planetname
                break
            if r_M < M_radius:
                body=moonname
                break
            
            #progress print statements
            if i==I1:
                print("25% complete...")
            elif i==I2:
                print("50% complete...")
            elif i==I3:
                print("75% complete...")
            elif i==I4:
                print("100% complete...")
            else:
                continue
        
        #setting up the arrow for plot later (depends on where i gets to in the loop in the case of rocket crashing)
        if body!='no crash':
            arrow = patches.FancyArrowPatch((xarray[i-2], yarray[i-2]), (xarray[i-1], yarray[i-1]), arrowstyle='simple', mutation_scale=0, color='r')
        else:
            arrow = patches.FancyArrowPatch((xarray[i-2], yarray[i-2]), (xarray[i-1], yarray[i-1]), arrowstyle='simple', mutation_scale=10, color='r')

        #energy conservation
        potinitial=pot_b(xarray[0],yarray[0],549054,EarthMass,M_mass)
        potfinal=pot_b(xarray[i],yarray[i],549054,EarthMass,M_mass)  
        kininitial=kin(vxarray[0],vyarray[0],549054)
        kinfinal=kin(vxarray[i],vyarray[i],549054)
        initialenergy=potinitial+kininitial
        finalenergy=potfinal+kinfinal
        energydiff=finalenergy-initialenergy
        print("\nPercentage energy difference (percentage uncertainty) =",abs(100*((finalenergy-initialenergy)/initialenergy)),"%")
    
        #finding suitable x and y ticks limits
        x_arraymax = xarray.max()
        x_arraymin = xarray.min()
        y_arraymax = yarray.max()
        y_arraymin = yarray.min()
        
        if body=='EARTH':
            ytickmin = xtickmin = -start-0.3*(start)    #rocket has crashed before leaving LEO
            xtickmax = ytickmax = start+0.3*(start)
            for a in vyarray:
                if a < 0:   #rocket has gone around moon before crashing
                    ytickmax = y_arraymax+5*10**7
                    ytickmin = y_arraymin-5*10**7
                    xtickmin = -0.5*(ytickmax-ytickmin)
                    xtickmax = 0.5*(ytickmax-ytickmin)
                    break
        else:   #rocket has not crashed OR rocket has crashed into Moon
            ytickmax = y_arraymax+5*10**7
            ytickmin = y_arraymin-5*10**7
            xtickmin = -0.5*(ytickmax-ytickmin)
            xtickmax = 0.5*(ytickmax-ytickmin)

        #"rocket has crashed into planet/moon" print statements
        timeofcol=(tarray[i]/3600)
        if xarray[numpoints-1]==0 or yarray[numpoints-1]==0:
            print("\n!! THE ROCKET HAS CRASHED INTO",body,"!!")
            print("Time of collision = %s" % (tarray[i]/3600),"hours into flight.")
            
        
        #if rocket has crashed, this will assign the remaining zeros in the y and x arrays to the final position vector of the rocket before crashing
        limitx = xarray[i]
        limity = yarray[i]
        while i != numpoints-1:
            i=i+1
            xarray[i]= limitx
            yarray[i]= limity    
        
        # plt.figure(dpi=100)
        fig, ax = plt.subplots(dpi=200)
        #setting x and y ticks limits
        ax.set_xlim([xtickmin, xtickmax])
        ax.set_ylim([ytickmin, ytickmax])
        
        #adding planet at origin:
        circle = patches.Circle((0, 0), radius=(P_Radius), facecolor=rgb, linewidth=1, edgecolor='black')
        ax.add_patch(circle)
        sizefont=23*(start+0.3*(start))/xtickmax
        ax.text(0,0,"Earth",horizontalalignment='center',color='w',fontsize=sizefont,fontname="Lucida Sans Unicode")
        
        #adding moon on y axis:
        circle = patches.Circle((0, earthmoon), radius=(M_radius), facecolor=(0.3,0.3,0.3), linewidth=0.75, edgecolor='black')
        ax.add_patch(circle)
        
        #labels:
        ax.set_xlabel('x (m)',fontname="Calibri")
        ax.set_ylabel('y (m)',fontname="Calibri")
        plt.title("Rocket Trajectory From Orbit Radius of %s km" % int(start/1000),weight='bold',fontname="Calibri")
        ax.xaxis.set_tick_params(labelsize=7)
        ax.yaxis.set_tick_params(labelsize=7)
        ax.xaxis.set_major_formatter(ticker.StrMethodFormatter('{x:.0e}'))
        ax.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:.0e}'))
        txt="(for initial tangential velocity = %s m/s, time of flight = %s hours, N = %s)" % (initialvelocity, round(timeofcol,3), numpoints)
        plt.figtext(0.5, -0.05, txt, wrap=True, horizontalalignment='center',verticalalignment='center', fontsize=7, style='italic')
        
        #plotting:
        ax.plot(xarray,yarray,color='red',linewidth=1,linestyle='dotted')
        ax.set_aspect('equal')
        ax.add_patch(arrow)  
        
        
        plt.show()
        
        
        
    elif MyInput != 'q':
        print('This is not a valid choice.\n')

print('\nYou have chosen to finish - goodbye!')