import math
#1. gas constants and ion concentrations
R = 8314.0                                      #gas constant                          
T = 306.15                                      #temperature
F = 96487.0                                     #farady constant
RTF_con = R * T / F                             #define as RTF constant   
FRT_con = F/(R*T)

Na_o = outside_sodium_concentration = 140.0     #outside_sodium_concentration      
Na_i = inside_sodium_concentration  =18.0       #inside_sodium_concentration

K_o = outside_potassium_concentration = 5.4     #outside_potassium_concentration
K_i = inside_potassium_concentration = 145.0    #inside_potassium_concentration

Ca_o = outside_calcium_concentration = 1.8      #outside_calcium_concentration
Ca_i = inside_calcium_concentration = 0.0002    #inside_calcium_concentration
#########################____________________finish

#2. calculate Equilibrium
def E(outside,inside):
    '''
    This block is to calculate the Ek of ions.
    Only for mono-valent ions.
    It takes in 2 arguments.
        (1).The first  arg is the EXtracellular concentration of the ion.
        (2).The second arg is the INtracellular concentration of the ion
    For example, in order to calculate the equilibrium potential of sodium channel, try:
        E(137,8)
    '''
    import math
    R = 8314.0
    T = 306.15
    F = 96487.0
    rtf_con = R*T/F
    rate_con = 1.0 * outside/inside
    Ek = rtf_con*math.log(rate_con)
    return Ek
#########################____________________finish

#3.stimulating parameters.
'''
this block try to build the Iapp
I_duration = deltaT
I_amplitude= I[i]            #lower_i is not current, but counter
'''
deltaT=0.01                  #deltaT=0.01 S
t = []                       #set t
for T in range(1000000):     #t = sep(0,1000,0.01)
    t.append(T/100.0)
I = []
for i in range(1000000):     #set I as t
    I.append(0)
for i in range(500):         #blank the beginning interval
    I[i] = 0
for i in range(10500,10700): #begin to stimulate
    I[i] = 20

for i in range(75500,75700): 
    I[i] = 100
#########################____________________finish

#4. Equilibrium Value and G_ion(conductance)
#4a: E_Na_complex   
gbar_Na = 23
E_Na    = E(Na_o,Na_i)

#4b: E_si_complex
gbar_si = 0.09
E_si    = 7.7 - 13.0827 * math.log(Ca_i)

#4c: E_K_complex
gbar_K  = 0.282  * math.sqrt(K_o/5.4)
PR_NaK = 0.01833
E_K_numerator   = K_o + PR_NaK * Na_o
E_K_denominator = K_i + PR_NaK * Na_i
E_K = RTF_con * math.log(E_K_numerator/E_K_denominator)

#4d: E_K1_complex
gbar_K1 = 0.6047 * math.sqrt(K_o/5.4)
E_K1    = E(K_o,K_i)

#4e: E_Kp_complex
gbar_Kp = 0.0183 
E_Kp    = E(K_o,K_i)

#4f: E_b_complex
gbar_b  = 0.03921
E_b     = -59.87
#########################____________________finish

#5. parameters of activations and inactivations
#5a. I_Na
def alpha_m(V):
    import math
    numerator   = 0.32 * (V + 47.13)
    denominator =   1  - math.e**(-0.1 * (V + 47.13) )
    return numerator / denominator
def beta_m(V):
    import math
    return 0.08 * math.e**( -V / 11.0)
def alpha_h(V):
    import math
    if V >= -40:
        return 0
    elif V < -40:
        return 0.135 * math.e**((80 + V) / -6.8)
def beta_h(V):
    import math
    if V >= -40:
        numerator_1   = 1.0 
        denominator_1 = 0.13 * (1 + math.e**((V + 10.66)/-11.1))
        return 1 / denominator_1
    elif V < -40:
        component_1 = 3.56 * math.e**(0.079 * V) 
        component_2 = 3.1 * 100000 * math.e**(0.35 * V)
        return component_1 + component_2
def alpha_j(V):
    import math
    if V >= -40:
        return 0
    elif V < -40:
        component_1 = -1.2714 * 100000  * math.e**(0.2444 * V)
        component_2 = 3.474   * 0.00001 * math.e**( - 0.04391 * V)
        denominator =     1   + math.e **(0.311 * (V + 79.23))
        result = (component_1 - component_2) * (V + 37.78)/ denominator
        return result
def beta_j(V):
    import math
    if V >= -40:
        numerator   = 0.3 * math.e**(-2.535 * 0.0000001 * V)
        denominator =   1 + math.e**(-0.1 * (V + 32.0))
        return numerator/denominator
    elif V < -40:
        numerator   = 0.1212 * math.e**(-0.01052 *  V)
        denominator =  1     + math.e**(-0.1378  * (V + 40.14))
        return numerator / denominator
#########################____________________finish
    
#5b. I_si
def alpha_d(V):
    import math
    numerator   = 0.095 * math.e**(-0.01  * (V - 5 )) 
    denominator =    1  + math.e**(-0.072 * (V - 5 ))
    return numerator / denominator
def beta_d(V):
    import math
    numerator   =0.07   * math.e**(-0.017 * (V + 44)) 
    denominator =  1    + math.e**(0.05   * (V + 44))
    return numerator / denominator    
def alpha_f(V):
    import math
    numerator   = 0.012 * math.e**(-0.008 * (V + 28)) 
    denominator =   1   + math.e**(  0.15 * (V + 28))
    return numerator / denominator    
def beta_f(V):
    import math
    numerator   =0.0065 * math.e**(-0.02 * (V + 30)) 
    denominator =   1   + math.e**(-0.2  * (V + 30))
    return numerator / denominator
#########################____________________finish

#5c. I_K
def alpha_X(V):
    import math
    numerator   = 0.0005 * math.e**(0.083 * (V + 50))
    denominator = 1      + math.e**(0.057 * (V + 50))
    return numerator / denominator

def beta_X(V): 
    import math
    numerator   = 0.0013 * math.e**(-0.06 * (V + 20))
    denominator = 1      + math.e**(-0.04 * (V + 20))
    return numerator / denominator
def Xi(V):
    import math
    numerator   = math.e**(0.04 * (V + 77)) - 1
    denominator = (V + 77) * math.e**(0.04 * (V + 35))
    if V > -100:
        try:
            return 2.837 * numerator / denominator
        except ZeroDivisionError:
            return 0.6088832786043272
    elif V <= -100:
        return 1
#########################____________________finish
    
#5d. I_K1
def alpha_K1(V):
    import math
    denominator =1 + math.e**(0.2385 * (V - E_K1 - 59.215)) 
    return 1.02 / denominator
def beta_K1(V):
    import math
    component_1 = 0.49124 * math.e**(0.08032 * (V - E_K1 + 5.476))
    component_2 =           math.e**(0.06175 * (V - E_K1 - 594.31))
    denominator =    1    + math.e**(-0.5143 * (V - E_K1 + 4.753))
    return (component_1 + component_2) / denominator
#########################____________________finish

#5e. I_Kp
def Kp(V):
    import math
    denominator = 1+ math.e**((7.488 - V) / 5.98)
    return 1.0 / denominator
#########################____________________finish

#5f. I_b
def I_b(V):
    return gbar_b * (V - E_b)
#########################____________________finish


