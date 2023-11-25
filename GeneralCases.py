import math
#import scipy
import scipy.integrate as integrate
import scipy.special as special
import sys
import matplotlib.pyplot as plt
import numpy as np


#function to quit out of a loop if check fails
def quit():
    sys.exit()
    return

#homebrew coth function
def coth(x):
    return 1/math.tanh(x)
def arccoth(x):
    return math.atanh(1/x)

#==============================================================================
#                  User input and universal tests on those inputs
#                  Tests performed here apply to all cases below
#==============================================================================

#store user inputs as an array of the form
#demo=[alpha,beta_1,beta_1_prime,beta_2,beta_2_prime,lambda_1]
demo=[0,5,11,1,1.5,2]
#demo=[.25*np.pi,5,10,1,1.5,1]
#demo=[.5*np.pi,5,10,1,1.5,.5]
#demo=[.75*np.pi,1,10,5,10,2]

#dispense user inputs: alpha, beta, lambda
alpha=demo[0]
beta_1=demo[1]
beta_1_prime=demo[2]
beta_2=demo[3]
beta_2_prime=demo[4]
lambda_1=demo[5]

#Declare and test r and q
def q(x):
    return math.cos(x) #x**2+x+1 #1/(math.cosh(x))**2 #x**2+x+1 #-100*x
def r(x):
    #incl test to see if r>0 for all x
    return 1+2*x #1+x+x**2
def sqrt_r(x):
    return math.sqrt(r(x))
def g(x):
    return integrate.quad(sqrt_r,0,x)


#check Delta_condition on Beta_j and Beta_j_prime 
delta_condition=(beta_1_prime*beta_2-beta_2_prime*beta_1>0)
if delta_condition == True:
    print("Success: delta condition passed.")
elif delta_condition == False:
    print("ERR: check delta condition.")
    quit()
else:
    print("ERR: check delta switch.")


#Determine Case 1&4, 2&5, or 3&6 based on alpha value
Case=0
if alpha == 0:
    Case = 1
    print("Success: Case 1")
elif alpha == .5*np.pi:
    Case = 2
    print("Success: Case 2")
elif alpha >= np.pi:
    print("ERR: alpha must be within [0,pi)")
    quit()
elif alpha < 0:
    print("ERR: alpha must be within [0,pi)")
    quit()
else:
    Case = 3
    print("Success: Case 3")


#Calculate phi and psi
phi=0
if Case == 1:
    print("phi= n/a in this case")
elif Case == 2:
    print("phi= n/a in this case")
elif Case == 3:
    phi=math.cos(alpha)/math.sin(alpha)
    print("phi=",phi)
else:
    print("ERR: Check phi switch")
    quit()

psi=(beta_1+lambda_1*beta_1_prime)/(beta_2+lambda_1*beta_2_prime)
print("psi=",psi)

#integrand containing q in positive lambda case
def integrand_q(s):
    #integral used in calc of p_1_lambda_pos
    return q(s)*math.sinh(2*math.sqrt(lambda_1)*g(s)[0])



#=========================================================================================================
# Case1&4 Case1&4 Case1&4 Case1&4 Case1&4 Case1&4 Case1&4 Case1&4 Case1&4 Case1&4 Case1&4 Case1&4 Case1&4
#=========================================================================================================
def P1_lambda_pos(x):
    a=1/(2*math.sqrt(lambda_1*r(x)))
    b=1/(math.cosh(math.sqrt(lambda_1)*g(x)[0])**2)
    c=psi*math.sinh(2*math.sqrt(lambda_1)*g(1)[0])-integrate.quad(integrand_q,x,1)[0]
    return a*b*c

def M_P1_pos_integrand(s):
    return q(s)*(math.sinh(math.sqrt(lambda_1)*g(s)[0])**2)

def M_P1():
    a=1/lambda_1
    b=psi*(math.sinh(math.sqrt(lambda_1)*g(1)[0])**2)
    c=integrate.quad(M_P1_pos_integrand,0,1)[0]
    M_P1=str(round(a*(b-c),3))
    return "M[P1]="+M_P1

def P1_positivity():
    if psi <= 0:
        print ("ERR: psi must be positive")
        quit()
    intermediate_array_x = [(1/100)*k for k in range(0,101)]
    intermediate_array_y = [P1_lambda_pos(k) for k in intermediate_array_x]
    a=-psi*math.sinh(2*math.sqrt(lambda_1)*g(1)[0])+integrate.quad(integrand_q,0,1)[0]
    for x in intermediate_array_x:
        b = integrate.quad(integrand_q,0,x)[0]
        truthValue = (b>a)
        if truthValue == True:
            continue
        elif truthValue == False:
            print("ERR: Sec 5 positivity condition failed at x near", x, ".")
            quit()
        else:
            print("ERR check truthValue of Sec 5.")
            quit()
    
    xcoord = intermediate_array_x
    ycoord = intermediate_array_y

    plt.plot(xcoord, ycoord)
    plt.xlabel("x in [0,1]")
    plt.ylabel("P1(x)")
    plt.title(M_P1())
    plt.show()
    
    return print("Success: positivity condition Sec 5.")
    

#=========================================================================================================
# Case2&5 Case2&5 Case2&5 Case2&5 Case2&5 Case2&5 Case2&5 Case2&5 Case2&5 Case2&5 Case2&5 Case2&5 Case2&5
#=========================================================================================================
def P2(x):
    a=1/(2*math.sqrt(lambda_1*r(x)))
    b=1/(math.sinh(math.sqrt(lambda_1)*g(x)[0])**2)
    c=integrate.quad(integrand_q,0,x)[0]
    return a*b*c

def M_P2_pos_integrand(s):
    return q(s)*(math.cosh(math.sqrt(lambda_1)*g(s)[0])**2)

def M_P2():
    a=1/(2*lambda_1)
    b=-coth(math.sqrt(lambda_1)*g(1)[0])*(integrate.quad(integrand_q,0,1)[0])
    c=2*integrate.quad(M_P2_pos_integrand,0,1)[0]
    M_P2=str(round(a*(b+c),3))
    return "M[P2]="+M_P2

def P2_positivity():
    intermediate_array_x = [(1/100)*k for k in range(0,101)]
    intermediate_array_y = [(1/100)*k for k in range(0,101)]
    intermediate_array_y[0]=q(0)/(2*lambda_1*r(0))
    for k in range(1,101):
        intermediate_array_y[k] = P2(intermediate_array_x[k])
    for x in intermediate_array_x[1:]:
        b = integrate.quad(integrand_q,0,x)[0]
        truthValue = (b>0)
        if truthValue == True:
            continue
        elif truthValue == False:
            #print("ERR: Sec 5 positivity condition failed at x near", x, ".")
            #quit()
            continue
        else:
            print("ERR check truthValue of Sec 5.")
            quit()
    
    #print(intermediate_array_x[0],intermediate_array_y[0])
    #print(intermediate_array_x[100],intermediate_array_y[100])
    xcoord = intermediate_array_x
    ycoord = intermediate_array_y

    plt.plot(xcoord, ycoord)
    plt.xlabel("x in [0,1]")
    plt.ylabel("P2(x)")
    plt.title(M_P2())
    plt.show()
    
    return print("Success: positivity condition Sec 5.")


#=========================================================================================================
# Case3&6 Case3&6 Case3&6 Case3&6 Case3&6 Case3&6 Case3&6 Case3&6 Case3&6 Case3&6 Case3&6 Case3&6 Case3&6 
#=========================================================================================================
#lambda_1 positive 
def integrand_qzeta1(s):
    #integral used in calc of zeta
    return q(s)*math.cosh(2*math.sqrt(lambda_1)*g(s)[0])

def integrand_qzeta2(s):
    #integral used in calc of zeta
    return q(s)*math.sinh(2*math.sqrt(lambda_1)*g(s)[0])

zeta=0
if Case == 3:
    a = phi + psi*math.cosh(2*math.sqrt(lambda_1)*g(1)[0]) - integrate.quad(integrand_qzeta1,0,1)[0]
    b = psi*math.sinh(2*math.sqrt(lambda_1)*g(1)[0]) - integrate.quad(integrand_qzeta2,0,1)[0]
    c = (-a)/b
    zeta = c
    print("zeta=",zeta)

def integrand_q36(s):
    #integral used in integrand of case 3 and 6 numerator
    return q(s)*math.sinh(2*math.sqrt(lambda_1)*g(s)[0]+arccoth(zeta))

#P3
def P3(x):
    a=1/(2*math.sqrt(lambda_1*r(x))*(math.cosh(math.sqrt(lambda_1)*g(x)[0]+.5*(math.atanh(1/zeta))))**2)
    b=-psi*math.sinh(arccoth(zeta))+integrate.quad(integrand_q36,0,x)[0]
    return a*b

def integrand_q3(s):
    return q(s)*(math.sinh(math.sqrt(lambda_1)*g(s)[0]+.5*math.atanh(1/zeta))**2)

def M_P3():
    a=math.sinh(math.sqrt(lambda_1)*g(1)[0])*(1/math.cosh(math.sqrt(lambda_1)*g(1)[0]+.5*arccoth(zeta)))
    b=2*lambda_1*math.cosh(.5*(arccoth(zeta)))
    c=a/b
    d=integrate.quad(integrand_q36,0,1)[0]-phi*math.sinh(arccoth(zeta))
    e=d*c
    f=-2*integrate.quad(integrand_q3,0,1)[0]
    M_P3=str(round(e+f,3))
    return "M[P3]="+M_P3

def P3_positivity():
    if psi > 0:
        if zeta < -1:
            print ("Success: psi pos and zeta < -1")
        else:
            print ("ERR: zeta must be <-1")
            quit()
    elif psi < 0:
        if zeta > 1:
            print ("Success: psi neg and zeta > 1")
        else:
            print ("ERR: zeta must be > 1")
            quit()
    intermediate_array_x = [(1/100)*k for k in range(0,101)]
    intermediate_array_y = [P3(k) for k in intermediate_array_x]
    a=psi*math.sinh(arccoth(zeta))
    for x in intermediate_array_x:
        b = integrate.quad(integrand_q36,0,x)[0]
        truthValue = (b>a)
        if truthValue == True:
            continue
        elif truthValue == False:
            print("ERR: Sec 5 positivity condition failed at x near", x, ".")
            quit()
        else:
            print("ERR check truthValue of Sec 5.")
            quit()
    
    xcoord = intermediate_array_x
    ycoord = intermediate_array_y

    plt.plot(xcoord, ycoord)
    plt.xlabel("x in [0,1]")
    plt.ylabel("P3(x)")
    plt.title(M_P3())
    plt.show()
    
    return print("Success: positivity condition Sec 5.")
    

#==============================================================================
#                              Run the program
#==============================================================================

if Case == 1:
    if lambda_1 > 0:
        P1_positivity()
    else:
        print("ERR: check lambda switch")
        quit()
elif Case == 2:
    if lambda_1 > 0:
        P2_positivity()
    else:
        print("ERR: check lambda switch")
        quit()
elif Case == 3:
    P3_positivity()
    print("")
else:
   print("ERR: Check Pk switch")
   quit()