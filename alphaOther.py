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
#demo=[45,1,10,5,10,2]
demo=[135,1,10,5,10,2]

#dispense user inputs: alpha, beta, lambda
alpha=demo[0]
beta_1=demo[1]
beta_1_prime=demo[2]
beta_2=demo[3]
beta_2_prime=demo[4]
lambda_1=demo[5]

#Declare and test r and q
def q(x):
    return -x
def r(x):
    #incl test to see if r>0 for all x
    return 1+2*x
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

#Calculate phi and psi
phi=math.cos(alpha)/math.sin(alpha)
print("phi=",phi)

psi=(beta_1+lambda_1*beta_1_prime)/(beta_2+lambda_1*beta_2_prime)
print("psi=",psi)

#integrand containing q in positive lambda case
def integrand_q(s):
    #integral used in calc of p_1_lambda_pos
    return q(s)*math.sinh(2*math.sqrt(lambda_1)*g(s)[0])

#calc zeta
def integrand_qzeta1(s):
    #integral used in calc of zeta
    return q(s)*math.cosh(2*math.sqrt(lambda_1)*g(s)[0])

def integrand_qzeta2(s):
    #integral used in calc of zeta
    return q(s)*math.sinh(2*math.sqrt(lambda_1)*g(s)[0])

a = phi + psi*math.cosh(2*math.sqrt(lambda_1)*g(1)[0]) - integrate.quad(integrand_qzeta1,0,1)[0]
b = psi*math.sinh(2*math.sqrt(lambda_1)*g(1)[0]) - integrate.quad(integrand_qzeta2,0,1)[0]
c = (-a)/b
zeta = c
print("zeta=",zeta)
zee=.5*arccoth(zeta)
print("zee=",zee)

#=========================================================================================================
# Case3&6 Case3&6 Case3&6 Case3&6 Case3&6 Case3&6 Case3&6 Case3&6 Case3&6 Case3&6 Case3&6 Case3&6 Case3&6 
#=========================================================================================================


def integrand_q36(s):
    #integral used in integrand of case 3 and 6 numerator
    return q(s)*math.sinh(2*math.sqrt(lambda_1)*g(s)[0]+2*zee)

def integrand_q63(s):
    #integral used in integrand of case 3 and 6 numerator
    return q(s)*math.cosh(2*math.sqrt(lambda_1)*g(s)[0]+2*zee)

#P3
def P3(x):
    a=1/(2*math.sqrt(lambda_1*r(x))*(math.cosh(math.sqrt(lambda_1)*g(x)[0]+zee)**2))
    b=-psi*math.sinh(2*zee)+integrate.quad(integrand_q36,0,x)[0]
    return a*b

def integrand_q3(s):
    return q(s)*(math.sinh(math.sqrt(lambda_1)*g(s)[0]+zee)**2)

def M_P3():
    a=math.tanh(math.sqrt(lambda_1)*g(1)[0]+zee)
    b=integrate.quad(integrand_q36,0,1)[0]-phi*math.sinh(zee)
    c=a*b
    d=-2*integrate.quad(integrand_q3,0,1)[0]+2*phi*(math.sinh(zee))**2
    e=c+d
    f=1/(2*lambda_1)
    h=f*e
    M_P3=str(round(h,3))
    return "M[p(x)]="+M_P3

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
    print(a)
    for x in intermediate_array_x:
        b = integrate.quad(integrand_q36,0,x)[0]
        truthValue = (b-a>0)
        if truthValue == True:
            continue
        elif truthValue == False:
            print(b-a)
            print("ERR: Sec 5 positivity condition failed at x near", x, ".")
            quit()
        else:
            print("ERR check truthValue of Sec 5.")
            quit()
    
    xcoord = intermediate_array_x
    ycoord = intermediate_array_y

    plt.plot(xcoord, ycoord)
    plt.xlabel("x in [0,1]")
    plt.ylabel("p(x)")
    plt.title(M_P3())
    plt.show()
    
    return print("Success: positivity condition Sec 5.")
    
#P6
def integrand_q6(s):
    return q(s)*(math.cosh(math.sqrt(lambda_1)*g(s)[0]+.5*arccoth(zeta))**2)

def P6(x):
    a=1/(2*math.sqrt(lambda_1*r(x))*(math.sinh(math.sqrt(lambda_1)*g(x)[0]+.5*(arccoth(zeta))))**2)
    b=-psi*math.sinh(arccoth(zeta))+integrate.quad(integrand_q36,0,x)[0]
    return a*b

def M_P6():
    a=coth(math.sqrt(lambda_1)*g(1)[0]+.5*arccoth(zeta))
    b=-integrate.quad(integrand_q36,0,1)[0]+phi*math.sinh(arccoth(zeta))
    c=a*b
    d=2*integrate.quad(integrand_q6,0,1)[0]-2*phi*(math.cosh(.5*arccoth(zeta)))**2
    e=c+d
    f=1/(2*lambda_1)
    h=f*e
    M_P6=str(round(h,3))
    return print("M[P6]="+M_P6)

def P6_positivity():
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
    intermediate_array_y = [P6(k) for k in intermediate_array_x]
    a=psi*math.sinh((arccoth(zeta)))
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
    plt.ylabel("P6(x)")
    plt.title(M_P6())
    plt.show()
    
    return print("Success: positivity condition Sec 5.")

def Compare_36():
    a=coth(2*math.sqrt(lambda_1)*g(1)[0]+arccoth(zeta))
    b=integrate.quad(integrand_q36,0,1)[0]-phi*math.sinh(arccoth(zeta))
    c=a*b
    d=-integrate.quad(integrand_q63,0,1)[0]+phi*math.cosh(arccoth(zeta))
    e=c+d
    comp=str(round(e,3))
    return print("sign=",comp)


#==============================================================================
#                              Run the program
#==============================================================================
P3_positivity()
#P6_positivity()




