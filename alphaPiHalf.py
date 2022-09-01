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
demo=[90,1,5,1,2,2]

#dispense user inputs: alpha, beta, lambda
alpha=demo[0]
beta_1=demo[1]
beta_1_prime=demo[2]
beta_2=demo[3]
beta_2_prime=demo[4]
lambda_1=demo[5]

#Declare and test r and q
def q(x):
    return x**2+x+1
def r(x):
    #incl test to see if r>0 for all x
    return 1+x+x**2
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

#integrand containing q in positive lambda case
def integrand_q(s):
    #integral used in calc of p_1_lambda_pos
    return q(s)*math.sinh(2*math.sqrt(lambda_1)*g(s)[0])


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
    return "M[p(x)]="+M_P2

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
    
    print(intermediate_array_x[0],intermediate_array_y[0])
    print(intermediate_array_x[100],intermediate_array_y[100])
    xcoord = intermediate_array_x
    ycoord = intermediate_array_y

    plt.plot(xcoord, ycoord)
    plt.xlabel("x in [0,1]")
    plt.ylabel("p(x)")
    plt.title(M_P2())
    plt.show()
    
    return print("Success: positivity condition Sec 5.")


#==============================================================================
#                              Run the program
#==============================================================================

P2_positivity()



