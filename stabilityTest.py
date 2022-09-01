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
demo=[0,0,1,1,0,2]

#dispense user inputs: alpha, beta, lambda
alpha=demo[0]
beta_1=demo[1]
beta_1_prime=demo[2]
beta_2=demo[3]
beta_2_prime=demo[4]
lambda_1=demo[5]

#Declare and test r and q
def q(x):
    return -x**2
def r(x):
    #incl test to see if r>0 for all x
    return 1
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

psi=(beta_1+lambda_1*beta_1_prime)/(beta_2+lambda_1*beta_2_prime)
print("psi=",psi)

#integrand containing q in positive lambda case
def integrand_q(s):
    #integral used in calc of p_1
    return q(s)*math.sinh(2*math.sqrt(lambda_1)*g(s)[0])


#=========================================================================================================
# Case1&4 Case1&4 Case1&4 Case1&4 Case1&4 Case1&4 Case1&4 Case1&4 Case1&4 Case1&4 Case1&4 Case1&4 Case1&4
#=========================================================================================================
def P1(x):
    a=1/(2*math.sqrt(lambda_1*r(x)))
    b=1/(math.cosh(math.sqrt(lambda_1)*g(x)[0])**2)
    c=psi*math.sinh(2*math.sqrt(lambda_1)*g(1)[0])-integrate.quad(integrand_q,x,1)[0]
    return a*b*c

intermediate_array_psi = [.5+(1/10)*k for k in range(0,11)]
#intermediate_array_psi = [(1/5)*k for k in range(-5,6)]
#intermediate_array_MP1 = [(1/10)*k for k in range(0,11)]
intermediate_array_MP1 = [k for k in range(0,12)]

def P1_mod(x,n):
    a=1/(2*math.sqrt(lambda_1*r(x)))
    b=1/(math.cosh(math.sqrt(lambda_1)*g(x)[0])**2)
    c=intermediate_array_psi[n]*psi*math.sinh(2*math.sqrt(lambda_1)*g(1)[0])-integrate.quad(integrand_q,x,1)[0]
    return a*b*c

def MP1_integrand(s):
    return q(s)*(math.sinh(math.sqrt(lambda_1)*g(s)[0])**2)

def MP1():
    a=1/lambda_1
    b=psi*(math.sinh(math.sqrt(lambda_1)*g(1)[0])**2)
    c=integrate.quad(MP1_integrand,0,1)[0]
    M_P1=str(round(a*(b-c),3))
    return "M[p(x)]="+M_P1

def MP1_mod(n):
    a=1/lambda_1
    b=n*psi*(math.sinh(math.sqrt(lambda_1)*g(1)[0])**2)
    c=integrate.quad(MP1_integrand,0,1)[0]
    M_P1=(round(a*(b-c),3))
    return M_P1

def P1_graph():
    intermediate_array_x = [(1/100)*k for k in range(0,101)]
    
    for n in range(0,11):
        intermediate_array_y = [P1_mod(k,n) for k in intermediate_array_x]
        xcoord = intermediate_array_x
        ycoord = intermediate_array_y
        plt.plot(xcoord, ycoord)
        intermediate_array_MP1[n]=MP1_mod(intermediate_array_psi[n])
    #print(intermediate_array_x[0],intermediate_array_y[0])
    #print(intermediate_array_x[100],intermediate_array_y[100])

    M_P1=str(np.ptp(intermediate_array_MP1))
    print(max(intermediate_array_MP1))
    print(min(intermediate_array_MP1))
    
    plt.xlabel("x")
    plt.ylabel("p(x)")
    plt.title("Range of MP1="+M_P1)
    plt.show()
    
    return print("Success")
    

#==============================================================================
#                              Run the program
#==============================================================================

P1_graph()

