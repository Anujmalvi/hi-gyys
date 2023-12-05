#Curve fitting in POWER equation: y = alpha*x^beta
# Name: Tanveer Nazir Shaikh, AY 2023-24, Sem-I
# Batch: B3, Roll No.: 352053, GR No.: 22220099
import numpy as np
import math as m
def TNS_CFPE(x,y):
  x = x.astype(float)
  y = y.astype(float)
  X = np.log(x)
  Y = np.log(y)
  a = np.array([[len(x),sum(X)],
               [sum(X),sum(X*X)]])
  d = np.array([[sum(Y)],
               [sum(X*Y)]])
  b = np.linalg.solve(a,d)
  alpha = m.exp(b[0])
  beta = b[1]
  print("y = %.4f * x ^ %.4f"%(alpha,beta))
  TNS_CFPE(np.array([61,26,7,2.6]), np.array([350,400,50,600]))

# Guass elimination method
# Name: Tanveer Nazir Shaikh, TY Sem I: 2023-24
# Batch: B3, Roll No.: 352053
import numpy as np
def TNS_GEM(a,d):
  a = np.array(a, dtype = float)
  d = np.array(d, dtype = float)
  n = len(d)
  print("Number of equations = ",n)
  for i in range(0,n,1):
    for k in range(i+1,n,1):
      fact = a[k,i]/a[i,i]
      for j in range(0,n,1):
          a[k,j] = a[k,j] - fact * a[i,j]
          d[k] = d[k]-fact*d[i]
  print("Upper triangular matrix = ")
  print(a)
  print(d)
  x = np.zeros(n)
  for i in range(n-1,-1,-1):
    temp = 0
    for j in range(i+1,n,1):
      temp = temp + a[i,j]*x[j]
    x[i] = (d[i]-temp)/a[i,i]
  for i in range(n):
    print("x(%i) = %.4f"%(i,x[i]))
TNS_GEM(np.array([[4,1,2,3],[3,4,1,2],[2,3,4,1],[1,2,3,4]]), np.array([40,40,40,40]))

# Guass seidel method
# Name: Tanveer Nazir Shaikh, TY Sem I: 2023-24
# Batch: B3, Roll No.: 352053
import numpy as np
#def TNS_GSM(a,d):
a = np.array([[4,1,2,3],[3,4,1,2],[2,3,4,1],[1,2,3,4]])
d = np.array([40,40,40,40])
a = np.array(a, dtype = float)
d = np.array(d, dtype = float)
n = len(d)
maxitr = 6
print("Number of equations = ",n)
x = np.zeros(n)
for itr in range(maxitr):
  for i in range(0,n,1):
    temp = 0
    for j in range(0,n,1):
      if (i!=j):
        temp = temp + a[i,j]*x[j]
    x[i] = (d[i]-temp)/a[i,i]
for i in range(0,n,1):
  print(x[i])
Number of equations = 4

# Lagrange's Interpolation
# Name: Tanveer Nazir Shaikh, AY 2023-24, Sem-I
# Batch: B3, Roll No.: 352053, GR No.: 22220099
import numpy as np
def TNS_LI(x,y,xr):
  x=x.astype(float)
  y=y.astype(float)
  n=len(x)
  yr=0
  for i in range(n):
    L=1
    for j in range(n):
      if i!=j:
        L=L*(xr-x[j])/(x[i]-x[j])
  yr=yr+y[i]*L
print("y at x = %.4f is equal to %.4f" %(xr,yr))
TNS_LI(np.array([5,7,11,13,17]),np.array([150,392,1452,2366,5202]),9)


# Newton Forward Difference Interpolation
# Name: Tanveer Nazir Shaikh, AY 2023-24, Sem-I
# Batch: B3, Roll No.: 352053, GR No.: 22220099
import math as m
import numpy as np
def TNS_NFDI(x,y,xr):
  x=x.astype(float)
  y=y.astype(float)
  n=len(x)
  delta=np.zeros((n-1,n-1))
  for j in range(n-1):
    for i in range((n-1)-j):
      if j==0:
        delta[i,j]=y[i+1]-y[i]
      else:
        delta[i,j]=delta[i+1,j-1]-delta[i,j-1]
  h=x[1]-x[0]
  u=(xr-x[0])/h
  term=0
  mult=1
for j in range(n-1):
  mult=mult*(u-j)
  yerm=term+delta[0,j]/m.factorial(j+1)*mult
yr=y[0]+term
print("y at x = %.4f is equal to %.4f" %(xr,yr))
TNS_NFDI(np.array([100,150,200,250,300,350,400]),np.array([10.63,13.03,15.04,16.81,18.42,19.90,21.27]),16)

# Newton raphson method
# Name: Tanveer Nazir Shaikh, TY Sem I: 2023-24
# Batch: B3, Roll No.: 352053
import math as m
def TNS_NRM(f,df,ddf,x0,maxitr,acc):
  while abs(f(x0)*ddf(x0)/df(x0)**2)>1:
    print("The initial guess is incorrect, provide new initial guess")
    x0 = float(input("Enter the new value of initial guess"))
  for i in range(maxitr):
    x1 = x0-f(x0)/df(x0)
    if abs(x1-x0)>acc:
      break
    else:
      x1 = x0
  print("The root of the iteration is %.4f" %x1)
TNS_NRM(lambda x:m.exp(x)-m.sin(x),lambda x:m.exp(x)-m.cos(x),lambda x:m.exp(x)+m.


# Quadratic/Parabolic curve fitting
# Name: Tanveer Nazir Shaikh, AY 2023-24, Sem-I
# Batch: B3, Roll No.: 352053, GR No.: 22220099
import numpy as np
def TNS_PCF(x,y):
  x = x.astype(float)
  y = y.astype(float)
  a = np.array([[len(x),sum(x),sum(x*x)],
               [sum(x),sum(x*x),sum(x*x*x)],
               [sum(x*x),sum(x*x*x),sum(x*x*x*x)]])
  d = np.array([[sum(y)],
               [sum(x*y)],
               [sum(x*x*y)]])
  b = np.linalg.solve(a,d)
  print("y = %.4f+(%.4f)*x+(%.4f)*x^2"%(b[0],b[1],b[2]))
TNS_PCF(np.array([0,1,2,3,4]),np.array([1,1.8,1.3,2,6.3]))

# Runge Kutta 2 method
# Name: Tanveer Nazir Shaikh, AY 2023-24, Sem-I
# Batch: B3, Roll No.: 352053, GR No.: 22220099
import math as m
def TNS_RK2(fun,x0,y0,xn,n):
    h=(xn-x0)/n
    for i in range(1,n+1):
        k1 = h*fun(x0,y0)
        k2 = h*fun(x0+h,y0+k1)
        ynew = y0+1/2*(k1+k2)
        x0 = x0 + h
        y0 = ynew
        print("x = %.4f"%x0,"y = %.4f"%y0)
TNS_RK2(lambda x,y:2-x*y,5,2,5.1,10)


 Simpson's 1/3rd Method
# Name: Tanveer Nazir Shaikh, TY Sem I: 2023-24
# Batch: B3, Roll No.: 352053
import math as m
def TNS_S13(fun,x0,xn,n):
  h = (xn-x0)/n
  y0 = fun(x0)
  yn = fun(xn)
  yodd = 0
  yeven = 0
  for i in range(1,n,2):
    yodd = yodd + fun(x0+i*h)
  for j in range(2,n,2):
    yeven = yeven + fun(x0+j*h)
  A = (1/3)*h*(y0+yn+4*yodd+2*yeven)
  print("Intergation = ", A)
TNS_S13(lambda x: 1/(1+x*x),0,6,12)


# Simpson's 3/8th Method
# Name: Tanveer Nazir Shaikh, TY Sem I: 2023-24
# Batch: B3, Roll No.: 352053
import math as m
def TNS_S38(fun,x0,xn,n):
  h = (xn-x0)/n
  y0 = fun(x0)
  yn = fun(xn)
  yr = 0
  ym3 = 0
  for i in range(1,n):
    yr = yr + fun(x0+i*h)
  for j in range(3,n,3):
    ym3 = ym3 + fun(x0+j*h)
  A = (3*h/8)*(y0+yn+3*(yr-ym3)+2*ym3)
  print("Intergation = ", A)
TNS_S38(lambda x: 1/(1+x*x),0,6,12)


# Successive approximation method
# Name: Tanverr Nazir Shaikh, TY Sem I: 2023-24
# Batch: B3, Roll No.: 352053
import math as m
def TNS_SAM(g,gf,x0,maxitr,acc):
  while abs (gf(x0))>1:
    print("The initial guess is incorrect, provide new initial guess")
    x0 = float(input("Enter the new value of initial guess"))
  for i in range(maxitr):
    x1 = g(x0)
    if abs(x1-x0)>acc:
      break
    else:
      x1 = x0
  print("The root of the iteration is %.4f" %x1)
TNS_SAM(lambda x:(x**3+8)/15, lambda x:x**2/5,1,15,0.001)


# Straight line curve fitting
# Name: Tanveer Nazir Shaikh, AY 2023-24, Sem-I
# Batch: B3, Roll No.: 352053, GR No.: 22220099
import numpy as np
def TNS_SLCF(x,y):
  x = x.astype(float)
  y = y.astype(float)
  a = np.array([[len(x),sum(x)],
               [sum(x),sum(x*x)]])
  d = np.array([[sum(y)],
               [sum(x*y)]])
  b = np.linalg.solve(a,d)
  print("y = %.4f + (%.4f) * x"%(b[0],b[1]))
TNS_SLCF(np.array(([6,7,7,8,8,8,9,9,10])),np.array([5,5,4,5,4,3,4,3,3]))

TNS_SLCF(np.array(([10,20,30,40,50,60,70,80,90])),np.array([5,5,4,5,4,3,4,3,3]))

# Trapezoidal rule
# Name: Tanveer Nazir Shaikh, TY Sem I: 2023-24
# Batch: B3, Roll No.: 352053
import math as m
def TNS_TR(fun, x0, xn, n):
  h = (xn-x0)/n
  y0 = fun(x0)
  yn = fun(xn)
  yr = 0
  for i in range(1,n):
    yr = yr + fun(x0+i*h)
  A = 0.5*h*(y0+yn+2*yr)
  print("Integration = ", A)
TNS_TR(lambda x: 2*x-x**2,0,3,6)


# Bisection method
# Name: Tanveer Nazir Shaikh, TY Sem I: 2023-24
# Batch: B3, Roll No.: 352053
import math as m
def TNS_BSM(f,x1,x2,itr,acc):
  while f(x1)*f(x2) > 0:
    print("The root does not lies in between x1 and x2, change initial guesses")
    x1=float(input("enter the value of new x1="))
    x2=float(input("enter the value of new x2="))
  for i in range(1,itr):
    x0=(x1+x2)/2
    if f(x1)*f(x0) < 0:
      x1=x1
      x2=x0
    elif f(x1)*f(x0)>0:
      x1=x0
      x2=x2
    else:
      break;
    if abs (x1-x0) <=acc: # checking for accuracy
      break;
  print("the root is ",x0)
TNS_BSM(lambda x:m.sin(10*x)+m.cos(3*x),3,6,100,0.001)

TNS_BSM(lambda x:m.e**x-x*m.sin(x),3,6,100,0.001)

TNS_BSM(lambda x:x*m.sin(x)-1,3,6,100,0.001)

# Eulers Method
# Name: Tanveer Nazir Shaikhi, AY 2023-24, Sem-I
# Batch: B3, Roll No.: 352053, GR No.: 22220099
import math as m
def TNS_EM(fun,x0,y0,xn,n):
    h=(xn-x0)/n
    for i in range(1,n+1):
        ynew = y0+h*fun(x0,y0)
        x0 = x0 + h
        y0 = ynew
        print("x = %.4f"%x0,"y = %.4f"%y0)
TNS_EM(lambda x,y:m.sqrt(x+m.sqrt(y)),2,4,2.5,10)
