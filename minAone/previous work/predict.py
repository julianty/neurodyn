'''
Usage: python /path/to/predict_minAone.py [predict-time-steps] [PATH file num]

For example, to predict after the end of D7_M3_PATH1.dat, run this in the same directory:
python ~/pythonscripts/predict_minAone.py 10000 1

OUTPUT: Will save two files: data.dat with contains D+2 columns: [times, x0, x1, ..., xD, inj]. It will also save the parameter values in params.dat

_________________________________________________________

Prediction script for minAone. Should read in the numbers of states, parameters, time points from  specs.txt; read the equations from equations.txt; and load the parameter values from one of the path files, and do an integration.

Current file must have sufficient data AFTER the estimation to use for prediction. You can change the first line of specs.txt after the estimtaion, and the prediction will start from the appropriate point part-way through the estimation.

For external functions (ex. ghk(V) in myfunctions.cpp) you will have to rewrite the function (derivatives not necessary) in python file called "myfunctions.py"

myfunctions.py should start with "from sympy.functions import *" if it uses any mathematical functions (i.e. use sympy version of exp(x), not numpy or other package). Any "if" statements you want to include myfunctions.py need to be written using the Piecewise sympy function instead. 

Currently can only handle 1 input currenet. I think need to change eqn_wrapper() with multiple currents.

'''
import matplotlib
import numpy as np
matplotlib.use('Agg')
import sys
import scipy as sp
from scipy.integrate import odeint
import sympy as sym
import matplotlib.pyplot as plt
import os
from scipy import *

try:
    import myfunctions
except:
    print "Alert: no myfunctions.py file in directory"

if len(sys.argv) < 2:
    raise ValueError("Num time steps not specified")
elif len(sys.argv)<3:
    raise ValueError("Path not specified")    

predict_steps = int(sys.argv[1])
pathnum = sys.argv[2]




### Load numbers of states, measurments and equations:

with open('equations.txt') as eqfile:
    count = 0
    eqns = []
    eqnstrings = []
    states = []
    params = []
    inj = []
    funcs = []
    funcstrings = []
    for i,line in enumerate(eqfile):
        if line[0] == '#': continue
        count += 1
        # skip first line, problem name
        if count == 1:  continue
        if count == 2: 
            nS, nP, nU, nI,nF,nM = [int(x) for x in line.strip().split(',')]
            
        #equations start at line 3
        elif count <= 2 + nS:
           # print "eqns lines = ", line
            eqnstrings.append(line)
        # Variables names at line 4+nS
        elif 3+ nS< count < 4 + 2*nS:
            # print "states lines = ", line
            states.append(sym.symbols(line))
        # Parameters start at line 4+2*nS
        elif 3+2*nS < count < 4+2*nS+nP:
            # print "param lines = ", line
            params.append(sym.symbols(line))
        # Injected current starts at 5+2*nS+nP+nU+nM
        elif 4+2*nS+nP+nU+nM<= count < 4+2*nS+nP+nU+nM+nI:
            #print "Iinj lines = ", line
            inj.append(sym.symbols(line))
        elif 4+2*nS+nP+nU+nM+nI <= count < 4+2*nS+nP+nU+nM+nI+nF:
            #print "Fcn lines = ", line
            fcnname, varnum = line.strip().split(',')

            funcstrings.append(fcnname)
            try:
                fcntmp = eval('myfunctions.'+fcnname)
            except: 
                print fcnname
                ValueError("Is function defined in myfunctions.py?")

                
            funcs.append(fcntmp)

data_files = []
I_files = []
with open('specs.txt') as specsfile:
    count = 0
    for i,line in enumerate(specsfile):
        if line[0] == '#': continue
        count += 1
        # skip first line, problem name
        if count == 1:
            nT = 2*int(line.strip())+1
        elif count == 2:
            skip = int(line.strip())
        elif count == 3:            
            dt = float(line.strip())/2.0
        elif 4 < count <= 4 + nM:
            data_files.append(line.strip())
        elif count >= 4 + nM:
            Ifile = line.strip()
            I_files.append(Ifile)
    

# Load in current file, ignore first 'skip' number of time steps
try:
  I1 = sp.loadtxt(I_files[0])[skip:]
except:
  I1 = np.zeros(2*nT + predict_steps)  
  
try:

	I2 = sp.loadtxt(I_files[1])[skip:]
except:
	I2 = np.zeros(2*len(I1))
 
 
if len(I1) < nT + predict_steps:
    raise ValueError("Current file too short. Only {} steps available for prediction".format(len(I1)-nT-1))

# Use data.dat file if that exists (this script as been run before). Else, read in param values from a path file, and save last path as data.dat
#try:
#    path = sp.loadtxt('data.dat', usecols=sp.arange(1,nS+1))
#    param_values = sp.loadtxt('params.dat', usecols=[1],dtype=float)
#    import ipdb; ipdb.set_trace()
#    print "Using existing data.dat file and params.dat file"
#except:
# Above was written to save some time loading path file, but I don't think it takes long neough to invest the effort to make it work.

print "loading Path %d, saving to data.dat and params.dat" % int(pathnum) 
# Take the last nP values, and the last row (w/ largest beta)
pathfile = 'D{}_M{}_IC{}.dat'.format(nS,nM,pathnum)
last_path_data = sp.loadtxt(pathfile)
last_path_data = last_path_data[-1,:]
param_values = last_path_data[-nP:]
path = last_path_data[3:3+nT*nS]
sp.ndarray.reshape
path = path.reshape((nT,nS))

data_times = (sp.arange(0,nT,1))*dt   

sp.savetxt('data.dat', sp.column_stack((data_times,path,I1[:nT])))

param_array = sp.array(zip([x.name for x in params],param_values))
sp.savetxt('params{0}.dat'.format(pathnum), param_array, fmt='%s',delimiter=' \t')
 
print "Path saved"

funcNameSpace = dict(zip(funcstrings,funcs))
paramNameSpace = dict(zip(params,param_values))
    
x = sym.symbols('x:'+str(nS+nI+1))
for eq in eqnstrings:
    eqfunc = sym.sympify(eq,locals=funcNameSpace)
    eqfunc = eqfunc.subs(paramNameSpace)
    eqlamb = sym.lambdify(x, eqfunc.subs(zip(states+inj,x)))      
    eqns.append(eqlamb)

def eqns_wrapper(x,t,I,dt):
    deriv = []
    Iinj = I[0][int(sp.floor(t/dt))] 
    I_time = I[1][int(sp.floor(t/dt))]   
    input = list(x)+[Iinj,I_time]
    for eq in eqns:
        deriv.append( eq(*input) )        
    return deriv        

init = path[-1,:]

print "Integrating..."
time = sp.linspace(nT*dt,(nT+predict_steps)*dt, predict_steps)

predict = odeint(eqns_wrapper, init, time, ([I1,I2],dt))

data0 = sp.loadtxt(data_files[0])
data0 = data0[skip:]

def inverse_h(data,offset,Rgate,Iref,eIr):
	return (data + offset)/(Rgate)/(Iref*(1+eIr))
	

try:
	offset1 = float(param_values[-3])
	offset2 = float(param_values[-7])
	offset3 = float(param_values[-6])
	offset4 = float(param_values[-5])
	Iref = float(param_values[-10])
	RgateM = 0.00153
	RgateH = 0.00153 
	RgateN = 0.00153 
	eIrM = float(param_values[1])
	eIrH = float(param_values[2]) 
	eIrN = float(param_values[5]) 
	data1 = sp.loadtxt(data_files[1])
	#data1 = inverse_h(data1[skip:],offset2,RgateM,Iref,eIrM)
	data2 = sp.loadtxt(data_files[2])
	data2 = inverse_h(data2[skip:],offset3,RgateH,Iref,eIrH)
	data3 = sp.loadtxt(data_files[3])
	data3 = inverse_h(data3[skip:],offset4,RgateN,Iref,eIrN)
except:
	pass
#
if len(data0) < nT+predict_steps:
    print "Warning: data file does not have enough data points for entire prediction range"
    data_length = len(data0)
else:
    data_length = nT+predict_steps

plt.figure(figsize=(20,10))
plt.subplot(2,3,1)
plt.plot(data_times,path[:,0], color = 'blue', label = "Estimate")
plt.plot(np.hstack((data_times,time)), data0[0:data_length]-offset1, color = 'black', label="Data", alpha = 0.7)
plt.plot(time,predict[:,0],color = 'red', label = 'Prediction')
plt.ylabel("Voltage 1", fontsize = 20)
plt.legend(loc=3)
try:
	plt.subplot(2,3,2)
	plt.plot(data_times,path[:,4], color = 'blue', label = "Estimate")
	plt.plot(np.hstack((data_times,time)), data1[0:data_length]-offset1, color = 'black', label="Data", alpha = 0.7)
	plt.plot(time,predict[:,4],color = 'red', label = 'Prediction')
	plt.ylabel("Volatage 2", fontsize = 20)
	plt.legend(loc=3)
	plt.subplot(2,3,3)
	plt.plot(data_times,path[:,2], color = 'blue', label = "Estimate")
	plt.plot(np.hstack((data_times,time)), data2[0:data_length], color = 'black', label="Data", alpha = 0.7)
	plt.plot(time,predict[:,2],color = 'red', label = 'Prediction')
	plt.ylabel("activation", fontsize = 20)
	plt.legend(loc=3)
	plt.subplot(2,3,4)
	plt.plot(data_times,path[:,3], color = 'blue', label = "Estimate")
	plt.plot(np.hstack((data_times,time)), data3[0:data_length], color = 'black', label="Data", alpha = 0.7)
	plt.plot(time,predict[:,3],color = 'red', label = 'Prediction')
	plt.ylabel("activation", fontsize = 20)
	plt.legend(loc=3)
except:
	pass

plt.subplot(2,3,5)
plt.plot(np.hstack((data_times,time)),I1[0:data_length], color = 'purple', label = "Injected Current")
plt.xlabel("time (ms)", fontsize = 20)
plt.ylabel("Current (nA)", fontsize = 20)
plt.legend(loc=3)
plt.savefig('path{0}.png'.format(pathnum))


predict = sp.column_stack((time,predict, I1[nT:nT+predict_steps]))
sp.savetxt('predict{0}.dat'.format(pathnum), predict)

#sp.savetxt('data_extracted\data_est.dat',path[:,0])
#sp.savetxt('data_extracted\data_pre.dat',predict[:,0])        
#sp.savetxt('data_extracted\data_rec.dat',data0[0:data_length]-offset1)    
#
#sp.savetxt('data_extracted\I.dat',I1[0:data_length])