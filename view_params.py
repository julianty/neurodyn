import numpy as np
import cPickle as  pickle

# Load path estimates and action curves
allpaths = np.load("neurodyn/paths.npy")
aerr = np.load("neurodyn/action_errors.npy")
params = np.load("neurodyn/params.npy")
with open('neurodyn/model_params.txt', "rb") as fp:
    model_params = pickle.load(fp)



def rmse(x,y, normalize=True): 
    """
    Computes the percent mse of two row vectors
    
    x: vector 1
    y: vector 2
    """
    x, y = (np.array(x), np.array(y))
    error = y - x
    if normalize == True:
        mask = x != 0
        error[:,mask] = np.divide(error[:,mask], x[mask])
    sq_error = np.square(error)
    mean_sq_error = sq_error.mean(axis=1)
    rmsq_error = np.sqrt(mean_sq_error)
    return rmsq_error

def param_error(model_params, params, breakdown=False):
    """
    Used to compute the RMSE of the two inputs.
    
    model_params: list or array of true model parameters
    params: list or array of parameters estimated using varanneal
    """
    #coerce params to np.arrays
    mp = np.array(model_params)
    
    #Thing messes up if you only give it one row, so here's a workaround
    try: 
        if params.shape[1] != 0:
            par = np.array(params)
    except IndexError: 
        par = np.array(params)[None,:]

    #0:3 are the conductances
    #3:6 are the reversal potentials
    #6:13 are the vBias values
    #13:20 are am[] values
    #20:27 are bm[] values
    #27:34 are ah[] values
    #34:41 are bh[] values
    #41:48 are an[] values
    #48:55 are bn[] values
    tot_rmse = []
    tot_rmse.append(rmse(mp[0], par[:,0:3]))
    tot_rmse.append(rmse(mp[1], par[:,3:6]))
    tot_rmse.append(rmse(mp[2], par[:,6:13]))
    tot_rmse.append(rmse(mp[3], par[:,13:20]))
    tot_rmse.append(rmse(mp[4], par[:,20:27]))
    tot_rmse.append(rmse(mp[5], par[:,27:34]))
    tot_rmse.append(rmse(mp[6], par[:,34:41]))
    tot_rmse.append(rmse(mp[7], par[:,41:48]))
    tot_rmse.append(rmse(mp[8], par[:,48:56]))

    tot_rmse = np.array(tot_rmse)
    if breakdown:
        return tot_rmse
    else:
        return tot_rmse.mean(axis=0)

last_beta = params.shape[0]
print("last beta was: %d" % last_beta)
print(param_error(model_params, params[-1,:], breakdown=True))
