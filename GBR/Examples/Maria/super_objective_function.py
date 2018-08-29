"""
Super Objective Function Module (SDEB)

Calculation of the square-root daily, exceedance and bias (SDEB) objective function used in current GBR Sacramento calibrations and Queensland Hydrology Group.

Objective function combines multiple terms to try to simultaneously achieve good fits in overall bias, daily flows, and daily exceedance statistics. It is well suited for rainfall-funoff model calibration and doesn't have any free parameters.

Module contained is designed to call on pairs of DataFrames (observed and predicted) with equivalent columns on a daily timestep. 
"""


"""
Coefficients
The two coefficients, lamda and alpha are used to balance the three terms within the objective function. 
If lamda is less than 1, the power transform has the effect of reducing the weight of errors on high flows, where the flow data are known to be less accurate.

Queensland Hydrology compared values of lamda randing from 0.05 to 1, and found that lamda=0.5 was the best compromise between high and low flow performance.

The weighting factor alpha was used to reduce the impact of timing errors on the objective function. This error can have a significant effect on the daily term where a misalignment of observed and modelled (pred) peak flow timing can result in large amplituded errors.
"""

lamda = 0.5
alpha = 0.1


def intersect(obs,pred):
    """
    Return the input pair of dataframes (obs,pred) with a common index made up of the intersection of
    the input indexes
    """
    if hasattr(obs,'intersect'):
        return obs.intersect(pred)
    idx = obs.index.intersection(pred.index)
    return obs.ix[idx],pred.ix[idx]

def bias_term(obs,pred):
    """
    The relative simulation bias
    """
    obs,pred = intersect(obs,pred)
    top = (obs-pred).sum()
    bottom = obs.sum()
    return (1 + abs(top/bottom))

def daily_term(obs,pred):
    """
    The sum of squared errors on power transformation of flow
    """
    obs,pred = intersect(obs,pred)
    square_error_transform = ((obs**lamda-pred**lamda)**2)
    return square_error_transform.sum()

def exceedance_term(obs,pred):
    """
    The sum of squared errors on power transformation of sorted observed and predicted flow series
    """
    sorted_obs = obs[:]
    sorted_pred = pred[:]
    sorted_obs.sort()
    sorted_pred.sort()
    sorted_square_error_transform = ((sorted_obs**lamda-sorted_pred**lamda)**2)
    return sorted_square_error_transform.sum()


"""
Super Objective Function which combines the three terms above
"""
def SDEB(daily_term,exceedance_term,bias_term):
    return (alpha*daily_term + (1-alpha)*exceedance_term)*bias_term
    