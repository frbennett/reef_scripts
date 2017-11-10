"""
*******************************************************************************
Global sensitivity analysis using a Sparse Random Sampling - High Dimensional 
Model Representation (HDMR) using the Group Method of Data Handling (GMDH) for 
parameter selection and Bayesian regression for parameter refinement
*******************************************************************************

author: 'Frederick Bennett'

"""
import pandas as pd
import math
import scipy.special as sp
from gmdhpy.gmdh import Regressor
from itertools import combinations
from sklearn.linear_model import Ridge
from sklearn.linear_model import Lasso
from sklearn.linear_model import ElasticNet
from sklearn.linear_model import RidgeCV 
from sklearn.linear_model import ARDRegression
from sklearn.linear_model import BayesianRidge
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
import matplotlib
from sklearn import metrics
from scipy.stats import linregress 
from numba import jit

class rshdmr():
    
    def __init__(self,data_file, poly_order=4, gmdh_ref_functions='linear_cov',**kwargs):
        self._seq_type='mode4_2'
        self._poly_order = poly_order 
        self._gmdh_ref_functions = 'linear_cov'
        self._admix_features = True
        self._alpha_ridge = 0.5
        self._alpha_lasso = 0.001
        self._epsilon = 0.01
        self._cutoff = 0.0001
        self._regression_type = 'bayesian'
        self._criterion_type='validate'
        for key, value in kwargs.items():
            setattr(self, "_"+key, value)
        self.read_data(data_file)

        
    def read_data(self,data_file):
        """
        dsdsd
        """
        if isinstance(data_file, pd.DataFrame):
            print(' found a dataframe')
            df = data_file
        if isinstance(data_file, str):
            df = pd.read_csv(data_file)
        self.Y = df['Y']
        self.X = df.drop('Y', axis=1)
        # we can clean up the original dataframe
        del df
        
    def shift_legendre(self,n,x):
        funct = math.sqrt(2*n+1) * sp.eval_sh_legendre(n,x)
        return funct
    
        
    def transform_data(self):
        self.X_T = pd.DataFrame()
        self.ranges = {}
        feature_names = list(self.X.columns.values)
        print(feature_names)
        for column in feature_names:
            max = self.X[column].max()
            min = self.X[column].min()
            print(column + " : min " + str(min) + " max " + str(max)) 
            self.X_T[column] = (self.X[column] - min) / (max-min)
            self.ranges[column] = [min,max]
        
            
    def legendre_expand(self):
        self.primitive_variables = []
        self.X_T_L = pd.DataFrame()
        for column in self.X_T:
            for n in range (1,self._poly_order+1):
                self.primitive_variables.append(column)
                column_heading = column + "_" + str(n)
                self.X_T_L[column_heading] = [self.shift_legendre(n, x) for x in self.X_T[column]]
        self.exp_feature_names = list(self.X_T_L.columns.values) 
        
    def gmdh_regression(self):
        self.gmdh_model = Regressor(ref_functions=(self._gmdh_ref_functions),
                      criterion_type= self._criterion_type,
                      feature_names=self.exp_feature_names,
                      criterion_minimum_width=5,
                      stop_train_epsilon_condition=self._epsilon,
                      layer_err_criterion='top',
                      l2=0.5,
                      seq_type= self._seq_type , 
                      max_layer_count= 30,
                      normalize=True,
                      admix_features = self._admix_features, 
                      n_jobs='max')
        self.gmdh_model.fit(self.X_T_L, self.Y)
        unselected_indicies = self.gmdh_model.get_unselected_features().split(", ")
        #for index in unselected_indicies:
            #del self.X_T_L[index]
            
        # calculate terms for variance
        self.variance_dict = {}
        selected_indices = self.gmdh_model.get_selected_features_indices()
        feature_count = len(self.exp_feature_names)
        prediction_coefficients = [0.0] * feature_count
        X0 = self.gmdh_model.predict([prediction_coefficients])[0]
        #print(X0)
        for index in selected_indices:
            prediction_coefficients[index] = 1.0
            variance = self.gmdh_model.predict([prediction_coefficients])[0] - X0
            self.variance_dict[self.exp_feature_names[index]]=[self.exp_feature_names[index], self.primitive_variables[index], variance]
            prediction_coefficients[index] = 0.0
            
        # calculate terms for covariance
        self.covariance_dict = {}
        for combo in combinations(selected_indices, 2): 
            #print(str(combo[0]) + ", " + str(combo[1]) + \
            #" : " + self.exp_feature_names[combo[0]] + ", " + self.exp_feature_names[combo[1]] + \
            #" : " + self.primitive_variables[combo[0]] + ", " + self.primitive_variables[combo[1]])
            if (self.primitive_variables[combo[0]] != self.primitive_variables[combo[1]]):
                X1akey = self.exp_feature_names[combo[0]]
                X1bkey = self.exp_feature_names[combo[1]]
                X1 = self.variance_dict[X1akey][2] + self.variance_dict[X1bkey][2]
                prediction_coefficients[combo[0]] = 1.0
                prediction_coefficients[combo[1]] = 1.0
                covariance = self.gmdh_model.predict([prediction_coefficients])[0] - X0 -X1
                #print("include this")
                key = self.exp_feature_names[combo[0]] + "_" + self.exp_feature_names[combo[1]]
                self.covariance_dict[key] = [self.exp_feature_names[combo[0]], self.exp_feature_names[combo[1]], \
                self.primitive_variables[combo[0]], self.primitive_variables[combo[1]], covariance]
                prediction_coefficients[combo[0]] = 0.0
                prediction_coefficients[combo[1]] = 0.0
                
    def select_features(self, cutoff):
        self.selected_features_dict = {}
        #sum up variances
        total_variance = 0
        for key in self.variance_dict:
            variance = self.variance_dict[key][2]
            total_variance += variance*variance
        for key in self.covariance_dict:
            variance = self.covariance_dict[key][4]
            total_variance += variance*variance
            
        for key in self.variance_dict:
            variance = self.variance_dict[key][2] * self.variance_dict[key][2] / total_variance * 100
            if (variance > cutoff):
                self.selected_features_dict[key] = self.variance_dict[key]
            
        for key in self.covariance_dict:
            variance = self.covariance_dict[key][4] * self.covariance_dict[key][4] / total_variance *100
            if (variance > cutoff):
                self.selected_features_dict[key] = self.covariance_dict[key]
    
 
    def ridge_regression(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, "_"+key, value)
        self.ridge_set = None
        self.ridge_set = pd.DataFrame()
        #build dataset
        for key in self.selected_features_dict:
            #covariance terms
            if (len(self.selected_features_dict[key])==5):
                variable_1 = self.selected_features_dict[key][0]
                variable_2 = self.selected_features_dict[key][1]
                self.ridge_set[key] = self.X_T_L[variable_1] * self.X_T_L[variable_2]

            #variance terms
            if (len(self.selected_features_dict[key])==3):
                variable_1 = self.selected_features_dict[key][0]
                self.ridge_set[key] = self.X_T_L[variable_1] 
            
            
            if self._regression_type == 'bayesian' :
                self.ridgereg = BayesianRidge()
            elif self._regression_type == 'lasso' :
                self.ridgereg = Lasso(alpha=self._alpha_lasso,normalize=True, max_iter=1e5) 
            elif self._regression_type == 'ridge' :
                self.ridgereg = Ridge(normalize=True,alpha=self._alpha_ridge) 
            elif self._regression_type == 'ard' :
                self.ridgereg = ARDRegression() 
            else :
                self.ridgereg = BayesianRidge()
                
                
            #self.ridgereg = Lasso(alpha=self._alpha_lasso,normalize=True, max_iter=1e5)
            #self.ridgereg.fit(self.ridge_set,self.Y) 
            #self.ridgereg = RidgeCV(alphas=[0.001, 0.005, 0.01, 0.05, 0.075, 0.1, 0.25, 0.5, 1.0, 5.0, 10.0, 15.0, 20.0])
            #self.myfit = self.ridgereg.fit(self.ridge_set,self.Y) 
            #self.ridgereg = BayesianRidge()
            #self.ridgereg = ARDRegression()
            #self.ridgereg = Ridge(normalize=True,alpha=self._alpha_ridge)
            #self.ridgereg.fit(self.ridge_set,self.Y)
            #self.ridgereg = ElasticNet(alpha=0.1, l1_ratio=0.7)
            #self.ridgereg.fit(self.ridge_set,self.Y)
            self.ridgereg.fit(self.ridge_set,self.Y)
            ridge_feature_names = list(self.ridge_set.columns.values)
            # now add new coefficients to dictionary
            self.ridge_coeffs = dict(zip(ridge_feature_names, self.ridgereg.coef_))
            # add the assignments to the dict in terms of the primitive variables
            for key in self.ridge_coeffs:
                if (len(self.selected_features_dict[key])==5):
                    variable_1 = self.selected_features_dict[key][2]
                    variable_2 = self.selected_features_dict[key][3]
                    ref_value = variable_1 + "_" + variable_2
                    self.ridge_coeffs[key] = [ref_value, self.ridge_coeffs[key]]
                
                if (len(self.selected_features_dict[key])==3):
                    variable_1 = self.selected_features_dict[key][1]
                    ref_value = variable_1
                    
                    self.ridge_coeffs[key] = [ref_value, self.ridge_coeffs[key]]
                    
    def eval_sobol_indices(self):
        self.sobol_indices = {}
        # finally evaluate the sobol indices from the ridge coefficients
        total_variance = 0
        for key in self.ridge_coeffs:
            total_variance += self.ridge_coeffs[key][1] * self.ridge_coeffs[key][1]
            
        for key in self.ridge_coeffs:
            variance = self.ridge_coeffs[key][1] * self.ridge_coeffs[key][1] / total_variance
            sobol_key = self.ridge_coeffs[key][0]
            if sobol_key in self.sobol_indices.keys():                
                self.sobol_indices[sobol_key] += variance
            else:
                self.sobol_indices[sobol_key] = variance
                
    def evaluate_func(self,X):
        sum = self.ridgereg.intercept_
        primitives = list(self.X.columns.values)
        X_expanded ={}
        for i in range(0, len(X)):
            # Transform input
            min = self.ranges[primitives[i]][0]
            max = self.ranges[primitives[i]][1]
            X_T = (X[i] - min) / (max-min)
            for j in range(1, self._poly_order+1):
                label = primitives[i] +'_' + str(j)
                legendre = self.shift_legendre(j,X_T)
                X_expanded[label] = [legendre]
    
        for key in self.ridge_coeffs:
            gmdh_coeff = self.selected_features_dict[key]
            ridge_coeff = self.ridge_coeffs[key][1]
            if len(gmdh_coeff)==3:    
                variable_term = X_expanded[gmdh_coeff[0]][0]
                sum += variable_term * ridge_coeff
            else:
                variable_term = X_expanded[gmdh_coeff[0]][0] * X_expanded[gmdh_coeff[1]][0]
                sum += variable_term * ridge_coeff         
        return sum
    
    def plot_hdmr(self):
        y_pred = self.ridgereg.predict(self.ridge_set)  
        matplotlib.pyplot.scatter(self.Y,y_pred)
        matplotlib.pyplot.ylabel('Predicted')
        matplotlib.pyplot.xlabel('Experimental')
        matplotlib.pyplot.show()
        
    def stats(self):
        y_pred = self.ridgereg.predict(self.ridge_set) 
        mse = metrics.mean_squared_error(y_pred,self.Y)
        mae = metrics.mean_absolute_error(y_pred,self.Y)
        evs = metrics.explained_variance_score(y_pred,self.Y)
        slope, intercept, r_value, p_value, std_err = linregress(self.Y, y_pred)
        print("mae error on test set   : {mae:0.3f}".format(mae=mae))
        print("mse error on test set   : {mse:0.3f}".format(mse=mse))
        print("explained variance score: {evs:0.3f}".format(evs=evs))
        print("===============================")
        print("slope     : ", slope)
        print("r value   : ", r_value)
        print("r^2       : ", r_value*r_value)
        print("p value   : ", p_value)
        print("std error : ", std_err)
        
    def print_sobol_indices(self):
        self.eval_sobol_indices()
        for key in self.sobol_indices:
            if self.sobol_indices[key] > 0.0001:
                print(key + " : " + str(self.sobol_indices[key]))
                
    def auto(self):
        self.transform_data()
        self.legendre_expand()
        print('====================================')
        self.gmdh_regression()
        print('====================================')
        self.select_features(self._cutoff)
        self.ridge_regression()
        

