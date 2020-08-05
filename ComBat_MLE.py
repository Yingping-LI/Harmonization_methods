from sklearn.exceptions import ConvergenceWarning
import warnings
import pandas as pd
import numpy as np



    
class ComBat_MLE(object):
    """
    ComBat_MLE: A model used to remove scanner effects from observed data.
    
    
    Attributes
    -------------
    
    gamma_s_dict: dict type,  consisting of S items corresponding to different scanner settings.
             key: scanner setting label s;  value:the addictive scanner effects of setting s, array of shape=(n_features,1). 
    
    Delta_s_dict: dict type,  consisting of S items corresponding to different scanner settings.
             key: scanner setting label s;  value:the multiplicative scanner effects of setting s, array of shape =(n_features,n_features). 
    
    mu: mean of clear data, array of shape=(n_features,1).

    Sigma: variance matrix of clear data, array of shape=(n_features,n_features). 
    
    """
    
    
    def __init__(self):
        pass
        
        
   
    def _vector_to_diagonalmatrix(self, v):
        """
        Convert a vector v to a matrix whose diagnal value is v.

        Parameter
        ----------
        v: array-like, shape=(n, 1)

        Return
        ---------
        M: array-like, shape=(n,n)
        """

        n=v.shape[0]
        M=np.zeros([n,n])
        for i in range(n):
            M[i,i]=v[i]

        return M


    
    def _calculate_inverse(self, M, reg_covar=1e-6):
        """
        Calculate the inverse of a given matrix M.

        Parameter
        ----------
        M: array-like matrix.
        
        reg_covar: float, defaults to 1e-6.  Non-negative regularization added to the diagonal of covariance.


        Return
        ---------
        inv_M: inverse of matrix M.

        """

        try:
            inv_M=np.linalg.inv(M)
        except:
            n_col=M.shape[1]
            M.flat[::n_col+1] += reg_covar
            inv_M=np.linalg.inv(M)

        return inv_M
    

    
    def _get_data_shape(self, X_s_dict):
        """Compute the total data number and feature dimensionality.
        
        Parameters
        -----------
        X_s_dict: dict type,  consisting of feature datas from different scanner settings. 
                  key: scanner setting label s, value: feature data corresponding to scanner setting s.
              
              
       Returns
       -----------
       n_samples: int, the total data number
       
       n_features: int, dimensionality of the feature.
       
        """
        
        n_samples=0
        for s, X_s in X_s_dict.items():
            N_s,n_features=X_s.shape
            n_samples+=N_s
            
        return n_samples, n_features
    
    
    def _check_input_data(self, data, setting_labels):
        """Check the input data and its corresponding setting labels.
        
        Parameters
        ----------
        data: feature data to be harmonized, array of shape=(n_samples,n_features).
       
        setting_labels: scanner setting labels corresponding to "data",  array of shape=(n_samples,).
        """
        
        # Check the shape of data and setting_labels.
        n_data_samples,_=data.shape
        n_label_samples=setting_labels.shape[0]
        if n_data_samples!=n_label_samples:
            raise ValueError("The data and the setting labels should have the same number of data samples!")
   
            
            
        
    def _estimate_model_parameters(self, X_s_dict):
        """estimate the model parameters.
        
        Parameters
        ---------------
        X_s_dict: dict type,  consisting of feature datas from different scanner settings. 
                  key: scanner setting label s, value: feature data corresponding to scanner setting s.
                    
        """
        
        N, n_features=self._get_data_shape(X_s_dict)
        
                
        #------------calculate mu-----------------
        mu_s_dict={}
        mu=np.zeros([n_features,1])
        for s, X_s in X_s_dict.items():
            X_s=X_s.values
            N_s, n_features=X_s.shape
            w_s=N_s/N
            
            mu_s=np.mean(X_s,axis=0)[:, np.newaxis]
            mu_s_dict[s]=mu_s
            mu+=w_s*mu_s


        
        #------------estimate Sigma-----------------
        Gamma=np.zeros((n_features,1))
        for s, X_s in X_s_dict.items():
            X_s=X_s.values
            N_s, n_features=X_s.shape
            w_s=N_s/N
            mu_s=mu_s_dict[s]

            Gamma_s=np.zeros((n_features,1))
            for n in range(N_s):
                x_n=X_s[n,:]   
                x_n=x_n[:,np.newaxis]
                Gamma_s+=np.multiply(x_n-mu_s, x_n-mu_s)
            Gamma_s=Gamma_s/(N_s-1)
            Gamma+=w_s*Gamma_s

        Sigma=self._vector_to_diagonalmatrix(Gamma)
        
       

        #------------calculate gamma_s_dict-----------------
        gamma_s_dict={} 
        for s, X_s in X_s_dict.items():
            X_s=X_s.values
            N_s, n_features=X_s.shape
            
            gamma_s_sum=np.zeros([n_features,1])
            for n in range(N_s):
                x_n=X_s[n,:]   
                x_n=x_n[:,np.newaxis]
                gamma_s_sum+=x_n-mu     
            
            gamma_s=(1/N_s)*gamma_s_sum  #shape=(n_features,1)
            gamma_s_dict[s]=gamma_s
            

        
        #------------calculate Delta_s_dict-----------------
        inv_Sigma=self._calculate_inverse(Sigma)
        Delta_s_dict={} 
        for s, X_s in X_s_dict.items():
            X_s=X_s.values
            N_s, n_features=X_s.shape
            gamma_s=gamma_s_dict[s]
            
            Lambda_s_sum=np.zeros([n_features,1])
            for n in range(N_s):
                x_n=X_s[n,:]   
                x_n=x_n[:,np.newaxis]
                Lambda_s_sum+=np.multiply(x_n-mu-gamma_s, x_n-mu-gamma_s)
                
            Lambda_s=1/(N_s)*np.dot(inv_Sigma,Lambda_s_sum)  #shape=(n_features,1)
            Delta_s_dict[s]=self._vector_to_diagonalmatrix(Lambda_s)   #shape=(n_features,n_features)
            
       
        #-----------------------save parameters as model attributes------------------------------------
        self.gamma_s_dict, self.Delta_s_dict, self.mu, self.Sigma = gamma_s_dict, Delta_s_dict, mu, Sigma


     
    def _fit(self, data, setting_labels):
        """Estimate model paramters using the given data and setting_labels.
        
        Parameters
        -----------        
        data: feature data to be harmonized, array of shape=(n_samples,n_features).
       
        setting_labels: scanner setting labels corresponding to "data",  array of shape=(n_samples,).
             
       
        Return
        -----------
        self
       
        """
        
        
        print('Begin fitting the model...')
        
        # group data by setting labels
        data=pd.DataFrame(data)
        setting_labels=pd.DataFrame(setting_labels,columns=['setting_labels'])
        X_s_dict=dict(list(data.groupby(setting_labels['setting_labels'])))
        print('There are %d scanner settings!'%(len(X_s_dict)))
        
        for s, X_s in X_s_dict.items():
            print('scanner setting s=%d, data shape=(%d, %d)'%(s, X_s.shape[0], X_s.shape[1]))
        
        
        self._estimate_model_parameters(X_s_dict)
        print('Finish fitting the model...')
        
        return self
    
    
        
    def _harmonize_data(self, data, setting_labels):
        """Remove scanner effects using the estimated model paramters.
        
        Parameters
        -----------        
        data: feature data to be harmonized, array of shape=(n_samples,n_features).
       
        setting_labels: scanner setting labels corresponding to "data",  array of shape=(n_samples,).
             
             
        Return
        -----------
        harmonized_data: harmonized data without scanner effects corresponding to "data", array of shape=(n_samples,n_features).

        """
        
        
        print('Begin harmonizing data...')
        
        inv_Delta_s_dict={}
        for s, Delta_s in self.Delta_s_dict.items():
            inv_Delta_s_dict[s]=self._calculate_inverse(Delta_s)

        n_samples, n_features=data.shape 
        
        harmonized_data=np.zeros([n_samples,n_features])
        for n in range(n_samples):
            x_n=data[n,:]
            x_n=x_n[:, np.newaxis]  
            s=setting_labels[n]
            
            gamma_s=self.gamma_s_dict[s]
            inv_Delta_s=inv_Delta_s_dict[s]
            
            harmonized_x_n= self.mu + np.dot(np.sqrt(inv_Delta_s), x_n-self.mu-gamma_s)
            harmonized_data[n,:]=harmonized_x_n.T
         
        print('Finish harmonizing data...')
        
        return harmonized_data
    
    
    def harmonize_data(self, data, setting_labels):
        """Remove scanner effects using the estimated model paramters.
        
        Parameters
        -----------        
        data: feature data to be harmonized, array of shape=(n_samples,n_features).
       
        setting_labels: scanner setting labels corresponding to "data",  array of shape=(n_samples,).
             
        Return
        -----------
        harmonized_data: harmonized data without scanner effects corresponding to "data", array of shape=(n_samples,n_features).
       
        """
        
        self._check_input_data(data, setting_labels)
   
        self._fit(data, setting_labels)
        harmonized_data=self._harmonize_data(data, setting_labels)
        
        return harmonized_data

            

        
    
    
    

    
    