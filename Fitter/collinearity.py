import numpy as np

def is_pos_def(x): # check whether the input matrix is positive definite or not
    return np.all(np.linalg.eigvals(x) > 0)

def check_collinearity(jacmat):
    covmat = np.cov(jacmat.iloc[:,:-1].transpose())
    try:
        eig_val_cov, eig_vec_cov = np.linalg.eig(covmat)
        cond_no = max(eig_val_cov)/min(eig_val_cov)
        if cond_no > 1000:
            print ("\nThe model may suffer from strong collinearity as the condition no. is %.2f, greater than 1000."%(cond_no))
        else:
            print ("\nThe collinearity of the model is tolerable as the condition no. is %.2f, not greater than 1000."%(cond_no))
    except:
        print ('\nLinAlgError: failure to calculate condition number.')