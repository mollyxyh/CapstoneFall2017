import numpy as np

def is_pos_def(x): # check whether the input matrix is positive definite or not
    '''@var x:input matrix'''
    return np.all(np.linalg.eigvals(x) > 0)

def check_collinearity(jacmat):
    '''@var jacmat:the Jacobian matrix of linear optimization we are using to check collinearity'''
    jacmat = jacmat.iloc[:,:-1]

    # print the jacobian matrix from the optimization
    jacmat = np.matrix(jacmat)
#    print "The Jacobian matrix is:\n", jacmat

    # construct the Hessian Matrix from the Jacobian matrix
    hess = np.dot(jacmat, jacmat.T)
#    print "\nThe Hessian Matrix is:\n", hess

#    if is_pos_def(hess):
#        print "\nThe Hessian Matrix is postive definite, thus the parameters can minimize the loss function"
#    else:
#        print "\nThe Hessian Matrix is not postive definite, thus the parameters may not minimize the loss function"

    # The condition no. is defined as: Max(eigenvalues)/Min(eigenvalues)
    cond_no = np.linalg.cond(hess)
    print "\nThe condition no. of the Hessian Matrix is:", cond_no

    if cond_no > 5000:
        print "\nThe model may suffer from strong collinearity as the condition no. is greater than 5000."
    else:
        print "\nThe collinearity of the model is tolerable as the condition no. is not greater than 5000."