import numpy as np
from scipy.optimize import minimize

class OptimizeRegularizer:

    # Tolerance can be changed to check for parameter convergence.
    def __init__(self, genotype, gene_expr, sigmax = 1, tol = 0.01, sigmabeta = -2):
        self._genotype = genotype
        self._expr = gene_expr
        self._sigmax = sigmax
        self._sigmareg = sigmabeta
        self._tolerance = tol
        self._niter = 0
        self._U, self._S, self._vt = self._svd(self._expr)
        
    #Decompose gene expression matrix only once   
    def _svd(self, expr):
        U, S, Vt = np.linalg.svd(np.transpose(expr),full_matrices=False)
        return (U, S, Vt)

    @property
    def sigmareg(self):
        return self._sigmareg


    @property
    def niter(self):
        return self._niter
    
    # used in log marginal likelihood calculation
    def _determinant_lemma(self, Y, sigmax, sigmabeta):
        ngenes = Y.shape[0]
        sigmabeta2 = sigmabeta * sigmabeta
        sigmax2 = sigmax * sigmax
        Yt = Y.T # shape N x G

        A = np.dot(Yt, Yt.T) * sigmabeta2 / sigmax2
        A[np.diag_indices(A.shape[0])] += 1
        logdetA = np.linalg.slogdet(A)
        logdet = 2 * ngenes * (np.log(sigmax) - np.log(sigmabeta)) + logdetA[1]
        return logdet

    # Log marginal likelihood calculation
    def _logml(self, logsigmab, *args):
        sigmabeta = np.e ** logsigmab
        # X, Y, S, U, sigmax):
        X, Y, S, U, sigmax = args
        nsnps      = X.shape[0]
        nsamples   = X.shape[1]
        ngenes     = Y.shape[0]
        sigmabeta2 = sigmabeta * sigmabeta
        sigmax2    = sigmax    * sigmax

        logdetA = self._determinant_lemma(Y , sigmax, sigmabeta) #A = np.dot(Yt.T, Yt)    
        
        Smod = np.diag(np.square(S) / (np.square(S) + sigmax2 / sigmabeta2))
        
        W = np.dot(U, np.dot(Smod, U.T))

        const_term = - 0.5 * ngenes * np.log(2 * np.pi * sigmabeta2) - 0.5 * logdetA

        snp_term = 0.5 * np.diag(np.dot(X, np.dot(W, X.T))) / sigmax2

        lml = nsnps * (const_term) + np.sum(snp_term)

        return -lml 

    # This is where optimization occurs.
    def update (self):
        sigmareg_old = np.e ** self._sigmareg
        N = 200 # expected number of true trans-eQTLs. Hard-coded
        iterate = True
        
        while iterate:
            fact = self._sigmax ** 2 / sigmareg_old ** 2
            diagw = np.square(self._S) / (np.square(self._S) + fact)
            Wsvd = np.dot(self._U, np.dot(np.diag(diagw), self._U.T))

            Qsvd  = np.diag(np.dot(self._genotype, np.dot(Wsvd,  self._genotype.T)))

            top_Q_indx = np.argsort(-Qsvd)[0:N]

            top_geno = self._genotype[top_Q_indx,:]

            #optimize using Nelder-mead simplex based simple optimization for single parameter.
            arguments = (top_geno, self._expr, self._S, self._U, self._sigmax)
            x = minimize(self._logml, [0], method="Nelder-Mead", args = arguments)#, jac=gradient_lmll)#, options=opts) #Initial guess is logsigmabeta = -2
            sigbeta_new = np.e**x.x[0]
            
            checksigma = self.check_convergence(sigbeta_new, sigmareg_old)

            if checksigma:
                iterate = False
                sigma_optimized = sigbeta_new
                
            self._niter += 1
            sigmareg_old = sigbeta_new 
            
        self._sigmareg = sigma_optimized
        
    # Check if the difference between sigmabeta old and new is just 1% which means it has conversed.
    def check_convergence(self, x, xold):
        check = False
        tol = self._tolerance
        diff = x - xold
        if diff == 0:
            check = True
        if not check and xold != 0.0:
            diff_percent = abs(diff) / xold
            if diff_percent < tol:
                check = True
            
        return check

