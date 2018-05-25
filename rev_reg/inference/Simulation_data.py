import numpy as np

#Simluated dataset. Undo when there is real data

class Simulation:
	def __init__(self, ngenes = 15000, nsamples = 600, nsnps = 1000, true_sigmabeta = 0.1):
		self._ngenes = ngenes
		self._nsamples = nsamples
		self._nsnps = nsnps
		self._true_sigmabeta = true_sigmabeta
		self._nmin = 50
		self._nmax = 250
		self._ncausal = np.random.randint(self._nmin, self._nmax, self._nsnps)
		self._Y = np.random.randn(ngenes * nsamples).reshape(ngenes, nsamples)
		self._X = np.zeros((nsnps, nsamples))

	@property
	def geno(self):
	    self._make_data()
	    return self._X

	@property
	def expr(self):
	    return self._Y

	def _make_data(self):
		for i in range(self._nsnps):
			choose = np.random.choice(self._ngenes, self._ncausal[i], replace=False)
			betas = np.random.normal(0, self._true_sigmabeta, self._ncausal[i])
			#betas = np.random.normal(0, true_sigmabeta, ngenes)
			#choose = np.arange(ngenes)
			self._X[i, :] = np.dot(self._Y[choose, :].T, betas) + np.random.normal(0, 1, self._nsamples)

		nsnps = 5000
		nsample = self._nsamples
		binom_geno = np.zeros((nsample, nsnps))
		binom_freq = np.random.uniform(0.2, 0.8, nsnps)  #0.199999, 0.8000001  #0.099999, 0.199999
		for i in range(nsnps):
			binom_geno[:, i] = np.random.binomial(2, binom_freq[i], nsample)
		#binom_freq[i] = sum(dosagelist) / 2 / nsample
		binom_geno = (binom_geno - (2 * binom_freq)) / np.sqrt(2 * binom_freq * (1 - binom_freq))    


		geno = np.vstack((self._X,binom_geno.T))
		self._X = geno.copy()
		return None