"""
Utilities for reading and writing parameters.
"""
import numpy as np


class FFDParameters():
	"""
	
	.. todo:: DOCS
	"""
	def __init__(self, n_control_points=[1, 1, 1]):
		self.conversion_unit = 1.
		
		self.a = 1.
		self.b = 1.
		self.c = 1.
		
		self.x0 = 0.
		self.y0 = 0.
		self.z0 = 0.

		self.aX = 0
		self.aY = 0
		self.aZ = 0

		self.Nx = n_control_points[0]
		self.Ny = n_control_points[1]
		self.Nz = n_control_points[2]
		
		self.muXX = np.zeros((self.Nx+1, self.Ny+1, self.Nz+1))
		self.muYY = np.zeros((self.Nx+1, self.Ny+1, self.Nz+1))
		self.muZZ = np.zeros((self.Nx+1, self.Ny+1, self.Nz+1))

		self.psi = np.diag([1./self.a, 1./self.b, 1./self.c])
		self.inv_psi = np.diag([self.a, self.b, self.c])

		self.rot_mat = np.eye(3)
		self.P0 = np.zeros(3)
		self.P1 = np.zeros(3)
		self.P2 = np.zeros(3)
		self.P3 = np.zeros(3)


	def print_info(self):
		"""
		This method prints all the FFD parameters on the screen.
		"""
		print 'conversion_unit = ' + str(self.conversion_unit) + '\n'
		print '(a, b, c) = (' + str(self.a) + ', ' + str(self.b) + ', ' + \
		str(self.c) + ')'
		print '(x0, y0, z0) = (' + str(self.x0) + ', ' + str(self.y0) + ', ' + \
		str(self.z0) + ')'
		print '(Nx, Ny, Nz) = (' + str(self.Nx) + ', ' + str(self.Ny) + ', ' + \
		str(self.Nz) + ')'
		print '\nmuXX ='
		print self.muXX
		print '\nmuYY ='
		print self.muYY
		print '\nmuZZ ='
		print self.muZZ
		print '\npsi ='
		print self.psi
		print '\ninv_psi ='
		print self.inv_psi
		print '\nrot_mat ='
		print self.rot_mat
		print '\nP0 ='
		print self.P0
		print '\nP1 ='
		print self.P1
		print '\nP2 ='
		print self.P2
		print '\nP3 ='
		print self.P3


