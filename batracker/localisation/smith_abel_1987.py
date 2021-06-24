# -*- coding: utf-8 -*-
"""
Smith & Abel 1987
=================
Smith & Abel 1987 develop solutions for N receivers. 


Reference
---------
Smith, J. O. and Abel, J. S., “Closed-form least-squares source location
estimation from range-difference measurements,” IEEE Trans. Acoust.
Speech 35, 1661–1669 (1987)
"""



# defined below 12
P_d = np.eye(d.shape[0]) - (d*d.T)/(d.T*d)

# eqn. 14
x_s = (1/2)*np.linalg.inv(S_T*P_d*W*P_d*S)*S_T*P_d*W*P_d*delta

