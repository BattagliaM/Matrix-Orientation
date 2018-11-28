#################################################################
# Name:     Orientation_Matrix.py                               #
# Authors:  Michael Battaglia                                   #
# Date:     8 June, 2017                                        #
# Function: Relate the Miller indices (hkl) in the reciprocal   #
#			space with its Cartesian coordinates in the         #
#			lab frame (XYZ) at zero goniometer position.        #
#################################################################

#imports
import numpy as np

##defining the orientation matrix R, where (XYZ) = R*(hkl)

#choose values for R
r11 = 1.7645050e-001
r12 = 1.1883173e-001
r13 = -4.5066652e-002
r21 = -1.1693675e-001
r22 = -1.0494427e-001
r23 = -6.1191586e-002
r31 = -1.5756760e-001
r32 = 2.1095556e-001
r33 = -5.0548963e-003 

R = np.array([[r11,r12,r13], [r21,r22,r23], [r31,r32,r33]])

##calculating lattice parameters of the reciprocal lattice

#the columns of R (vectors), ie the three axes of the reciprocal unit cell
a_star = np.array([r11,r21,r31])
b_star = np.array([r12,r22,r32])
c_star = np.array([r13,r23,r33])

#lengths of the reciprocal unit cell vectors
a_star_length = np.linalg.norm(a_star)
b_star_length = np.linalg.norm(b_star)
c_star_length = np.linalg.norm(c_star)

#the angles between the reciprocal space axes. Returns angles in radians.
alpha_star = np.arccos((np.dot(b_star,c_star))/(np.linalg.norm(b_star)*np.linalg.norm(c_star)))
beta_star = np.arccos((np.dot(a_star,c_star))/(np.linalg.norm(a_star)*np.linalg.norm(c_star)))
gamma_star = np.arccos((np.dot(a_star,b_star))/(np.linalg.norm(a_star)*np.linalg.norm(b_star)))

#volume of the reciprocal unit cell
Volume_star = a_star*b_star*c_star*np.sqrt(1-(np.cos(alpha_star))**2 -(np.cos(beta_star))**2
			 -(np.cos(gamma_star))**2 +2*np.cos(alpha_star)*np.cos(beta_star)*np.cos(gamma_star))

##obtaining the coordinates of the crystallographic axes (a,b,c) in the lab frame (XYZ)

a_0 = np.cross(b_star,c_star)/(np.linalg.norm(b_star)*np.linalg.norm(c_star)*np.sin(alpha_star))
b_0 = np.cross(c_star,a_star)/(np.linalg.norm(c_star)*np.linalg.norm(a_star)*np.sin(beta_star))
c_0 = np.cross(a_star,b_star)/(np.linalg.norm(a_star)*np.linalg.norm(b_star)*np.sin(gamma_star))

##calculating the orientation matrix R_D that relates the crystallographic axes (a,b,c)
##with the lab frame (XYZ).

phi_aX = np.arccos(np.dot(a_0,np.array([1,0,0])))/(np.linalg.norm(a_0))
phi_aY = np.arccos(np.dot(a_0,np.array([0,1,0])))/(np.linalg.norm(a_0))
phi_aZ = np.arccos(np.dot(a_0,np.array([0,0,1])))/(np.linalg.norm(a_0))
phi_bX = np.arccos(np.dot(b_0,np.array([1,0,0])))/(np.linalg.norm(b_0))
phi_bY = np.arccos(np.dot(b_0,np.array([0,1,0])))/(np.linalg.norm(b_0))
phi_bZ = np.arccos(np.dot(b_0,np.array([0,0,1])))/(np.linalg.norm(b_0))
phi_cX = np.arccos(np.dot(c_0,np.array([1,0,0])))/(np.linalg.norm(c_0))
phi_cY = np.arccos(np.dot(c_0,np.array([0,1,0])))/(np.linalg.norm(c_0))
phi_cZ = np.arccos(np.dot(c_0,np.array([0,0,1])))/(np.linalg.norm(c_0))

R_D_radians = np.array( [[phi_aX,phi_aY,phi_aZ], [phi_bX,phi_bY,phi_bZ], [phi_cX,phi_cY,phi_cZ]] )

R_D = np.degrees(R_D_radians)

print "The orientation matrix is:"
print R_D