# # Developing a Coordinate Transformation Routine for Solving Hamilton's EOM in Spherical Coordinates
# 
# ### Author: Siddharth Mahesh
# 
# ## Step : Import Modules

import numpy as np
import sympy as sp


# 
# ## Step : Transforming Cartesian Coordinates to Spherical Polar

def TransformCartesianToSpherical(values,off,sep):
    x , y , z = values[0] , values[1] , values[2]
    px , py , pz = values[off+sep*0] , values[off+sep*1] , values[off+sep*2]
    r = np.sqrt(x*x + y*y + z*z)
    phi = np.arctan2(y,x)
    theta = np.arccos(z,r)
    c_phi , s_phi = np.cos(phi) , np.sin(phi)
    c_theta , s_theta = np.cos(theta) , np.sin(theta)
    sph_values = np.zeros(len(values))
    sph_values[0] = r
    sph_values[1] = phi
    sph_values[2] = theta
    sph_values[off+sep*0] = px*c_phi*s_theta + py*s_phi*s_theta + pz*c_theta
    sph_values[off+sep*1] = -px*r*s_phi*s_theta + py*c_phi*s_theta
    sph_values[off+sep*2] = px*r*c_phi*c_theta + py*r*s_phi*c_theta - pz*r*s_theta
    return sph_values   


# ## Step : Transforming Spherical Polar Coordinates to Cartesian


def TransformSphericalToCartesian(values,off,sep):
    r , phi , theta = values[0] , values[1] , values[2]
    pr , pphi , ptheta = values[off+sep*0] , values[off+sep*1] , values[off+sep*2]
    c_phi , s_phi = np.cos(phi) , np.sin(phi)
    c_theta , s_theta = np.cos(theta) , np.sin(theta)
    cart_values = np.zeros(len(values))
    cart_values[0] = r*c_phi*s_theta
    cart_values[1] = r*s_phi*s_theta
    cart_values[2] = r*c_theta
    cart_values[off+sep*0] = pr*c_phi*s_theta + (-pphi*s_phi*s_theta + ptheta*c_phi*c_theta)/r
    cart_values[off+sep*1] = pr*s_phi*s_theta + (pphi*c_phi*s_theta + ptheta*s_phi*c_theta)/r
    cart_values[off+sep*2] = pr*c_theta - ptheta*s_theta/r
    return sph_values


# ## Step : The Jacobian Matrix of Transformation

def JacobianSphericalbyCart(values):
    r , phi , theta = values[0] , values[1] , values[2]
    pr , pphi , ptheta = values[off+sep*0] , values[off+sep*1] , values[off+sep*2]
    c_phi , s_phi = np.cos(phi) , np.sin(phi)
    c_theta , s_theta = np.cos(theta) , np.sin(theta)
    rinv = 1/r
    csc_theta = 1/s_theta
    jacobiansphericalbycart = np.array([[c_phi*s_theta, s_theta*s_phi, c_theta,0, 0, 0], [(csc_theta*s_phi)*rinv, -((c_phi*csc_theta)*rinv), 0, 0, 0,0], [-((c_theta*c_phi)*rinv), -((c_theta*s_phi)*rinv), s_theta*rinv, 0, 0, 0], [(ptheta*c_theta*c_phi - pphi*s_theta*s_phi)*(rinv*rinv), (pphi*c_phi*s_theta + ptheta*c_theta*s_phi)*(rinv*rinv), -((ptheta*s_theta)*(rinv*rinv)), c_phi*s_theta, s_theta*s_phi, c_theta], [(pphi*c_phi*s_theta + ptheta*c_theta*s_phi + pr*r*s_theta*s_phi)*rinv, (-c_phi*(ptheta*c_theta + pr*r*s_theta) + pphi*s_theta*s_phi)*rinv,0, -r*s_theta*s_phi, r*c_phi*s_theta, 0], [(ptheta*c_phi*s_theta - c_theta*(pr*r*c_phi + pphi*s_phi))*rinv, (pphi*c_theta*c_phi - pr*r*c_theta*s_phi + ptheta*s_theta*s_phi)*rinv, (ptheta*c_theta)*rinv + pr*s_theta, r*c_theta*c_phi, r*c_theta*s_phi, -r*s_theta]])
    return jacobiansphericalbycart


# ## Step : Construct the Ibar matrix for getting the EOM Rhs

def ConstructIbar():
    return np.block([[np.zeros([3,3]) , np.identity(3)],[-np.identity(3) , np.zeros([3,3])]])


# ## Step : The EOM for Spherical Polar Coordinates in Symplectic Form


def getEOMsphericalRHS(coords_spherical,hamiltonian_deriv_funcs):
    coords_cart = TransformSphericaltoCartesian(coords_spherical,0,1)
    J_matrix = JacobianSphericalbyCart(coords_spherical)
    Ibar_matrix = ConstructIbar()
    cart_derivs = np.zeros(len(coords_spherical))
    for i in range(len(coords_spherical)):
        cart_derivs[i] = hamiltonian_deriv_funcs[i](coords_cart)
    return np.dot(J_matrix,np.dot(Ibar_matrix,cart_derivs))

