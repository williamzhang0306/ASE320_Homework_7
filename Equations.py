import math
import numpy as np

class System():
    # sets up a system of equations that solves for q, the strenghts of the sources/sinks

    def __init__(self,u_inf,D):
        self.u_inf = u_inf
        self.D = D
        self.q_xloc = [-3*D, -2*D, -1*D, 0 ,D]

    def velocity(self,x,y,q):
        # returns the flow velocity at (x,y) given the strenghts of sources q

        u = self.u_inf
        for i in range(0,5):
            top = x - self.q_xloc[i]
            bottom = 2*math.pi*( (x - self.q_xloc[i])**2 + y**2)
            u += q[i]*top/bottom

        v = 0
        for i in range(0,5):
            top = y
            bottom = 2*math.pi*( (x - self.q_xloc[i])**2 + y**2)
            v += q[i]*top/bottom

        return np.array([u,v])

    def streamfunction(self,x,y,q):
        # returns the stream function value at (x,y) given the strenghts of sources q
        
        value = self.u_inf*y
        for i in range(0,5):
            x_rel = x-self.q_xloc[i]
            value += q[i]*math.atan2(y,x_rel)/ (2*math.pi)

        return value

    def equation_1(self,q):
        # velocity is parallel to surface at (-2d,1/4)
        # v(-2,0.25) dot <-1,8> = 0

        x,y = (-2*self.D, self.D/4)
        normal = np.array([-1,8])
        return np.dot(normal,self.velocity(x,y,q))

    def equation_2(self,q):
        # velocity is parallel to surface at (0,1/4)
        # v(0,0.5) dot <1,4> = 0

        x,y = (0, self.D/2)
        normal = np.array([1,4])
        return np.dot(normal,self.velocity(x,y,q))

    def equation_3(self,q):
        # constant stream function value on the surface
        # psi(-2,1/4) = 0. We choose 0 because it is a known value of psi at (2,0), which is on the surface

        x,y = (-2*self.D, self.D/4)
        return self.streamfunction(x,y,q)

    def equation_4(self,q):
        # constant stream funciton value on the surface
        # similar justification to equation_3

        x,y = (0, self.D/2)
        return self.streamfunction(x,y,q)

    def equation_5(self,q):
        # net strenghts of the sources and sinks is 0, since mass is conserved
        return sum(q)
    
    def system_constraint_error(self,q):
        # calculates the error in the system given a set of strengths, q.
        # In the way we arranged the equations (constraints), equation_1...equation_5 = 0

        # all constraints in system are satisifed if all equations = 0
        # error measures how far away a given q is from satisfying all constraints

        error = 0
        expected_values = [0,0,0,0,0]

        for i, equation in enumerate([self.equation_1,self.equation_2,self.equation_3,self.equation_4,self.equation_5]):
            error += (expected_values[i]-equation(q))**2

        return error
