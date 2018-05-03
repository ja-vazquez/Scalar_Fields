import numpy as np
from scipy.integrate import odeint
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


class DarkM:
    def __init__(self, mquin=1, name='Alberto'):

        self.name   = name
        self.cte    = 3/2.0


        self.NP = 50000 #000
        self.Ni = np.log(1.0e-0)
        self.Nf = np.log(1.0e-6)
        self.d  = (self.Nf- self.Ni)/self.NP
        self.t = [np.exp(self.Ni+self.d*i) for i in np.arange(self.NP)]


    def something(self):
        print self.name


    def ode_int_rk(self, func, y_0, t, args={}):
        # Inicia el arreglo de las aproximaciones
        y = np.zeros([len(t), len(y_0)])
        y[0] = y_0
        for i, t_i in enumerate(t[:-1]):

            h = self.d #t[i+1] - t_i
            k_1 = func(t_i, y[i], args)
            k_2 = func(t_i+h/2., y[i]+h/2.*k_1, args)
            k_3 = func(t_i+h/2., y[i]+h/2.*k_2, args)
            k_4 = func(t_i+h, y[i]+h*k_3, args)

            y[i+1] = y[i] + h/6.*(k_1 + 2.*k_2 + 2.*k_3 + k_4) # RK4 step

        return y


    def solver(self):
        y0       = np.array([np.sqrt(0.2295), np.sqrt(0.00043), np.sqrt(0.000043), np.sqrt(0.73), 1.0E3])
        y_result = self.ode_int_rk(self.RHS, y0, self.t)
        return y_result


    def RHS(self, t, y, args):
        x0, x1, x2, x3, x4 = y
        Pe = 2.0*x0*x0+4.0*x2*x2/3.0
        return np.array([-3*x0-x1*x4+self.cte*Pe*x0, x0*x4+self.cte*Pe*x1, self.cte*(Pe-(4.0/3.0))*x2,
                self.cte*Pe*x3, self.cte*Pe])


    def plot(self):
        x, u, z, l, s = self.solver().T
        fig = plt.figure(figsize=(9,10))
        ax3 = fig.add_subplot(111)

        i =0
        tiempo = []
        phi = []
        dphi = []
        for t, p, d  in zip(self.t, x, u):
            if i%100 ==0:
               tiempo.append(t)
               phi.append(p)
               dphi.append(d)
            i+=1
        ax3.semilogx(tiempo, np.array(phi)*np.array(phi), 'blue')
        ax3.semilogx(tiempo, np.array(dphi)*np.array(dphi), 'red')
        plt.xlabel('$\ln a$', fontsize=20)
        plt.ylabel('$\Omega(a)$', fontsize=20)
        plt.savefig('Omega_dark_matter.pdf')
        plt.show()







if __name__ == '__main__':
    DM = DarkM(mquin= 2.)
    print DM.plot()


