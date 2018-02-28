
#from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
#from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.integrate import quad
from scipy.integrate import odeint
#import os
import matplotlib as mpl
import scipy.optimize as optimize

# Build a function for interpolation
# Can be chosen LCDM or SFDM

class Early:
    def __init__(self, lA=0.14, lB=2, lam=4, name='Alberto'):

        self.name   = name

        self.H0     = 0.7
        self.Ocb    = 0.23
        self.Orad   = 0.0001
        self.Ode    = 1-self.Ocb-self.Orad

        self.lna    = np.linspace(-20, 0, 300)
        self.z      = np.exp(-self.lna) - 1.
        self.zvals  = np.linspace(0, 3, 100)

        self.cte    = 3*self.H0**2

        self.lB   = lB
        self.lA   = lA
        self.lam = lam
        self.vp0 = 1.

        self.find_Ode(num=0)

    def something(self):
        print self.name




    def Vearly(self, x, select):
        #0-Potential, 1-Derivative
        lB = self.lB
        lA= self.lA/self.lam**2
        func1 = (x- lB)**2 + lA
        func2 = 2.0*(x- lB)- self.lam*func1

        if select == 0:
            return func1*np.exp(-self.lam*(x- lB))/lA*self.vp0
        if select == 1:
            return func2*np.exp(-self.lam*(x- lB))/lA*self.vp0





    def hubble(self, lna, x_vec=None, SF = False):
        a = np.exp(lna)
        if SF:
            quin, dotquin= x_vec
            Ode_quin =  0.5*dotquin**2 + self.Vearly(quin, 0)
            Ode      = Ode_quin
        else:
            Ode = self.Ode

        return np.sqrt(self.Ocb/a**3 + self.Orad/a**4 + Ode)


    def logatoz(self, func):
        # change function from lna to z
        tmp     = interp1d(self.lna, func)
        functmp = tmp(self.lna)
        return  np.interp(self.zvals, self.z[::-1], functmp[::-1])



    def mu(self, z, hub):
        tmp   = interp1d(self.zvals, 1./hub)
        xi = np.array([quad(tmp, 0.0, i)[0] for i in z])
        dL = (1+z)*xi
        return 5*np.log10(dL) + 52.5 #Checar este numero


    def RHS(self, x_vec, lna):
        quin, dotquin= x_vec
        hubble = self.hubble(lna, x_vec, SF=True)
        return [np.sqrt(3.0)*dotquin/hubble, -3*dotquin - np.sqrt(3.0)*self.Vearly(quin, 1)/hubble]


    def phidot(self, x0):
        return np.sqrt(4*self.Vearly(x0, 0))


    def solver(self, quin0, vp0):
        self.vp0 = vp0
        y0       = [quin0, self.phidot(quin0)]
        y_result = odeint(self.RHS, y0, self.lna, h0=1E-10)
        return y_result


    #Cambiar
    def find_phi0(self, x):
        return (self.Vearly(x, 1)**2/self.Vearly(x, 0)/4 - self.Vearly(x,0))*np.exp(4*self.lna[0])*(3./self.Orad) - 1.



    def find_Ode(self, num=1):
        lowphi, highphi = 0, 2
        tol, tol1       = 101, 100
        Ttol            = 5e-3
        count           = 0

        quin0 = optimize.bisect(self.find_phi0, -20, 20, xtol= 0.001) if num != 0 else 0

        while (np.abs(tol)> Ttol):
            mid = (lowphi + highphi)/2.0
            sol = self.solver(quin0, mid).T

            quin, dotq= sol

            rho_quin =   0.5*dotq[-1]**2 + self.Vearly(quin[-1], 0)
            Ode = rho_quin/self.hubble(0.0, [quin[-1], dotq[-1]], SF=True)**2
            tol = (1- self.Ocb- self.Orad) - Ode
            #print Ode, mid, self.Ode
            if(np.abs(tol) < Ttol):
                print 'reach tolerance', 'Ode=', Ode, 'phi_0=', mid, 'error=', np.abs(tol)
                break
            else:
                if(tol<0):
                    highphi = mid
                else:
                    lowphi  = mid

            count+= 1
            if (count > 10):
                mid = 0
                print 'No solution found!'
                break

        self.hub_SF   = interp1d(self.lna, self.hubble(self.lna, sol, SF=True))
        self.hub_SF_z = self.logatoz(self.hubble(self.lna, sol, SF=True))
        self.solution = sol
        return mid


    #   Plots
    #--------------------------------


    def plot_Vomegas(self, num=1):
        a    = np.exp(self.lna)
        X =[]
        Y =[]
        P =[]
        zz=[]
        hh=[]
        ww = []

        min, max = (1., 10.2)
        step     = (max-min)/20.
        mymap    = mpl.colors.LinearSegmentedColormap.from_list('mycolors',['blue','red'])
        Z        = [[0,0],[0,0]]
        levels   = np.arange(min,max+step,step)
        CS3      = plt.contourf(Z, levels, cmap=mymap)

        for i in [1]: #np.arange(min, max, step):
                self.find_Ode()
                scalars = self.solution #self.solver(phi0, phan0).T #self.solution #
                quin, dquin = scalars

                hub_SF_z = self.logatoz(self.hubble(self.lna, scalars, SF=True))
                rho_quin =   0.5*dquin**2 + self.Vearly(quin, 0)
                Omega = rho_quin/self.hubble(self.lna, scalars, SF=True)**2
                w1 = 0.5*dquin**2 - self.Vearly(quin, 0)
                w2 = 0.5*dquin**2 + self.Vearly(quin, 0)

                # Plotting what I actually want
                X.append(self.lna)
                Y.append(Omega)
                P.append(i)
                zz.append(self.zvals)
                hh.append(100*hub_SF_z)
                #ww.append(self.logatoz(w1/w2))
                ww.append(w1/w2)

        fig = plt.figure(figsize=(9,10))
        ax3 = fig.add_subplot(311)
        for x,y,z in zip(X,ww,P):
            r    = (float(z)-min)/(max-min)
            g, b = 0, 1-r
            ax3.plot(x,y, color=(r,g,b))
        plt.axhline(y=-1, color='k', linestyle='--')
        ax3.legend(loc='lower right', frameon=False)
        #plt.title('$m_\phi$ = %.1f, $m_\psi$=%.1f'%(self.mquin, self.mphan))



        ax2 = fig.add_subplot(312)
        ax2.plot(self.lna, self.Ocb/a**3/self.hubble(self.lna)**2,   'k-')
        ax2.plot(self.lna, self.Orad/a**4/self.hubble(self.lna)**2, 'k-')
        for x,y,z in zip(X,Y,P):
            r = (float(z)-min)/(max-min)
            g, b = 0, 1-r
            ax2.plot(x,y,color=(r,g,b))
        ax2.plot(self.lna, self.Ode/(self.hubble(self.lna))**2, 'o',  markersize=2)
        plt.ylabel('$w_{\phi, \psi}$')

        cbaxes = fig.add_axes([0.91, 0.1, 0.03, 0.8])
        cbar   = plt.colorbar(CS3, cax = cbaxes)
        cbar.set_label('$\psi_0$', rotation=0, fontsize=20)

        ax1 = fig.add_subplot(313)
        for x,y,z in zip(zz,hh,P):
            r = (float(z)-min)/(max-min)
            g,b = 0, 1-r
            ax1.plot(x,y,color=(r,g,b))
        ax1.plot(self.zvals, 100*self.logatoz(self.hubble(self.lna)), 'o',  markersize=2)
        dataHz = np.loadtxt('Hz_all.dat')
        redshifts, obs, errors = [dataHz[:,i] for i in [0,1,2]]
        plt.errorbar(redshifts, obs/0.7, errors, xerr=None,
                color='purple', marker='o', ls='None',
                elinewidth =2, capsize=3, capthick = 1)
        ax1.legend(loc='lower right', frameon=False)
        plt.xlabel('redshift $z$')
        #plt.savefig('Omega_phi_%.1f_psi_%.1f.pdf'%(self.mquin, self.mphan))
        plt.show()


    def plot_potential(self, select=0):
        quin = np.arange(-10, 10, 0.1)
        phan = np.arange(-10, 10, 0.1)
        quin, phan = np.meshgrid(quin, phan)
        Pot = self.Vquin(quin, select) + self.Vphan(phan, select)

        fig = plt.figure()
        ax = fig.gca(projection='3d')
        surf = ax.plot_surface(quin, phan, Pot, cmap=cm.bwr,
                       linewidth=0, antialiased=False)

        plt.ylabel('$\phi$', fontsize=20)
        plt.xlabel('$\psi$', fontsize=20)
        plt.title('$m_\phi$ = %.1f, $m_\psi$=%.1f'%(self.mquin, self.mphan))
        fig.colorbar(surf, shrink=0.5, aspect=5)
        plt.savefig('Potential_Phan_phi_%.1f_psi_%.1f.pdf'%(self.mquin, self.mphan))
        plt.show()



    def plot_hubble(self):
        hub_CDM    = self.logatoz(self.hubble(self.lna))
        plt.plot(self.zvals, self.hub_SF_z)
        plt.plot(self.zvals, hub_CDM, 'o',  markersize=2)

        dataHz = np.loadtxt('Hz_all.dat')
        redshifts, obs, errors = [dataHz[:,i] for i in [0,1,2]]
        plt.errorbar(redshifts, obs, errors, xerr=None,
                color='purple', marker='o', ls='None',
                elinewidth =2, capsize=5, capthick = 1, label='$Datos$')
        plt.xlabel(r'$z$')
        plt.ylabel(r'$H(z) [km/s Mpc^{-1}]$')
        plt.show()




    def plot_SN(self):
        names = ['name', 'z', 'mu', 'error']
        result = pd.read_table('sn_z_mu_dmu_union2.txt', sep='\s+', names=names, index_col='z')
        result=result.sort_index()
        plt.figure(figsize=(14,7))
        result['mu'].plot(yerr=result['error'], linestyle='None', label = 'SN')


        hub_CDM    = self.logatoz(self.hubble(self.lna))
        mu = self.mu(self.zvals, hub_CDM)
        plt.plot(self.zvals, mu, 'o',  markersize=2, label = '$\Omega_{DM} = 0.24, \Omega_\Lambda=0.76$')

        mu_SF = self.mu(self.zvals, self.hub_SF_z)
        plt.plot(self.zvals, mu_SF, label = 'SF', color = 'r')


        plt.xlabel(r'Corrimiento al rojo - z', fontsize = 20)
        plt.ylabel(r'Distancia modular - $\mu$', fontsize = 20)
        plt.title('Supernova Tipo Ia')
        plt.legend(loc='lower right', frameon=False)
        plt.savefig('SN_models.pdf')
        plt.show()


if __name__ == '__main__':
    Q = Early()

    #Q.find_Ode()
    #Q.plot_potential()
    Q.plot_Vomegas()
    #print Q.Vearly(1,0), Q.Vearly(1,1)
    #Q.find_Ode()
    #Q.plot_hubble()
    #Q.plot_SN()

