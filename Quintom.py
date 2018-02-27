
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

# Build a function for interpolation
# Can be chosen LCDM or SFDM

class Quintom:
    def __init__(self, mquin=1, mphan=1, phan0=0, name='Alberto'):

        self.name   = name

        self.mquin  = mquin
        self.mphan  = mphan
        self.phan0  = phan0

        self.H0     = 0.7
        self.Ocb    = 0.24
        self.Omrad  = 0.0001
        self.Ode    = 1-self.Ocb-self.Omrad

        self.lna    = np.linspace(-10, 0, 300)
        self.z      = np.exp(-self.lna) - 1.
        self.zvals  = np.linspace(0, 3, 100)

        self.cte    = 3*self.H0**2

        if name == 'Early':
            self.B   = 1.
            self.A   = 1.
            self.lam = 1.

    def something(self):
        print self.name



    def Vphan(self, x, select):
        #0-Potential, 1-Derivative
        if select == 0:     return 0.5*(x*self.mphan)**2
        if select == 1:     return x*self.mphan**2


    def Vquin(self, x, select):
        #0-Potential, 1-Derivative
        if self.name == 'Early':
            diff = x-self.B
            V0   = 0.5*self.mquin**2
            V = V0*(diff**2+ self.A)*np.exp(-self.lam*x)
            if select == 0:
                    return V
            if select == 1:
                    return 2*V0*diff*np.exp(-self.lam*x) -self.lam*V
        else:
            if select == 0:
                    return 0.5*(x*self.mquin)**2
            if select == 1:
                    return x*self.mquin**2



    def hubble(self, lna, x_vec=None, SF = False):
        a = np.exp(lna)
        if SF:
            quin, dotquin, phan, dotphan = x_vec
            Ode_quin =  0.5*dotquin**2 + self.Vquin(quin, 0)/self.cte
            Ode_phan = -0.5*dotphan**2 + self.Vphan(phan, 0)/self.cte
            Ode      = Ode_quin + Ode_phan
        else:
            Ode = self.Ode

        return self.H0*np.sqrt(self.Ocb/a**3 + self.Omrad/a**4 + Ode)


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

        quin, dotquin, phan, dotphan = x_vec
        hubble = self.hubble(lna, x_vec, SF=True)
        return [np.sqrt(3.0)*self.H0*dotquin/hubble, -3*dotquin - self.Vquin(quin, 1)/(np.sqrt(3.0)*self.H0*hubble),
                np.sqrt(3.0)*self.H0*dotphan/hubble, -3*dotphan + self.Vphan(phan, 1)/(np.sqrt(3.0)*self.H0*hubble)]



    def solver(self, quin0, dotquin0, phan_0, dotphan0):
        y0       = [quin0, dotquin0, phan_0, dotphan0]
        y_result = odeint(self.RHS, y0, self.lna, h0=1E-10)
        return y_result


    def find_Ode(self, phan0=0.3):
        lowphi, highphi = 0, 5
        tol, tol1       = 101, 100
        Ttol            = 1e-3
        count           = 0

        while (np.abs(tol)> Ttol):
            mid = (lowphi + highphi)/2.0
            if self.name == 'Early':
                dotphi0 = np.sqrt(4*self.Vquin(mid, 0))
            else:
                dotphi0 = 0
            sol = self.solver(mid, dotphi0, phan0, 0).T
            quin, dotq, phan, dotp = sol

            rho_quin =   0.5*dotq[-1]**2 + self.Vquin(quin[-1], 0)/self.cte
            rho_phan =  -0.5*dotp[-1]**2 + self.Vphan(phan[-1], 0)/self.cte

            Ode = (rho_quin + rho_phan)*self.H0**2/self.hubble(0.0, [quin[-1], dotq[-1], phan[-1], dotp[-1]], SF=True)**2
            tol = (1- self.Ocb- self.Omrad) - Ode
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


    def plot_Vomegas(self):
        a    = np.exp(self.lna)
        #cte  = 3*self.H0**2
        X =[]
        Y =[]
        P =[]
        zz=[]
        hh=[]
        ww = []

        min, max = (1.1, 1.2)
        step     = (max-min)/20.
        mymap    = mpl.colors.LinearSegmentedColormap.from_list('mycolors',['blue','red'])
        Z        = [[0,0],[0,0]]
        levels   = np.arange(min,max+step,step)
        CS3      = plt.contourf(Z, levels, cmap=mymap)

        for i in np.arange(min, max, step):
            phan0 = i #self.phan0
            phi0  = self.find_Ode(phan0)
            print '***', phi0, phan0
            if phi0 == 0:
                continue
            else:
                scalars = self.solution #self.solver(phi0, 0, phan0, 0).T
                quin, dquin, phan, dphan = scalars

                hub_SF_z = self.logatoz(self.hubble(self.lna, scalars, SF=True))
                rho_quin =   0.5*dquin**2 + self.Vquin(quin, 0)/self.cte
                rho_phan =  -0.5*dphan**2 + self.Vphan(phan, 0)/self.cte

                Omega = (rho_quin + rho_phan)*self.H0**2/self.hubble(self.lna, scalars, SF=True)**2
                w1 = 0.5*dquin**2 - self.Vquin(quin, 0)/self.cte - 0.5*dphan**2  - self.Vphan(phan, 0)/self.cte
                w2 = 0.5*dquin**2 + self.Vquin(quin, 0)/self.cte - 0.5*dphan**2  + self.Vphan(phan, 0)/self.cte

                # Plotting what I actually want
                X.append(self.lna)
                Y.append(Omega)
                P.append(i)
                zz.append(self.zvals)
                hh.append(100*hub_SF_z)
                ww.append(self.logatoz(w1/w2))

        fig = plt.figure(figsize=(9,10))
        ax3 = fig.add_subplot(311)
        for x,y,z in zip(zz,ww,P):
            r    = (float(z)-min)/(max-min)
            g, b = 0, 1-r
            ax3.plot(x,y, color=(r,g,b))
        plt.axhline(y=-1, color='k', linestyle='--')
        ax3.legend(loc='lower right', frameon=False)
        plt.title('$m_\phi$ = %.1f, $m_\psi$=%.1f'%(self.mquin, self.mphan))



        ax2 = fig.add_subplot(312)
        ax2.plot(self.lna, self.H0**2*self.Ocb/a**3/self.hubble(self.lna)**2,   'k-')
        ax2.plot(self.lna, self.H0**2*self.Omrad/a**4/self.hubble(self.lna)**2, 'k-')
        for x,y,z in zip(X,Y,P):
            r = (float(z)-min)/(max-min)
            g, b = 0, 1-r
            ax2.plot(x,y,color=(r,g,b))
        ax2.plot(self.lna, self.H0**2*self.Ode/(self.hubble(self.lna))**2, 'o',  markersize=2)
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
        plt.errorbar(redshifts, obs, errors, xerr=None,
                color='purple', marker='o', ls='None',
                elinewidth =2, capsize=3, capthick = 1)
        ax1.legend(loc='lower right', frameon=False)
        plt.xlabel('redshift $z$')
        plt.savefig('Omega_phi_%.1f_psi_%.1f.pdf'%(self.mquin, self.mphan))
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
    Q = Quintom(mquin=4.0, mphan=0., phan0=0., name='Early')

    #Q.plot_potential()
    Q.plot_Vomegas()
    #Q.find_Ode()
    #Q.plot_hubble()
    #Q.plot_SN()

