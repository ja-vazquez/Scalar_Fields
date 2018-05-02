
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
from scipy import optimize
import matplotlib as mpl



# TODO-me Can be chosen LCDM or SFDM
# TODO-me Build a function for interpolation

class Quintom:
    def __init__(self, mquin=1, mphan=1, quin0=0, phan0=0, name='Alberto'):

        self.name   = name

        self.mquin  = mquin
        self.mphan  = mphan
        self.quin0  = quin0
        self.phan0  = phan0

        self.H0     = 0.68
        self.Ocb    = 0.3
        self.Omrad  = 0.0001
        self.Odeobs = 1-self.Ocb-self.Omrad

        self.lna    = np.linspace(-15, 0, 300)
        self.z      = np.exp(-self.lna) - 1.
        self.zvals  = np.linspace(0, 3, 300)

        self.cte    = 3*self.H0**2
        self.n      = 2
        self.m      = 4


    def something(self):
        print self.name



    def Vphan(self, x, select):
        "0-Potential, 1-Derivative"
        if select == 0:     return 0.5*(x*self.mphan)**self.m
        if select == 1:     return 0.5*(x*self.mphan)**(self.m-1)*self.mphan*self.m


    def Vquin(self, x, select):
        if select == 0:     return 0.5*(x*self.mquin)**self.n
        if select == 1:     return 0.5*(x*self.mquin)**(self.n-1)*self.mquin*self.n



    def hubble(self, lna, x_vec=None, SF = False):
        a = np.exp(lna)
        if SF:
            quin, dotquin, phan, dotphan = x_vec
            Ode_quin =  0.5*dotquin**2 + self.Vquin(quin, 0)/self.cte
            Ode_phan = -0.5*dotphan**2 + self.Vphan(phan, 0)/self.cte
            Ode      = Ode_quin + Ode_phan
        else:
            Ode      = self.Odeobs

        return self.H0*np.sqrt(self.Ocb/a**3 + self.Omrad/a**4 + Ode)


    def logatoz(self, func):
        "change functions from lna to z "
        tmp     = interp1d(self.lna, func)
        functmp = tmp(self.lna)
        return  np.interp(self.zvals, self.z[::-1], functmp[::-1])



    def mu(self, z, hub):
        " Useful for plotting SN "
        tmp = interp1d(self.zvals, 1./hub)
        xi  = np.array([quad(tmp, 0.0, i)[0] for i in z])
        dL  = (1+z)*xi
        return 5*np.log10(dL) + 52.5 #Checar este numero


    def RHS(self, x_vec, lna):
        sqrt3H0  = np.sqrt(self.cte)
        quin, dotquin, phan, dotphan = x_vec
        hubble = self.hubble(lna, x_vec, SF=True)
        return [sqrt3H0*dotquin/hubble, -3*dotquin - self.Vquin(quin, 1)/(sqrt3H0*hubble),
                sqrt3H0*dotphan/hubble, -3*dotphan + self.Vphan(phan, 1)/(sqrt3H0*hubble)]



    def solver(self, quin0, dotquin0, phan_0, dotphan0):
        y0       = [quin0, dotquin0, phan_0, dotphan0]
        y_result = odeint(self.RHS, y0, self.lna, h0=1E-10)
        return y_result




    def calc_Ode(self, mid):
        phan0 = self.phan0
        if   self.mquin == 0.:   a, b = 0, mid
        elif self.mphan == 0.:   a, b = mid, 0
        else:
            a,b = np.sqrt(2*3*self.H0**2*self.Odeobs + (self.mphan*mid)**2)/self.mquin, mid

        sol = self.solver(a, 0.0, b, 0.0).T

        quin, dotq, phan, dotp = sol
        rho_quin =   0.5*dotq[-1]**2 + self.Vquin(quin[-1], 0)/self.cte
        rho_phan =  -0.5*dotp[-1]**2 + self.Vphan(phan[-1], 0)/self.cte

        Ode = (rho_quin + rho_phan)*(self.H0/self.hubble(0.0, [quin[-1], dotq[-1], phan[-1], dotp[-1]], SF=True))**2
        tol = self.Odeobs - Ode

        return tol, Ode


    def bisection(self):
        "Search for intial condition of phi such that \O_DE today is 0.7"
        lowphi, highphi = -100, 100
        Ttol            = 1E-3
        mid =(lowphi + highphi)/2.0
        while (highphi - lowphi )/2.0 > Ttol*0.01:
            ode_mid = self.calc_Ode(mid)[0]
            ode_low = self.calc_Ode(lowphi)[0]
            if(np.abs(ode_mid) < Ttol):
                print 'reach tolerance',  'phi_0=', mid, 'error=', self.calc_Ode(mid)
                return mid
            elif ode_low*ode_mid<0:
                highphi  = mid
            else:
                lowphi   = mid
            #print mid, self.calc_Ode(mid) ,(1- self.Ocb- self.Omrad)
            mid = (lowphi + highphi)/2.0

        print 'No solution found!', mid, ode_mid
        mid = 0
        return mid



    def prueba(self):
        #min, max = (0, .1)
        #step     = (max-min)/10.
        self.phan0 = 1.8
        #for i in np.arange(min, max, step):
        self.bisection()


    def search_ini(self):
        mid = self.bisection()
        if   self.mquin == 0:   a, b = 0, mid
        elif self.mphan == 0:   a, b = mid, 0
        else:
            a,b = np.sqrt(2*3*self.H0**2*self.Odeobs + (self.mphan*mid)**2)/self.mquin, mid
        sol = self.solver(a, 0.0, b, 0.0).T
        self.hub_SF   = interp1d(self.lna, self.hubble(self.lna, sol, SF=True))
        self.hub_SF_z = self.logatoz(self.hubble(self.lna, sol, SF=True))
        self.solution = sol
        return mid





    #   Plots
    #================================================================================


    def plot_Vomegas(self):
        a    = np.exp(self.lna)
        X =[]
        Y =[]
        P =[]
        zz=[]
        hh=[]
        ww = []
        Oma=[]
        vphi = []
        dphi = []

        min, max = (.1, 4)
        step     = (max-min)/10.

        mymap    = mpl.colors.LinearSegmentedColormap.from_list('mycolors',['blue','red'])
        Z        = [[0,0],[0,0]]
        levels   = np.arange(min, max+step, step)
        CS3      = plt.contourf(Z, levels, cmap=mymap)

        for i in np.arange(min, max, step):
            #self.mphan = i
            self.mquin = i

            phi0  = self.search_ini()
            if phi0 == 0:
                continue
            else:
            #if True:
                #self.solution = self.solver(phi0, 0, 0, 0).T
                scalars = self.solution
                quin, dquin, phan, dphan = scalars

                #hub_SF_z = self.logatoz(self.hubble(self.lna, scalars, SF=True))
                hub_SF_z = self.hub_SF_z
                rho_quin =   0.5*dquin**2 + self.Vquin(quin, 0)/self.cte
                rho_phan =  -0.5*dphan**2 + self.Vphan(phan, 0)/self.cte

                Omega = (rho_quin + rho_phan)*(self.H0/self.hubble(self.lna, scalars, SF=True))**2
                w1 = 0.5*dquin**2 - self.Vquin(quin, 0)/self.cte - 0.5*dphan**2  - self.Vphan(phan, 0)/self.cte
                w2 = 0.5*dquin**2 + self.Vquin(quin, 0)/self.cte - 0.5*dphan**2  + self.Vphan(phan, 0)/self.cte


                Om = self.Ocb/a**3*(self.H0/self.hubble(self.lna, self.solution, SF=True))**2

                # Plotting what I actually want
                X.append(self.lna)
                Y.append(Omega)
                P.append(i)
                zz.append(self.zvals)
                hh.append(100*hub_SF_z)
                ww.append(self.logatoz(w1/w2))
                Oma.append(Om)
                vphi.append(self.logatoz(w1))
                dphi.append(self.logatoz(w2))


        fig = plt.figure(figsize=(9,10))
        ax1 = fig.add_subplot(311)
        for x,y,z in zip(zz,ww,P):
            r    = (float(z)-min)/(max-min)
            g, b = 0, 1-r
            ax1.plot(x, y , color=(r,g,b))
        plt.axhline(y=-1, color='k', linestyle='--')
        ax1.legend(loc='lower right', frameon=False)
        plt.title('Quintom, $m_{\phi}$=1', fontsize=20)
        plt.ylabel('$w_{\phi, \psi}$', fontsize=20)



        cbaxes = fig.add_axes([0.9, 0.1, 0.03, 0.8])
        cbar   = plt.colorbar(CS3, cax = cbaxes)
        cbar.set_label('$m_\phi$', rotation=0, fontsize=20)


        ax2 = fig.add_subplot(312)
        for x,y,z in zip(zz,hh,P):
            r = (float(z)-min)/(max-min)
            g,b = 0, 1-r
            ax2.plot(x,y,color=(r,g,b))
        ax2.plot(self.zvals, 100*self.logatoz(self.hubble(self.lna)), 'o',  markersize=2)
        dataHz = np.loadtxt('Hz_all.dat')
        redshifts, obs, errors = [dataHz[:,i] for i in [0,1,2]]
        ax2.errorbar(redshifts, obs, errors, xerr=None,
                color='purple', marker='o', ls='None',
                elinewidth =2, capsize=3, capthick = 1)
        ax2.legend(loc='lower right', frameon=False)
        plt.ylabel('$H(z)$', fontsize=20)
        plt.xlabel('redshift $z$', fontsize=20)


        ax3 = fig.add_subplot(313)
        ax3.plot(self.lna, self.Ocb/a**3*(self.H0/self.hubble(self.lna))**2,   'o',  markersize=2)
        ax3.plot(self.lna, self.Omrad/a**4*(self.H0/self.hubble(self.lna))**2, 'o',  markersize=2)
        ax3.plot(self.lna, self.Odeobs*(self.H0/self.hubble(self.lna))**2, 'o',  markersize=2)
        for x,y,z in zip(X,Y,P):
            r = (float(z)-min)/(max-min)
            g, b = 0, 1-r
            ax3.plot(x,y,color=(r,g,b))

        for x,y,z in zip(X,Oma,P):
            r = (float(z)-min)/(max-min)
            g, b = 0, 1-r
            ax3.plot(x,y,color=(r,g,b))

        plt.xlabel('$\ln a$', fontsize=20)
        plt.ylabel('$\Omega(a)$', fontsize=20)
        plt.savefig('Omega_phi_%.1f_psi_%.1f.pdf'%(self.mquin, self.mphan))
        plt.show()


    def plot_potential(self, select=0):
        "Only the potentials in 3D"
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
        "Plot only Hubble data"
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
        "Plot only SN data"
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


    def contours(self):
        dir ='/Users/josevazquezgonzalez/Desktop/Desktop_Jose/work/Papers/Quintom/chains/'
        name_in = 'Quintom_phy_Planck_15+BBAO_'

        file = open('Quintom_contour.txt','w')
        ztmp = [0, 0.5, 1, 1.5, 2, 2.5]

        for j in range(4):
            i = 0
            for line in reversed(open(dir + name_in + '%s.txt'%(j+1)).readlines()):
                if i% 10==1:
                    a= line.split(' ')
                    #print a[5:7]
                    self.mquin = float(a[5]) #map(float, a[5:8])
                    self.mphan = float(a[6])
                    self.H0  = float(a[4])
                    self.Ocb = float(a[2])

                    phi0  = self.search_ini()
                    #print phi0, '****'
                    if phi0 == 0:
                        continue
                    else:
                        scalars = self.solution
                        quin, dquin, phan, dphan = scalars

                        w1 = 0.5*dquin**2 - self.Vquin(quin, 0)/self.cte - 0.5*dphan**2  - self.Vphan(phan, 0)/self.cte
                        w2 = 0.5*dquin**2 + self.Vquin(quin, 0)/self.cte - 0.5*dphan**2  + self.Vphan(phan, 0)/self.cte

                        bb = interp1d(self.zvals,  self.logatoz(w1/w2), kind='linear')

                        file.write('%s %s %s\n'%(a[0] , ' '.join(map(str, ztmp)) , ' '.join(map(str, bb(ztmp)))))

                    print j+1, i
                i+=1
                if i == 10000: break
                file.flush()
        file.close()


if __name__ == '__main__':
    Q = Quintom(mquin= 2., mphan= 2., quin0=0., phan0=1.)

    #print Q.find_Ode()
    #print Q.prueba()
    #Q.plot_potential()
    Q.plot_Vomegas()
    #Q.contours()
    #Q.plot_hubble()
    #Q.plot_SN()

