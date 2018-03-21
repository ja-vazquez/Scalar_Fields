

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import matplotlib as mpl

class Horndeski:
    def __init__(self, name='Alberto'):

        self.name = name

        self.w0 = -1.3
        self.w1 =  0.5

        self.H0     = 0.7
        self.Ocb    = 0.24
        self.Omrad  = 0.0001
        self.Ode    = 1-self.Ocb-self.Omrad

        self.lna    = np.linspace(-3, 0, 50)
        self.aa     = np.exp(self.lna)
        self.z      = np.exp(-self.lna) - 1.
        self.zvals  = np.linspace(0, 3, 20)

    def something(self):
        print self.name


    def logatoz(self, func):
        # change function from lna to z
        tmp     = interp1d(self.lna, func)
        functmp = tmp(self.lna)
        return  np.interp(self.zvals, self.z[::-1], functmp[::-1])


    def eos(self, a):
        return -1.24 + np.sin((1-a)*6/self.w1)**2#(1-a)*self.w1



    def cdet(self, a):
        return -0.24 + np.sin((1-a)*6/self.w1)**2 #(1+ self.eos(a))/(1- self.eos(a))


    def rho_de(self):
        return  (1+ self.eos(self.aa))/self.aa


    def rho_int(self, a):
        tmp     = interp1d(self.aa, self.rho_de())
        rho = np.array([quad(tmp, i, 1)[0] for i in a])
        return np.exp(3*rho)


    def lambdadet(self, a):
        return np.cos((1-a)*6/self.w1)**2 #2*(self.w1)*np.ones(len(a)) #*np.cos((1-a)*6*self.w1) #self.Ode*np.ones(len(a)) #0.5*(1- self.eos(a))*self.rho_int(a)


    def hubble(self, lna, HD=True):
        a = np.exp(lna)
        if HD:
            self.lam = self.lambdadet(a)
            self.cde = self.cdet(a)
            Ode = self.lam + self.cde #*self.lam

        else:
            Ode = self.Ode

        return self.H0*np.sqrt(self.Ocb/a**3 + self.Omrad/a**4 + Ode)


    def plot_hub(self):
        f = plt.figure(figsize=(9,10))
        ax1 = f.add_subplot(221)
        ax2 = f.add_subplot(222)
        ax3 = f.add_subplot(223)
        ax4 = f.add_subplot(224)

        #f, axes = plt.subplots(nrows=2, ncols=2, sharex=False, figsize=(9,10))
        #((ax1, ax2), (ax3, ax4)) = axes
        dataHz = np.loadtxt('Hz_all.dat')
        redshifts, obs, errors = [dataHz[:,i] for i in [0,1,2]]
        ax1.errorbar(redshifts, obs, errors, xerr=None,
                color='purple', marker='o', ls='None',
                elinewidth =2, capsize=3, capthick = 1, label='$Datos$')
        ax1.plot(self.zvals, 100*self.logatoz(self.hubble(self.lna, HD=False)),  'k-')

        min, max = (0, 1)
        step     = (max-min)/10.

        hh =[]
        cc =[]
        ll =[]
        zz =[]
        ww =[]
        P  =[]

        for i in np.arange(min, max, step):
            self.w1 =  i
            hh.append(100*self.logatoz(self.hubble(self.lna)))
            cc.append(self.logatoz(self.cde))
            ll.append(self.logatoz(self.lam))
            ww.append(self.logatoz(self.eos(self.aa)))
            zz.append(self.zvals)
            P.append(i)


        mymap    = mpl.colors.LinearSegmentedColormap.from_list('mycolors',['blue','red'])
        Z        = [[0,0]]
        levels   = np.arange(min,max+step,step)
        CS3      = plt.contourf(Z, levels, cmap=mymap)
        cbaxes = f.add_axes([0.91, 0.1, 0.03, 0.8])
        cbar   = plt.colorbar(CS3, cax = cbaxes)
        cbar.set_label('$w_1$', rotation=0, fontsize=20)

        for x,y,z in zip(zz,hh,P):
            r    = (float(z)-min)/(max-min)
            g, b = 0, 1-r
            ax1.plot(x,y, color=(r,g,b))

        for x,y,z in zip(zz,cc,P):
            r    = (float(z)-min)/(max-min)
            g, b = 0, 1-r
            ax2.plot(x,y, color=(r,g,b))


        for x,y,z in zip(zz,ll,P):
            r    = (float(z)-min)/(max-min)
            g, b = 0, 1-r
            ax3.plot(x,y, color=(r,g,b))


        for x,y,z in zip(zz,ww,P):
            r    = (float(z)-min)/(max-min)
            g, b = 0, 1-r
            ax4.plot(x,y, color=(r,g,b))

        ax3.set_xlabel('redshift $z$')
        ax4.set_xlabel('redshift $z$')
        #ax4.set_xticks( np.arange(3) )
        #ax4.set_xticklabels( np.arange(3) )


           # ax1.plot(self.zvals, 100*self.logatoz(self.hubble(self.lna)))
           # ax2.plot(self.zvals, self.logatoz(self.cde))
           # ax3.plot(self.zvals, self.logatoz(self.lam))
           # ax4.plot(self.zvals, self.logatoz(self.eos(self.aa)))

        ax1.set_ylabel('$H(z)$')
        ax2.set_ylabel('$c(a)$')
        ax3.set_ylabel('$L(a)$')
        ax4.set_ylabel('$w(a)$')
        f.subplots_adjust(hspace=0.1)
        plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
        plt.legend(loc= 'best', frameon=False)
        plt.savefig('Horn_2.pdf')
        plt.show()

if __name__ == '__main__':
    H = Horndeski()
    print H.plot_hub()