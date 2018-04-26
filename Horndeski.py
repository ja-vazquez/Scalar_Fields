

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl


class Horndeski:
    def __init__(self, name='Alberto'):

        self.name = name

        self.H0     = 0.7
        self.Ocb    = 0.24
        self.Omrad  = 0.0001
        self.Ode    = 1-self.Ocb

        self.z_fin  = 12.
        self.z_cut  = 3.

        self.w0     = 0.7
        self.wa     = 0.71
        self.wb     = 0.67
        self.wc     = 0.7

        self.lna    = np.linspace(-3, 0, 50)
        self.aa     = np.exp(self.lna)
        self.z      = np.exp(-self.lna) - 1.
        self.zvals  = np.linspace(0, self.z_fin, 25)

        #Initialize, compute integral
        self.inter()


    def something(self):
        print self.name



    def rhow(self, z):
        if z<= self.z_cut:
            x =[0.0, 0.8, 1.6, 2.4, self.z_cut]
            y =[self.w0,  1-self.Ocb, 1-self.Ocb, 1.-self.Ocb, 1.-self.Ocb]
            #self.wa, self.wb, self.wc, 1] #.-self.Ocb]
            f =interp1d(x, y, kind='linear')
            #f= UnivariateSpline(x,y)
            rhow = f(z)
        else:
            rhow = 1.-self.Ocb
        return rhow


    def func_L(self, z):
        return self.rhow(z)*(1.+z)**(-7.)


    def Integral_L(self, z):
        integral = quad(self.func_L, 0, z)[0]
        calc     = -6*(1.+z)**6*(integral - self.test)
        return calc


    def inter(self):
        self.test = quad(self.func_L, 0, self.z_fin)[0]
        self.L_int  =[self.Integral_L(i) for i in self.zvals]
        self.interp      = interp1d(self.zvals, self.L_int, kind='linear')
        return 1


    def hubble(self, lna, HD=True):
        a = np.exp(lna)
        z = 1./a - 1.0


        if (z>=0  and z< self.z_cut):
                Ode = self.interp(z)
        else:
                Ode = 1.0 - self.Ocb #self.test #self.interp(self.z_fin)

        return self.H0*np.sqrt(self.Ocb/a**3 + self.Omrad/a**4 + Ode)


    def logatoz(self, func):
        "change functions from lna to z "
        tmp     = interp1d(self.lna, func)
        functmp = tmp(self.lna)
        return  np.interp(self.zvals, self.z[::-1], functmp[::-1])


 #   Plots
    #================================================================================


    def plot_c(self):
        plt.plot(self.zvals, [self.rhow(i) for i in self.zvals], 'r-')
        int_co = [self.Integral_L(i)  for i in self.zvals]
        plt.plot(self.zvals, int_co )
        #plt.plot(self.zvals, int_co-self.rhow(self.zvals), 'g-' )
        plt.xlim(0, 4)
        plt.ylim(0.6, 0.8)
        plt.show()


    def contours(self):
        #dir = '/Users/josevazquezgonzalez/Desktop/SimpleMC/SimpleMC/trunk/chains/'
	dir ='/Users/josevazquezgonzalez/Desktop/Desktop_Jose/work/Papers/Horndeski/Horn_chains/'
	which = '_1p_c'
	name_in = 'Horn'+which+'_phy_Planck+BBAO_'

        file = open('Lambda'+which+'.txt','w')
        file2 = open('rho'+which+'.txt','w')
        file3 = open('c'+which+'.txt','w')

        zv = [0.0, 0.8, 1.6, 2.4, self.z_cut]

        for j in range(2):
            i = 0
            for line in reversed(open(dir + name_in + '%s.txt'%(j+1)).readlines()):
                if i% 10==1:
                    a= line.split(' ')
                    #print a[5:7]
                    self.w0 = float(a[5]) #map(float, a[5:8])
                    self.H0  = float(a[4])
                    self.Ocb = float(a[2])
                    self.inter()

                    ls = [self.w0,  1-self.Ocb, 1-self.Ocb, 1- self.Ocb, 1- self.Ocb] #]self.wa, self.wb, self.wc, 1] #-self.Ocb]
                    rs = self.interp(zv)
                    cs = rs - ls

                    c= ' '.join(map(str, zv))

                    file.write('%s %s %s\n'%(a[0] ,c , ' '.join(map(str, ls))))
                    file2.write('%s %s %s\n'%(a[0] , c , ' '.join(map(str, rs))))
                    file3.write('%s %s %s\n'%(a[0] , c , ' '.join(map(str, cs))))
                    print j+1, i
                i+=1
                #if i == 1000: break
                file.flush()
                file2.flush()
                file3.flush()
        file.close()
        file2.close()
        file3.close()

        #with open(dir+'Horn_phy_Planck_15+SN+BBAO_2.txt') as f:
        #    print  [line.rstrip('\n') for line in f][0]

       # file = open('testfile.txt','w')

       # for i in range(5000):
       #     a , self.w0,  self.wa, self.wb, self.wc = df.iloc[i].values
       #     self.inter()
       #     b= ' '.join(map(str, self.interp([0.0,    1.0,        2.0,    self.z_fin])))
       #     file.write('%e %s %f %f %f %f\n'%(a , b, self.w0,  self.wa, self.wb, self.wc))
       # file.close()

#        self.w0     = 0.1
#        self.wa     = -0.1
#        self.wb     = -0.1
#        self.wc     = 0.1
#        self.inter()
        #print self.interp([0.0,    1.0,        2.0,    self.z_fin])

        #df.insert(1,0,0)
        #df.insert(2,1.0,1.0)
        #df.insert(3,2.0,2.0)
        #df.insert(4,3.0,self.z_fin)
        #df.to_csv('horn.csv', sep=' ', header=False, index=False)
        #df.insert([1,2,3,4],[1,2,3,4],[0,1,2,self.z_fin])
        #print df

if __name__ == '__main__':
    H = Horndeski()
    #H.plot_c()
    H.contours()
