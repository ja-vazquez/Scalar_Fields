
from scipy.interpolate import interp1d


class Conotours:
    def __init__(self, name='Alberto'):

        self.name = name

        self.H0     = 0.7
        self.Ocb    = 0.24
        self.Omrad  = 0.0001
        self.Ode    = 1-self.Ocb

        self.w0     = 1
        self.wa     = 1

        self.dir     = 'forJose_fgivenx/chains/'
        self.name_in = 'owaCDM_phy_BBAO+Planck_'
        self.name_out= self.dir + 'pto_1'


    def rhow(self, z):
        if z<= 3:
            x =[0.0, 1.0, 2.0, 3.0]
            y =[1.0,  self.w0, self.wa, 1.0]
            f =interp1d(x, y, kind='linear')
            #f= UnivariateSpline(x,y)
            rhow = f(z)
        else:
            rhow = 1.-self.Ocb
        return rhow


    def files(self):
        file = open(self.name_out + '.txt','w')
         #number of chains
        for j in range(2):
            i = 0
            for line in reversed(open(self.dir + self.name_in + '%s.txt'%(j+1)).readlines()):
                if i% 2==1:
                    a= line.split(' ')

                     #Parameteres in the chains
                    self.wa  = float(a[6])
                    self.w0  = float(a[5]) #map(float, a[5:8])
                    self.H0  = float(a[4])
                    self.Ocb = float(a[2])
                    prob     = float(a[0])

                    zvals = [0.0,  1.0, 2.0, 3.0]
                    rhows = map(self.rhow, zvals)

                    file.write('%s %s %s\n'%(prob ,' '.join(map(str, zvals)) , ' '.join(map(str, rhows))))

                #if i == 1000: break
                #print j+1, i
                i+=1
                file.flush()
            print 'cadena', j+1
        file.close()


if __name__ == '__main__':
    C = Conotours()
    C.files()