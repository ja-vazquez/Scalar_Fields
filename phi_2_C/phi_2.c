


// Articulo phi_2 como Materia Oscura
// Considerando potenciales efectivos para DM  

#include<stdio.h>
#include<math.h>


#define NP 5000000                                  //Numero de puntos                     
#define Ni log(1.0e-0)                              //Valor inicial para a         
#define Nf log(1.0e-6)                              //Valor final para a
#define d (Nf-Ni)/NP                                //Intervalo
#define c 3.0/2.0 


#define n 7                                          //Numero de ecuaciones
#define p 5                                           //Numero de Omegas
#define m 8                           

main()
{
  

                                                      //Definiendo variables a utilizar 
  int i,j,k;
  double t,x[n],O[n],f[n][m];

                                                      //Definiedo variables
  double rhs( double x[n], int i);
  void Runge( double x[n]);
  void ABM( int k, double x[n],double f[][m]);
  void Omegas( double O[p], double x[n]);

                                                     //Abriendo archivos para imprimir 
    FILE *fp1;
    FILE *fp2;    
    FILE *fp3;
    FILE *fp4;
    FILE *fp5;
	FILE *fp6;

      fp1=fopen("Odm_phi2.dat","w");
      fp2=fopen("Or_phi2.dat","w");
      fp3=fopen("Ode_phi2.dat","w");
      fp4=fopen("On_phi2.dat","w"); 
      fp5=fopen("Ob_phi2.dat","w");    
	  fp6=fopen("F.dat","w");                   
                                                       //Condiciones iniciales 

        f[0][0]=x[0]=sqrt(0.2295);         //x
        f[1][0]=x[1]=sqrt(0.00043);        //u
        f[2][0]=x[2]=sqrt(0.000043);        //z
        f[3][0]=x[3]=sqrt(0.73);           //l
        f[4][0]=x[4]=1.0e3;                //s      
        f[5][0]=x[5]=sqrt(0.000027);       //n
        f[6][0]=x[6]=sqrt(0.04);           //b



      Omegas(O,x); 

      fprintf(fp1,"%22.16e\t%22.16e\n",exp(Ni),O[0]);
      fprintf(fp2,"%22.16e\t%22.16e\n",exp(Ni),O[1]);
      fprintf(fp3,"%22.16e\t%22.16e\n",exp(Ni),O[2]);
      fprintf(fp4,"%22.16e\t%22.16e\n",exp(Ni),O[3]);
      fprintf(fp5,"%22.16e\t%22.16e\n",exp(Ni),O[4]);
      fprintf(fp6,"%22.16e\t%22.16e\n",exp(Ni),O[5]);
                                                      //Inicializar ABM con R.K       
    for(i=1;i<=3;i++){

      t=exp(Ni+d*i);       

                       Runge( x);  
      
      for(j=0;j<n;j++)     f[j][i]=d* rhs( x, j);

                      Omegas(O, x);

     fprintf(fp1, "%22.16e\t%22.16e\n",t , O[0]);
     fprintf(fp2, "%22.16e\t%22.16e\n",t , O[1]);
     fprintf(fp3, "%22.16e\t%22.16e\n",t , O[2]);
     fprintf(fp4, "%22.16e\t%22.16e\n",t , O[3]);
     fprintf(fp5, "%22.16e\t%22.16e\n",t , O[4]);
	 fprintf(fp6, "%22.16e\t%22.16e\n",t , O[5]);
    }
    
                                                         //Calcular valores usando ABM
     for(i=4;i<=NP;i++){

       t=exp(Ni+d*i);

                        k=i % 4+ 3;

                       ABM(k, x ,f);

                       Omegas(O, x);
 
      
      if(i%30000==0)  printf("F= %22.16e\n",O[0]+O[1]+O[2]+O[3]+O[4]);  
                                                       //Imprimir
      if(i%5000==0){ 
      fprintf(fp1, "%22.16e\t%22.16e\n",t , O[0]);
      fprintf(fp2, "%22.16e\t%22.16e\n",t , O[1]);
      fprintf(fp3, "%22.16e\t%22.16e\n",t , O[2]);
      fprintf(fp4, "%22.16e\t%22.16e\n",t , O[3]);
      fprintf(fp5, "%22.16e\t%22.16e\n",t , O[4]);
	  fprintf(fp6, "%22.16e\t%22.16e\n",t , O[5]);
    }
}
                                                      //Cerrar archivos impresos  
      fclose(fp1);
      fclose(fp2);    
      fclose(fp3);
      fclose(fp4);
      fclose(fp5);
	  fclose(fp6);
}





                                                    //Funcion__Parte derecha de las ecuaciones   
      double rhs( double x[], int i){
 
       double Pe;       
       
       Pe=2.0*x[0]*x[0]+4.0*x[2]*x[2]/3.0+4.0*x[5]*x[5]/3.0+x[6]*x[6];
	   
      if(i==0) return(-3*x[0]-x[1]*x[4]+c*Pe*x[0]);
      if(i==1) return(x[0]*x[4]+c*Pe*x[1]);  
      if(i==2) return(c*(Pe-(4.0/3.0))*x[2]);
      if(i==3) return(c*Pe*x[3]);
      if(i==4) return(c*Pe);
      if(i==5) return(c*(Pe-(4.0/3.0))*x[5]);
      if(i==6) return(c*(Pe-1.0)*x[6]);

	
	}


                                                    //Funcion__Runge-Kuta  4o Orden
      void Runge( double x[]){
  
      int j;
      double A1[n],A2[n],A3[n],A4[n],k1[n],k2[n],k3[n],k4[n];
      

      for(j=0;j<n;j++)    A1[j]= x[j]+ 0.5*(k1[j]=  d* rhs( x , j));
      for(j=0;j<n;j++)    A2[j]= x[j]+ 0.5*(k2[j]=  d* rhs( A1, j));
      for(j=0;j<n;j++)    A3[j]= x[j]+     (k3[j]=  d* rhs( A2, j));
      for(j=0;j<n;j++)                      k4[j]=  d* rhs( A3, j) ;

       for(j=0;j<n;j++)    x[j]+= (k1[j]+2.0*k2[j]+2.0*k3[j]+k4[j])/6.0;
    }




                                                   //Funcion__Metodo Numerico ABM  
    void ABM (int k, double x[],double f[][m]){
 
	int j; 
        double y[n],g[n][m];

       for(j=0;j<n;j++)   y[j]= x[j]+ (55.0*f[j][k]-59.0*f[j][k-1]+37.0*f[j][k-2]-9.0*f[j][k-3])/24.0;
         
       for(j=0;j<n;j++)   g[j][k+1]= d* rhs( y, j);
	
       for(j=0;j<n;j++)   x[j]= x[j]+ (9.0*g[j][k+1]+19.0*f[j][k]-5.0*f[j][k-1]+f[j][k-2])/24.0;

       for(j=0;j<n;j++)   f[j][k+1]=d* rhs( x, j);

                          for(j=0;j<n;j++)   f[j][k-3]=f[j][k+1];
    
     }
   

  
                                                   //Funcion_Parametros de densidad
void Omegas( double O[], double x[]){

       O[0]=x[0]*x[0]; //+x[1]*x[1];
       O[1]=x[2]*x[2];
       O[2]=x[3]*x[3];
	   O[3]=x[5]*x[5];
       O[4]=x[1]*x[1]; //x[6]*x[6];
	   O[5]=O[0]+O[1]+O[2]+O[3]+O[4];

       }

  
