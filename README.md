# Scalar_Fields
Scalar fields in the universe

[**Aqui los calculos**](https://nbviewer.jupyter.org/github/ja-vazquez/Scalar_Fields/blob/master/Cosmo_rhos.ipynb)

**Articulos de referencia**

[Reconstruccion w(z)](d)
 
[Quintom](https://arxiv.org/abs/astro-ph/0410654)

La evolucion de la ecuacion de estado, se parecen?

La idea es basica:
		
		Un campo escalar, con un potencial V(Phi)
		puede ser el principal responsable de la
		Energia oscura?
		
En este repo agregue la cosmologia del fondo
con un campo escalar $V(\phi)=V_0 e^{-\lam \phi}$
que sustituye a la constante cosmologica.
		
Teoria: 
		
		En la literaruta existen muchos candidatos a energia 
		oscura. Dos candidatos en particular: un fluido perfecto 
		con ecuacion de estado w(z) y un campo escalar con potencial V(phi). 
		Existe algun potential preferido, puede un solo
		campo o son necesarios dos? -- Phantom
		
<img src="https://github.com/ja-vazquez/Scalar_Fields/blob/master/Omegas_LCDM.jpg" widt="100p" height="200"/>
<img src="https://github.com/ja-vazquez/Scalar_Fields/blob/master/Omegas_V.jpg" widt="200p" height="200"/>
		
Reconstruccion 

		Se eligio un conjunto de parametros para 
		\lam, V_0, y con base a estos se construyo un conjunto de datos
		para realizar un forecast, simulando experimentos idealizados en el 
		futuro.
		Es posible recuperar la forma de este potencial? :
		Celia E. va a hacer un blind test.
		
Comparacion con las observaciones.

		Para realizar la estimacion de parametros necesitamos codigos
		eficientes, y por supuesto, datos cosmologicos actualizados.
		Para esto recurriremos en primera instancia al codigo simpleMC.
		Ademas, si las queremos considerar perturbaciones tendremos que considerar
		MontePython o CosmoMC (tambien para calcular evidencia).
		 



<img src="https://github.com/ja-vazquez/Scalar_Fields/blob/master/rhos_V.jpg" widt="200p" height="200"/>

<img src="https://github.com/ja-vazquez/Scalar_Fields/blob/master/Hz_V.jpg" widt="200p" height="200"/>
