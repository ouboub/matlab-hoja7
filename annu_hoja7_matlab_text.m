%% Practicas de Matlab
%% Metodos adaptativos
%% Hoja 7
% *Nombre:*
% 
% *Apellido:*
% 
% *DNI:*
% 
% *Email:*
%% 
% %% *Par encajado*
% $h_{op}$ está dado por 
% 
% $$      h_{opt}=\mbox{mín}      \left(        \mbox{{ HMAX}},h_{n}\mbox{mín}\left(\mbox{{FACMAX,FAC}}          
% \left(            \frac{TOL    h_{n}}{\mathrm{Error}}          \right)^{\frac{1}{p+1}}        
% \right)      \right)$$
% 
% _FACMAX_ se suele tomar entre $1.5$ y $5$. 
% 
% %% 
%% Práctica 1: (Runge-Kutta Fehlberg-4-5) 
% La tabla de RKF-45 está dada por
% 
% $$    \begin{array}{c|cccccc}      0 &  & & & & & \\[0.2cm]      \frac{1}{4}& 
% \frac{1}{4} & & & & & \\[0.2cm]      \frac{3}{8}& \frac{3}{32}& \frac{9}{32} 
% &  & & &\\[0.2cm]      \frac{12}{13} & \frac{1932}{2197} & -\frac{7200}{2197} 
% &      \frac{7296}{2197} & & & \\[0.2cm]      1 & \frac{439}{216} & -8 & \frac{3680}{513} 
% & -\frac{845}{4104} & &\\[0.2cm]      \frac{1}{2}& -\frac{8}{27}& 2 & -\frac{3544}{2565} 
% & \frac{1859}{4104}      & -\frac{11}{40} & \\[0.2cm]      \hline      & \frac{25}{216} 
% & 0 & \frac{1408}{2565} & \frac{2197}{4104} &      -\frac{1}{5}      & 0\\[0.2cm]      
% & \frac{16}{135} & 0 & \frac{6656}{12825} & \frac{28561}{56430} & -\frac{9}{50}        
% & \frac{2}{55}\\[0.2cm]      \hline      & \frac{1}{360}& 0 & -\frac{128}{4275} 
% & - \frac{2197}{75240}&      \frac{1}{50} & \frac{2}{55}    \end{array}$$
% 
% con el error:
% 
% $$    ERR = h \left| \frac{1}{360} k_1  - \frac{128}{4275} k_3 -      \frac{2197}{75240} 
% k_4 + \frac{1}{50} k_5 + \frac{2}{55} k_6\right| . $$
% 
% Implementa dicho método llamando la funcion *mirk45fehlberg* con la sintaxis 
%%
% 
%  function  [t,y,ev,hchng_vec,err_vec]=mirk45fehlberg(f,intv,y0,TOL,hmin,hmax)
%
%% 
% * hmin= $10^{-5}$
% * $hmax=\displaystyle\frac{(intv(2)-intv(1))}{50}$
% * TOL=0.01;
%% 
%% Práctica 2 (Euler mejorado-Euler (2-1) (método de extrapolación) 
% Consideramos el siguiente método de extrapolación local con el tablero:
% 
% $$    \begin{array}{l|l}      \begin{array}{l}        {0}	\\        {1}      
% \end{array}      &	       \begin{array}{ll}        0 & 0 \\        1 & 0      
% \end{array}      \\      \hline      y_{n+1}&      \begin{array}{cc}        
% \frac{1}{2}& \frac{1}{2}	      \end{array}      \\      \hline      \hat{y}_{n+1}&      
% \begin{array}{cc}        1& 0      \end{array}    \end{array}$$
% 
% que
%% 
% * avanza con el método de Euler mejorado y
% * estima con el método de Euler.
% * en este caso el método de avance es de orden 2, pero a cambio hay que hacer 
% dos evaluaciones de función por paso.
%% 
% Implementa dicho método llamando la funcion *mieuler21* con la sintaxis 
%%
% 
%  function  [t,y,ev,hchng_vec,err_vec]=mieuler21(f,intv,y0,TOL,hmin,hmax)
%
%% 
% * hmin= $10^{-5}$
% * $hmax=\displaystyle\frac{(intv(2)-intv(1))}{50}$
% * TOL=0.01;
%% Práctica 3: FSAL (Euler-Euler-mejorado (1-2))
% $$      \begin{array}{l|l}        \begin{array}{l}          {0}	\\          
% {1}        \end{array}        &	         \begin{array}{ll}          0 & 0 \\          
% 1 & 0        \end{array}        \\        \hline        y_{n+1}&        \begin{array}{cc}          
% 1& 0        \end{array}        \\        \hline        \hat{y}_{n+1}&        
% \begin{array}{cc}          \frac{1}{2}& \frac{1}{2}        \end{array}      
% \end{array}$$
%% 
% * avanza con el método de Euler y
% * estima con el método de Euler mejorado.
%% 
% Implementa dicho método llamando la funcion *mieuler12fsal* con la sintaxis 
%%
% 
%  function  [t,y,ev,hchng_vec,err_vec]=mieuler12fsal(f,intv,y0,TOL,hmin,hmax)
%
%% 
% * hmin= $10^{-5}$
% * $hmax=\displaystyle\frac{(intv(2)-intv(1))}{50}$
% * TOL=0.01
% * Usad la propiededad *FSAL*
%% Práctica 4 (Método adaptativo de Dormand-Prince FSAL)
% Su tabla está dada por 
% 
% $$\begin{array}{c|ccccccc}0 \\\displaystyle \frac{1}{5} & \displaystyle \frac{1}{5}               
% \\	 \displaystyle \frac{3}{10}     & \displaystyle \frac{3}{40} & \displaystyle 
% \frac{9}{40} \\	 \displaystyle \frac{4}{5} &       \displaystyle \frac{44}{45} 
% & -\displaystyle \frac{56}{15} & \displaystyle \frac{32}{9} \\ \displaystyle 
% \frac{8}{9} & \displaystyle \frac{19372}{6561}&  -\displaystyle \frac{25360}{2187} 
% & \displaystyle \frac{64448}{6561}& -\displaystyle \frac{212}{729} \\1 &      
% \displaystyle \frac{9017}{3168}& -\displaystyle \frac{355}{33}& \displaystyle 
% \frac{46732}{5247}&     \displaystyle \frac{49}{176}&  -\displaystyle \frac{5103}{18656}  
% \\	1 &      \displaystyle \frac{35}{384}&  0&  \displaystyle \frac{500}{1113}& 
% \displaystyle \frac{125}{192}&     -\displaystyle \frac{2187}{6784}& \displaystyle 
% \frac{11}{84}	\\\hline y_{1}&	  \displaystyle \frac{35}{384}& 0& \displaystyle 
% \frac{500}{1113}& \displaystyle \frac{125}{192}& -\displaystyle \frac{2187}{6784}& 
% \displaystyle \frac{11}{84}&  0	\\\hline\widehat  y_{1}&  \displaystyle \frac{5179}{57600}&0&  
% \displaystyle \frac{7571}{16695} & \displaystyle \frac{393}{640} &\displaystyle 
% \frac{92097}{339200}  & \displaystyle \frac{187}{2100} &  \displaystyle \frac{1}{40}	
% \\\end{array}$$
% 
% Implementa dicho método llamando la funcion *midp45sal* con la sintaxis 
%%
% 
%  function  [t,y,ev,hchng_vec,err_vec]=midp45fsal(f,intv,y0,TOL,hmin,hmax)
%
%% 
% * hmin= $10^{-5}$
% * $hmax=\displaystyle\frac{(intv(2)-intv(1))}{50}$
% * TOL=0.01
% * Usad la propiededad *FSAL*
%% Aplicación
% Práctica 5 (Solución que explota)
% Considera el PVI
% 
% $$  \begin{cases}    x'(t)&=x^2(t)\\     x(0)&=1    \end{cases} $$
% 
% La solución exacta es
% 
% $$   x(t)=\frac{1}{1-t}$$
% 
% que es no acotada cuando $t \to 1$.
%% 
% * Usando el método de _Euler explicito_  resuelve el problema en el intervalo  
% $[0 \quad 2]$.
% * Utiliza ahora los 4 metodos adaptativos. 
% * ¿Qué sucede cerca de la discontinuidad que aparece en $t=1$?
%% 
% *Solución*



%% 
% % Práctica 6 (Ecucacion rigida) 
% Considerar el siguiente sistema 
% 
% $$  y^{\prime}(t)  =  Ay(t) + B(t) \quad t\in [0,10]$$
% 
% $$  \left(   A=  \begin{array}{cc}    -2 & 1\\    998 & -999   \end{array}   
% \right)  \quad  B(t)=\left(   \begin{array}{c}    2\sin(t)\\    999(\cos(t)-\sin(t))  
% \end{array} \right)$$
% 
% $$          y(0)=          \left(             \begin{array}{c}              
% 2\\              3            \end{array}          \right)$$
% 
% La solución exacta es:
% 
% $$  y=2e^{-t}   \left(     \begin{array}{c}      1\\      1    \end{array}  
% \right)  +  \left(     \begin{array}{c}      \sin(t)\\      \cos(t)    \end{array}  
% \right)$$
% 
% Haz un diagrama de eficiencia en la misma manera como en la *hoja1,* con las 
% siguientes diferencias.
%% 
% * Para hacer un diagrama de eficiencia para un método adaptativo cambia la 
% tolerancia, empezando con $TOL_{initial}=0.01$ y repite el calculo con $TOL_{nuevo}=  
% TOL/2$.
% * comparando el método (con paso fijo) del trapecio (con Newton) con  *mieuler12.m* 
% y *mieuler21.m*.
%% 
% *Solucion*


%% Apendice: Las funciones

function [t,y,ev,hchng_vec,err_vec,err_vec2]=mirkfehlb45(f,intv,y0,TOL,hmin,hmax,fac,facmax)
disp('H7: file: UB')
end

%%
function [t,y,ev,hchng_vec,err_vec]=mieuler21(f,intv,y0,TOL,hmin,hmax,fac,facmax)
disp('H7: file: UB')
end
%% 
% 
function [t,y,ev,hchng_vec,err_vec]=mieuler12fsal(f,intv,y0,TOL,hmin,hmax,fac,facmax)
disp('H7: file: UB')
end
%% 
% 


function [t,y,ev,hchng_vec,err_vec,err_vec2]=midp45fsalf,intv,y0,TOL,hmin,hmax,fac,facmax)
disp('H7: file: UB')
end