# Las siguientes reglas son obligatorias

1.  El alumno encuentra  para cada práctica un documento en formato **mlx**

2.  Tiene que bajar dicho fichero a su ordenador y ponerlo en una
    carpeta: cada práctica en una carpeta.

3.  Las carpetas y los ficheros no deben tener **espacio** ni
    **acentos** ni **-**, (**\_** están permitidos).

4.  Tiene que poner Nombre, Apellido y DNI en los campos correspondientes.

5.  El alumno tiene que responder a los ejercicios poniendo los
    comandos de matlab relevantes, en un campo, usando el comando
    *insert code* (control+e). Siempre al principio del campo tiene que poner
    
    `disp('Hoja X: Ej Y Nombre y apellido')`

6.  Cada campo tiene que terminar un un **section break** (Control+Alt+Enter).

7.  El alumno tiene que implementar ciertos algoritmos, cuyo
    pseudocódigo se encuentra en el CV. Dicha implementación debe
    ser una **function** que tiene que **siempre** estar al final del
    fichero mlx, en el apéndice. La sintaxis es algo así:
    
    `[t,y]=mieuler(f,intv,y0,N)`
    
    `end`
    
    El end al final es muy importante. También dicha función debe
    tener
    `disp('Hoja X: func: Euler Nombre y Apellido')`

8.  El alumno tiene que entregar el fichero **mlx** pero también
    salvar el fichero con la extensión **m** pero con el sufijo text,
    por ejemplo: hoja7.mlx y hoja7<sub>text.m</sub>.

9.  Las entregas de prácticas no cuentan para la nota pero, si el
    alumno ha entregado prácticas a tiempo, las puede usar en el
    examen. Si no entrega prácticas no puede llevar nada al examen.


# Prácticas de Matlab (Distribución temporal de las hojas de Matlab)
<table border="2" cellspacing="0" cellpadding="6" frame="border" rules="all">
<colgroup>
<col  class="org-right" />
</colgroup>
<colgroup>
<col  class="org-left" />
</colgroup>
<colgroup>
<col  class="org-left" />
</colgroup>
<colgroup>
<col  class="org-left" />
</colgroup>
<colgroup>
<col  class="org-left" />
</colgroup>
<thead>
<tr>
<th scope="col" class="org-right">Hoja</th>
<th scope="col" class="org-left">Contenido</th>
<th scope="col" class="org-left">Publicación</th>
<th scope="col" class="org-left">Objetivo</th>
<th scope="col" class="org-left">Entrega</th>
</tr>
<tr>
<th scope="col" class="org-right">&#xa0;</th>
<th scope="col" class="org-left">&#xa0;</th>
<th scope="col" class="org-left">aproximadamente</th>
<th scope="col" class="org-left">&#xa0;</th>
<th scope="col" class="org-left">&#xa0;</th>
</tr>
</thead>
<tbody>
<tr>
<td class="org-right">1</td>
<td class="org-left">Ode45</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp"> 24.01.2023 </span></span></td>
<td class="org-left">aprender ode45, sprintf, title</td>
<td class="org-left">opcional</td>
</tr>
</tbody>
<tbody>
<tr>
<td class="org-right">2</td>
<td class="org-left">Bucles simples</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp"> 24.01.2023 </span></span></td>
<td class="org-left">aprender a usar bucles sin índices</td>
<td class="org-left">opcional</td>
</tr>
</tbody>
<tbody>
<tr>
<td class="org-right">3</td>
<td class="org-left">Euler/RK</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp"> 30.01.2023 </span></span></td>
<td class="org-left">Implementar métodos monopaso</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp"> 13.02.2023 </span></span></td>
</tr>
<tr>
<td class="org-right">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">explícitos</td>
<td class="org-left">&#xa0;</td>
</tr>
</tbody>
<tbody>
<tr>
<td class="org-right">4</td>
<td class="org-left">Diagrama de eficiencia</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp"> 06.02.2023 </span></span></td>
<td class="org-left">Construir diagramas para calcular</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp"> 20.02.2023 </span></span></td>
</tr>
<tr>
<td class="org-right">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">el orden de un método</td>
<td class="org-left">&#xa0;</td>
</tr>
</tbody>
<tbody>
<tr>
<td class="org-right">5</td>
<td class="org-left">Métodos implícitos</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp"> 20.02.2023 </span></span></td>
<td class="org-left">Implementar métodos monopaso</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp"> 20.03.2023 </span></span></td>
</tr>
<tr>
<td class="org-right">&#xa0;</td>
<td class="org-left">por Newton y punto fijo</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">implícitos mediante</td>
<td class="org-left">&#xa0;</td>
</tr>
<tr>
<td class="org-right">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">dos esquemas: punto fijo y Newton;</td>
<td class="org-left">&#xa0;</td>
</tr>
<tr>
<td class="org-right">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">diagramas de</td>
<td class="org-left">&#xa0;</td>
</tr>
<tr>
<td class="org-right">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">eficiencia</td>
<td class="org-left">&#xa0;</td>
</tr>
<tr>
<td class="org-right">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>
</tbody>
<tbody>
<tr>
<td class="org-right">6</td>
<td class="org-left">Métodos Multipaso</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp"> 20.03.2023 </span></span></td>
<td class="org-left">Implementar un método multipaso</td>
<td class="org-left">&#xa0;</td>
</tr>
<tr>
<td class="org-right">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">implícito mediante</td>
<td class="org-left">&#xa0;</td>
</tr>
<tr>
<td class="org-right">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">Newton</td>
<td class="org-left">&#xa0;</td>
</tr>
</tbody>
<tbody>
<tr>
<td class="org-right">7</td>
<td class="org-left">Adaptativos</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp"> 27.03.2023 </span></span></td>
<td class="org-left">Implementar métodos monopaso con</td>
<td class="org-left">&#xa0;</td>
</tr>
<tr>
<td class="org-right">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">paso adaptativo</td>
<td class="org-left">&#xa0;</td>
</tr>
<tr>
<td class="org-right">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">usando h<sub>opt</sub> y pares encajados</td>
<td class="org-left">&#xa0;</td>
</tr>
</tbody>
<tbody>
<tr>
<td class="org-right">8</td>
<td class="org-left">Locus</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp"> 10.04.2023 </span></span></td>
<td class="org-left">usando cálculo simbólico, pintar</td>
<td class="org-left">opcional</td>
</tr>
<tr>
<td class="org-right">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">regiones de</td>
<td class="org-left">&#xa0;</td>
</tr>
<tr>
<td class="org-right">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">estabilidad</td>
<td class="org-left">&#xa0;</td>
</tr>
</tbody>
<tbody>
<tr>
<td class="org-right">9</td>
<td class="org-left">Disparo</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp"> 17.04.2023 </span></span></td>
<td class="org-left">usando RK, implementar el disparo</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp"> 28.04.2023 </span></span></td>
</tr>
<tr>
<td class="org-right">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">lineal para</td>
<td class="org-left">&#xa0;</td>
</tr>
<tr>
<td class="org-right">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">condiciones de Dirichlet y Neuman</td>
<td class="org-left">&#xa0;</td>
</tr>
</tbody>
</table>
