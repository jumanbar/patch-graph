Motivo:
---------

Funciones para hacer un análisis de percolación de paisaje basado en la ubicación y distancias entre parches.

Función principal:
------------------

__patchCluster__


Requiere:
---------

* vegan
* igraph

Ejemplo:
# Para analizar un grupo de parches distribuidos uniformemente en el paisaje.
x <- patchCluster(100)

Salida gráfica:
---------------

Ejemplo:

![](https://github.com/jumanbar/patch-graph/blob/master/runif1000.png)

La salida gráfica muestra (x panel):
* Los parches
* El Minimum Spaning Tree (MST) obtenido a partir de los parches (en base a un grafo cuyos links tienen un peso == distancia geográfica).
* Distribución de los pesos (distancias) de los links en el MST
* Gráfica de percolación:
    a. Número de componentes del grafo (línea gris) construido así:
           G(i,j) = 1 <==> d(i,j) <= d_movimiento (d_movimiento = eje x).
    b. Valor esperado de la cantidad de parches a los que tiene acceso un individuo con camacidad de
       movimiento = d_movimiento.
