Motivo:
---------

Funciones para hacer un análisis de percolación de paisaje basado en la ubicación y distancias entre parches.

__Función principal:__

patchCluster

__Requiere:__

* vegan
* igraph

__Ejemplo:__

Para analizar un grupo de parches distribuidos uniformemente en el paisaje.

```R
x <- patchCluster(100)
```

El objeto x es una lista con varias cosillas.

__Salida gráfica:__


![](https://github.com/jumanbar/patch-graph/raw/master/runif1000.png)

Para generar esto:
```R
patchCluster(100)
```

La salida gráfica muestra (x panel):

1.  Los parches
2.  El Minimum Spaning Tree (MST) obtenido a partir de los parches (en base a un grafo cuyos links tienen un peso == distancia geográfica).
3.  Distribución de los pesos (distancias) de los links en el MST
4.  Gráfica de percolación:
    *  Número de componentes del grafo (línea gris) construido así:
       <G(i,j) = 1 <==> d(i,j) <= d_movimiento (d_movimiento = eje x).>
    * Valor esperado de la cantidad de parches a los que tiene acceso un individuo con camacidad de
       movimiento = d_movimiento.
