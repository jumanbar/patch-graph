## Motivo:

Análisis de percolación de paisaje basado en la ubicación y distancias entre parches.

### Función principal:

> patchCluster

### Requiere: 

* vegan
* igraph

## Ejemplo

Para analizar un grupo de parches distribuidos uniformemente en el paisaje.

```R
x <- patchCluster(100)
```

El objeto x es una lista con varias cosillas.

### Salida gráfica:


![](https://github.com/jumanbar/patch-graph/raw/master/runif1000.png)

Para generar esto, correr:

```R
patchCluster(1000)
```

La salida gráfica muestra (por panel):

1.  Los parches
2.  El Minimum Spaning Tree (MST) obtenido a partir de los parches (en base a un grafo cuyos links tienen un peso == distancia geográfica).
3.  Distribución de los pesos (distancias) de los links en el MST
4.  Gráfica de percolación:
    *   Eje x: distancia de movimiento (d_mov)
    *   Eje y:
        1.  Línea gris: número de componentes del grafo G construido así:

            > G(i,j) = 1 <==> d(i,j) <= d_mov

        2.  Valor esperado de la cantidad de parches a los que tiene acceso un individuo.

