### LEC_QHULL

Solves the largest empty circle (LEC) problem for a set of 2-d points.
This consists of finding the largest circle that contains none of the
2-d points and is also centred inside the convex hull of the points.
This solution uses the IDL routine QHULL to determine the Delaunay
triagulation and the Veronoi diagram of the set of points.

Example:
```IDL
IDL> x = randomn(seed, 20)
IDL> y = randomn(seed, 20)
IDL> plot, psym=6, x, y, /isotropic, /ynozero
IDL> result = lec_qhull(x, y, /plot_tri, /plot_vor)
IDL> tvcircle, lec.radius, lec.x, lec.y, /data
```
