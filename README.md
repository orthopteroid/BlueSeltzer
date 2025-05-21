# BlueSeltzer
(refreshing code for old skool shape recognition... No AI in here!)

This little demo loads an image and tries to identify the shapes by doing some
simple edge filtering, point-cloud collection, and finally by constructing and
interpreting some basic metrics. Those metrics are:

1. Angular Spectrum (sidedness) - From the center of the bounding-box the point-cloud is
split radially into a number of pie-segments and the points in each segment are
counted and analyzed against points in adjacent segments and overall. This is actually done using an FFT
which gives a radial frequency spectrum. The goal here is to guess a geometric shape from the
radial frequencies [1,2].

2. Circularity Measure - This is calculated using the
the average point distance
from the center of a point-cloud's bounding-box, 
and the variance of that distance. Small values for this
variance would indicate a higher chance we're looking at a cloud of points that
resemble a circle.

3. Radial Density Factor - This is the slope of the radial density of the
point cloud. The first draft code does a simple line through the endpoints
of the function but a better approach would be to do a least-squares
fitting. Numbers near zero show uniformity, +ve numbers show increasing
density towards the outside edge and negative numbers show decreasing
density.

4. Radial Spectrum (dartboard) - This analysis derives a count-spectrum for each angle
around each object and then show the average frequency circularly around each object
in each angle 'wedge' and the stdev of the frequencies around the object. This analysis
is probably only useful in analyzing radially periodic patterns, like spirals. In this case
the average frequency would probably relate to the ratio between the size of the radial bins
and the size of the object under analysis, so there should be some tuning here possibly. For
patterns that are strongly radially periodic and radially symmetrical the stdev coefficient
should likely be below 1 and hopefully close to 0.

And now for some tests!

## A Circle?

![s](circle.jpg)
```
circle.jpg width 474 height 332
object  x       y       sided   circ    rdens   rfavg   rfdev
1       228     168     6       0.927   +0.060  2.000   1.925
```

## A Triangle?

![s](tri.jpg)
```
tri.jpg width 474 height 332
object  x       y       sided   circ    rdens   rfavg   rfdev
1       260     157     3       0.505   +0.270  3.062   1.708
```

## Shapes Everywhere!

It looks like object 2 isn't recognized. I'm not sure what is happening there...

![s](circ-tri-squ-pent.jpg)
```
circ-tri-squ-pent.jpg width 474 height 332
object  x       y       sided   circ    rdens   rfavg   rfdev
1       98      86      5       0.831   +0.003  2.438   1.830
3       269     89      4       0.793   +0.676  2.625   1.876
4       411     56      2       0.901   +0.000  2.000   2.121
5       416     119     6       0.928   +0.000  2.000   2.121
6       358     234     3       0.563   +0.058  3.000   1.655
7       242     197     3       0.768   +0.000  2.000   2.121
8       130     252     5       0.862   +1.000  4.000   1.768
```

## Tricky shapes...

These are a work in progress. The first stage has been in ensuring the point-clouds
are unified. Later stages to the analysis might be a radial-density,
linear-segment-length or fractal shape analysis to disambiguate stars from spirals.

![s](uu.jpg)
```
uu.jpg width 474 height 332
object  x       y       sided   circ    rdens   rfavg   rfdev
1       204     164     2       0.288   +0.039  3.469   -nan(ind)
```

![s](stars.jpg)
```
stars.jpg width 474 height 332
object  x       y       sided   circ    rdens   rfavg   rfdev
1       429     39      4       0.675   +0.000  1.562   -nan(ind)
2       43      45      4       0.666   +0.000  1.750   -nan(ind)
4       256     170     14      0.452   -0.258  1.938   -nan(ind)
8       77      257     5       0.690   -0.840  1.750   -nan(ind)
```

![s](spiral.jpg)

For the sprial, the radial frequency deviation is below 1. This is probably indicating
a consistient radial frequency.

```
spiral.jpg width 474 height 332
object  x       y       sided   circ    rdens   rfavg   rfdev
1       212     160     2       0.605   -0.083  2.000   0.445
```

## Dubious References

1. "Radial frequency patterns describe a small and perceptually distinct subset of all possible planar shapes" https://www.sciencedirect.com/science/article/pii/S0042698918302219
2. "The role of local features in shape discrimination of contour- and surface-defined radial frequency patterns at low contrast" https://www.sciencedirect.com/science/article/pii/S0042698911003555
3. "Boundaries and Coastlines: The Fractals Paradox" https://www.georgeszpiro.com/22-fractal-coastline