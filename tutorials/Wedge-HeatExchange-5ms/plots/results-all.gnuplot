set datafile separator ','
set ylabel "Temp, K" # label for the Y axis
set xlabel 'Axial Distance' # label for the X axis
set xrange [0.0:0.5]
set x2range [0.5:0.0]
set ytics 100, 50
set yrange [300:1000]
set title 'DRM: 1D vs 3D, Temp vs. axial distance 5 m/s'   
set key outside
set key font ",8"
#set key right #top
set term pngcairo dashed
#set termoption dashed
set output "DRM-HXOnly-DUO-vs-OpenFOAM-vs-PBR-Wedge-U5.0ms.png"
plot 'PBR-molefracs00001.csv' using 1:2  with lines axis x1y1 lc rgb "dark-blue" title 'Temp-PBR', \
'Daten-Mittelwerte-U5.0.csv' using 1:8 with linespoint axis x2y1 lc rgb "dark-blue" title 'T-Average-DUO', \
'Daten-Symmetrieachse-U5.0.csv' using 1:2 with lines axis x2y1 dt 2 lc rgb "dark-blue" title 'T-Center-DUO', \
'results.csv'using 1:2 with lines axis x1y1 dt 2 lc rgb "dark-red" title 'T-CHTMRF-2D', \
