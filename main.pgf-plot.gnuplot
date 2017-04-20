set table "main.pgf-plot.table"; set format "%.5f"
set format "%.7e";; set contour base; set cntrparam levels 15 set style data lines; set sgrid3d 50, 50; splot "test.csv" u 1:2:3 w l ; 
