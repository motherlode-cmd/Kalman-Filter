splot "plot/filter.txt" u 1:2:3 w d t "filter" , \
    "plot/mls.txt" u 1:2:3 w d t "mls" , \
    "plot/real.txt" u 1:2:3 w lines t "real"
pause -1