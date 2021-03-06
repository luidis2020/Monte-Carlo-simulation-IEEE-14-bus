%
% https://labs.ece.uw.edu/pstca/pf14/ieee14cdf.txt
%
marcador = 1;
           
%load_variation_ieee14;
load_variation_ieee14BETA;

barras_ipfc =  [  

 1  1   1  1  3  1.060  0.000      0         0        232.4   -16.90   0.0  1.00     0.0      0.0   0.0  0.00  0;
 2  2   1  1  2  1.045  0.000    21.7       12.7      40.00    42.40   0.0  1.00    50.0    -40.0   0.0  0.00  0;
 3  3   1  1  2  1.010  0.000    94.2       19.0       0.00    23.40   0.0  1.00    40.0      0.0   0.0  0.00  0;
 4  4   1  1  0  1.000  0.000    47.8       -3.9       0.00     0.00   0.0  1.00     0.0      0.0   0.0  0.00  0;
 5  5   1  1  0  1.000  0.000     7.6        1.6       0.00     0.00   0.0  1.00     0.0      0.0   0.0  0.00  0;
 6  6   1  1  2  1.070  0.000    11.2        7.5       0.00    12.00   0.0  1.00    24.0     -6.0   0.0  0.00  0;
 7  7   1  1  0  1.000  0.000      0         0         0.00     0.00   0.0  1.00     0.0      0.0   0.0  0.00  0;
 8  8   1  1  2  1.090  0.000      0         0         0.00    17.40   0.0  1.00    24.0     -6.0   0.0  0.00  0;
 9  9   1  1  0  1.000  0.000    29.5       16.6       0.00     0.00   0.0  1.00     0.0    000.0   0.0  0.19  0;
10 10   1  1  0  1.000  0.000     9.0        5.8       0.00     0.00   0.0  1.00     0.0    000.0   0.0  0.00  0;
11 11   1  1  0  1.000  0.000     3.5        1.8       0.00     0.00   0.0  1.00     0.0    000.0   0.0  0.00  0;
12 12   1  1  0  1.000  0.000     6.1        1.6       0.00     0.00   0.0  1.00     0.0    000.0   0.0  0.00  0;
13 13   1  1  0  1.000  0.000    13.5        5.8       0.00     0.00   0.0  1.00     0.0    000.0   0.0  0.00  0;
14 14   1  1  0  1.000  0.000    14.9        5.0       0.00     0.00   0.0  1.00     0.0    000.0   0.0  0.00  0;
];

barras_ipfc(:,8:9)=vector_PQ;
ramos_ipfc = [  
	1	2  1  1 1 0	0.03876	0.11834	0.0264	0  0  0  0 0  1.0   0.0  0.0  0.0  0.0  0.0  0.0;
	1	2  1  1 1 0	0.03876	0.11834	0.0264	0  0  0  0 0  1.0   0.0  0.0  0.0  0.0  0.0  0.0;
	1	5  1  1 1 0	0.05403	0.22304	0.0492	0  0  0  0 0  1.0   0.0  0.0  0.0  0.0  0.0  0.0;
	2	3  1  1 1 0	0.04699	0.19797	0.0438	0  0  0  0 0  1.0   0.0  0.0  0.0  0.0  0.0  0.0;
	2	4  1  1 1 0	0.05811	0.17632	0.034	0  0  0  0 0  1.0   0.0  0.0  0.0  0.0  0.0  0.0;
	2	5  1  1 1 0	0.05695	0.17388	0.0346	0  0  0  0 0  1.0   0.0  0.0  0.0  0.0  0.0  0.0;
	3	4  1  1 1 0	0.06701	0.17103	0.0128	0  0  0  0 0  1.0   0.0  0.0  0.0  0.0  0.0  0.0;
	4	5  1  1 1 0	0.01335	0.04211	0	    0  0  0  0 0  1.0   0.0  0.0  0.0  0.0  0.0  0.0;
	4	7  1  1 1 0	0	    0.20912	0	    0  0  0  0 0  0.950  0.0  0.0  0.0  0.0  0.0  0.0;
	4	9  1  1 1 0	0	    0.55618	0	    0  0  0  0 0  0.961  0.0  0.0  0.0  0.0  0.0  0.0;
	5	6  1  1 1 0	0	    0.25202	0	    0  0  0  0 0  0.932  0.0  0.0  0.0  0.0  0.0  0.0;
	6	11 1  1 1 0	0.09498	0.1989	0	    0  0  0  0 0  1.0   0.0  0.0  0.0  0.0  0.0  0.0;
	6	12 1  1 1 0	0.12291	0.25581	0	    0  0  0  0 0  1.0   0.0  0.0  0.0  0.0  0.0  0.0;
	6	13 1  1 1 0	0.06615	0.13027	0	    0  0  0  0 0  1.0   0.0  0.0  0.0  0.0  0.0  0.0;
	7	8  1  1 1 0	0	    0.17615	0	    0  0  0  0 0  1.0   0.0  0.0  0.0  0.0  0.0  0.0;
	7	9  1  1 1 0	0	    0.11001	0	    0  0  0  0 0  1.0   0.0  0.0  0.0  0.0  0.0  0.0;
	9	10 1  1 1 0	0.03181	0.0845	0	    0  0  0  0 0  1.0   0.0  0.0  0.0  0.0  0.0  0.0;
	9	14 1  1 1 0	0.12711	0.27038	0	    0  0  0  0 0  1.0   0.0  0.0  0.0  0.0  0.0  0.0;
	10	11 1  1 1 0	0.08205	0.19207	0	    0  0  0  0 0  1.0   0.0  0.0  0.0  0.0  0.0  0.0;
	12	13 1  1 1 0	0.22092	0.19988	0	    0  0  0  0 0  1.0   0.0  0.0  0.0  0.0  0.0  0.0;
	13	14 1  1 1 0	0.17093	0.34802	0	    0  0  0  0 0  1.0   0.0  0.0  0.0  0.0  0.0  0.0;
];

