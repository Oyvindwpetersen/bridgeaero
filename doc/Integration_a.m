%%

clc
clear all
close all

syms x L q1 q2 a11 a12 a13 a21 a22 a23 a31 a32 a33

x_vec_poly=[1 x x^2 x^3].'

H1=[1 0/L -3/L^2 2/L^3]*x_vec_poly
H2=[0 -1 2/L -1/L^2]*x_vec_poly
H3=[0 0/L 3/L^2 -2/L^3]*x_vec_poly
H4=[0 0/L 1/L -1/L^2]*x_vec_poly

N1=1-x/L;
N2=x/L;

N_r_lin=[...
N1 0 0 , N2 0 0 ;... % U2
0 N1 0 , 0 N2 0 ;... % U3
0 0 N1 , 0 0 N2 ;... % UR1
]

N_f_lin=[...
0 0 0 ;... % U1
N1 0 0 ;... % U2
0 N1 0 ;... % U3
0 0 N1 ;... % UR1
0 0 0 ;... % UR2
0 0 0 ;... % UR3
...
0 0 0 ;... % U1
N2 0 0 ;... % U2
0 N2 0 ;... % U3
0 0 N2 ;... % UR1
0 0 0 ;... % UR2
0 0 0 ;... % UR3
]

a=[...
    a11 a12 a13 ; ...
    a21 a22 a23 ; ...
    a31 a32 a33 ; ...
]

% T_int=int(N_f_lin*N_r_lin,x,[0 L])
T_int=int(N_f_lin*a*N_r_lin,x,[0 L])

%%
clc

N_f_herm=[...
0 0 0 ;... % U1
H1 0 0 ;... % U2
0 H1 0 ;... % U3
0 0 H1 ;... % UR1
0 H2 0 ;... % UR2
-H2 0 0 ;... % UR3
...
0 0 0 ;... % U1
H3 0 0 ;... % U2
0 H3 0 ;... % U3
0 0 H3 ;... % UR1
H4 0 0 ;... % UR2
0 -H4 0 ;... % UR3
]

% T_int2=int(N_mat_2*N_mat_1,x,[0 L])
T_int2=int(N_f_herm*a*N_r_lin,x,[0 L])


