$set matout "'matsol2.gdx', X, JMPC, Quality";
OPTION LIMROW=10
OPTION LIMCOL=10

SETS
T 'time steps to end of batch'
M 'measured variables from mult-model'
P 'variables predicted using mult-model'
L 'lags in the ARX model'
Q 'qualities'
max_min 'maximum or minimum for constraint'
U 'inputs (set used only for constraints)'
PC 'number of principal components in quality model' ,
PD 'past data used in SPE calculations for quality model';



$gdxin "QBMPCv2sets.gdx"
$load T
$load M
$load P
$load L
$load Q
$load max_min
$load U
$load PC
$load PD
$gdxin


alias(M, M1);

PARAMETERS
data(T, M)
Coeff_ARX(M,P)
qual_coeffs(P, Q, T)
q_meas_mean(T,P)
q_meas_std(T,P)
q_pred_mean(Q)
q_pred_std(Q)
contrib_past(Q)
qual_set(Q)
input_const(U, max_min)
;

$gdxin "MPCinitialize.gdx"
$load data
$load Coeff_ARX
$load qual_coeffs
$load q_meas_mean
$load q_meas_std
$load q_pred_mean
$load q_pred_std
$load contrib_past
$load qual_set
$load input_const

$gdxin


VARIABLES
X(T,M) 'All trajectories measured by the dynamic model and used by the quality model'
Y(T,P) 'Measurements'
Quality(Q) 'qualities'
JMPC 'MPC objecive'
Dummy


;

* Initialize variables
X.l(T,M) = data(T,M);

* fix the values corresponding to curent batch condition (ie current state and output)
x.fx('1','1') = data('1', '1');

x.fx('1','2') = data('1', '2');

x.fx('1','3') = data('1', '3');

x.fx('1','4') = data('1', '4');




loop (T,
    x.lo(T, '5') = input_const('1', '1');
    x.up(T, '5') = input_const('1', '2');

    x.lo(T, '6') = input_const('2', '1');
    x.up(T, '6') = input_const('2', '2');

    x.lo(T, '7') = input_const('3', '1');
    x.up(T, '7') = input_const('3', '2');
);



loop ( T $(ord(T)>1),

    x.l(T,'1') =  sum(M, x.l(T-1,M)*Coeff_ARX(M, '1'));

    x.l(T,'2') =  sum(M, x.l(T-1,M)*Coeff_ARX(M, '2'));

    x.l(T,'3') =  sum(M, x.l(T-1,M)*Coeff_ARX(M, '3'));

    x.l(T,'4') =  sum(M, x.l(T-1,M)*Coeff_ARX(M, '4'));

)

loop (T,

   y.l(T,'1')= x.l(T,'1');
   y.l(T,'2')= x.l(T,'2');
   y.l(T,'3')= x.l(T,'3');
   y.l(T,'4')= x.l(T,'4');
);




loop (Q,

         Quality.l(Q) = (contrib_past(Q) + sum(T, sum(P,
                 (y.l(T,P)-q_meas_mean(T,P))/q_meas_std(T,P)
                 *qual_coeffs(P, Q, T))));


);



JMPC.l =sum(Q, (quality.l(Q) - qual_set(Q))*(quality.l(Q) - qual_set(Q)));




EQUATIONS

Y1_Mult_Model(T)
Y2_Mult_Model(T)
Y3_Mult_Model(T)
Y4_Mult_Model(T)

Quality_eqn(Q)

Objective
;

Y1_Mult_Model(T) $(ord(T)>1) .. x(T, '1') =e= sum(M, x.l(T-1,M)*Coeff_ARX(M, '1'));

Y2_Mult_Model(T) $(ord(T)>1) .. x(T, '2') =e= sum(M, x.l(T-1,M)*Coeff_ARX(M, '2'));

Y3_Mult_Model(T) $(ord(T)>1) .. x(T, '3') =e= sum(M, x.l(T-1,M)*Coeff_ARX(M, '3'));

Y4_Mult_Model(T) $(ord(T)>1) .. x(T, '4') =e= sum(M, x.l(T-1,M)*Coeff_ARX(M, '4'));

Quality_eqn(Q) .. Quality(Q) =e= (contrib_past(Q) + sum(T, sum(P, (y.l(T,P)-q_meas_mean(T,P))/q_meas_std(T,P)*qual_coeffs(P, Q, T))));



Objective .. JMPC =e= sum(Q, (quality.l(Q) - qual_set(Q))*(quality.l(Q) - qual_set(Q)))
+ .000000*sum(T $(ord(T)>1), (x(T,'5') - 0.04)*(x(T,'5')-0.04)+(x(T,'6') - 8)*(x(T,'6') - 8)+(x(T,'7') - 30)*(x(T,'7') - 30));



MODEL QBMPC /all/;

solve QBMPC using nlp minimizing JMPC;

execute_unload %matout%;
