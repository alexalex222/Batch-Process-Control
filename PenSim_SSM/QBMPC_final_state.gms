$set matout "'matsol2.gdx', x, JMPC, Quality";
*OPTION LIMROW=10
*OPTION LIMCOL=10

SETS
T 'time steps to end of batch'
M 'measured variables from mult-model'
P 'variables predicted using mult-model'
L 'lags in the ARX model'
Q 'qualities'
max_min 'maximum or minimum for constraint'
U 'inputs (set used only for constraints)'




$gdxin "QBMPCv2sets.gdx"
$load T
$load M
$load P
$load L
$load Q
$load max_min
$load U

$gdxin




PARAMETERS
data(T, M)
Coeff_ARX(M, P, L)
Coeff_Quality(M, Q, T)
qual_set(Q)
*input_const(U, max_min)
input_const(T,U)
q_pred_mean(Q)
q_pred_std(Q)

;

$gdxin "MPCinitialize.gdx"
$load data
$load Coeff_ARX
$load Coeff_Quality
$load qual_set
$load input_const
$load q_pred_mean
$load q_pred_std
$gdxin


VARIABLES
X(T,M) 'All trajectories measured by the dynamic model and used by the quality model'
*Y(T,P) 'Measurements'
Quality(Q) 'qualities'
JMPC 'MPC objecive'
Dummy


;

* Initialize variables
X.l(T,M) = data(T,M);

* fix the initial two measurements

x.fx('1','1') = data('1', '1');

x.fx('1','2') = data('1', '2');

x.fx('1','3') = data('1', '3');

x.fx('1','4') = data('1', '4');

x.fx('1','5') = data('1', '5');

x.fx('1','6') = data('1', '6');

x.fx('2','1') = data('2', '1');

x.fx('2','2') = data('2', '2');

x.fx('2','3') = data('2', '3');

x.fx('2','4') = data('2', '4');

x.fx('2','5') = data('2', '5');

x.fx('2','6') = data('2', '6');


* fix the initial inputs

x.fx('1','7') = data('1','7');

x.fx('1','8') = data('1','8');

x.fx('1','9') = data('1','9');

x.fx('1','10') = data('1','10');

* the contraints for inputs
*loop (T $(ord(T)>1),

*    x.lo(T, '7') = input_const('1', '1');
*    x.up(T, '7') = input_const('1', '2');

*    x.lo(T, '8') = input_const('2', '1');
*    x.up(T, '8') = input_const('2', '2');

*    x.lo(T, '9') = input_const('3', '1');
*    x.up(T, '9') = input_const('3', '2');

*    x.lo(T, '10') = input_const('4', '1');
*    x.up(T, '10') = input_const('4', '2');
*);

    x.lo(T, '7') = input_const(T, '1');
    x.up(T, '7') = input_const(T, '5');

    x.lo(T, '8') = input_const(T, '2');
    x.up(T, '8') = input_const(T, '6');

    x.lo(T, '9') = input_const(T, '3');
    x.up(T, '9') = input_const(T, '7');

    x.lo(T, '10') = input_const(T, '4');
    x.up(T, '10') = input_const(T, '8');


* Predict future ouputs

loop ( T $(ord(T)>2),

    x.l(T,'1') =  sum(M, x.l(T-1,M)*Coeff_ARX(M, '1', '1')) + sum(M, x.l(T-2,M)*Coeff_ARX(M, '1', '2'));

    x.l(T,'2') =  sum(M, x.l(T-1,M)*Coeff_ARX(M, '2', '1')) + sum(M, x.l(T-2,M)*Coeff_ARX(M, '2', '2'));

    x.l(T,'3') =  sum(M, x.l(T-1,M)*Coeff_ARX(M, '3', '1')) + sum(M, x.l(T-2,M)*Coeff_ARX(M, '3', '2'));

    x.l(T,'4') =  sum(M, x.l(T-1,M)*Coeff_ARX(M, '4', '1')) + sum(M, x.l(T-2,M)*Coeff_ARX(M, '4', '2'));

    x.l(T,'5') =  sum(M, x.l(T-1,M)*Coeff_ARX(M, '5', '1')) + sum(M, x.l(T-2,M)*Coeff_ARX(M, '5', '2'));

    x.l(T,'6') =  sum(M, x.l(T-1,M)*Coeff_ARX(M, '6', '1')) + sum(M, x.l(T-2,M)*Coeff_ARX(M, '6', '2'));

);

*loop (T,

*   y.l(T,'1') = x.l(T,'1');
*   y.l(T,'2') = x.l(T,'2');
*   y.l(T,'3') = x.l(T,'3');
*   y.l(T,'4') = x.l(T,'4');
*   y.l(T,'5') = x.l(T,'5');
*   y.l(T,'6') = x.l(T,'6');
*);




*loop (Q,

        quality.l(Q) = sum((T ,M), x.l(T, M)*Coeff_Quality(M, Q, T));



*);





JMPC.l =sum(Q, ((quality.l(Q) - qual_set(Q))/q_pred_std(Q))* ((quality.l(Q) - qual_set(Q))/q_pred_std(Q)));




EQUATIONS


Y1_Mult_Model(T)
Y2_Mult_Model(T)
Y3_Mult_Model(T)
Y4_Mult_Model(T)
Y5_Mult_Model(T)
Y6_Mult_Model(T)

Quality_eqn(Q)

Quality_low(Q)



Objective
;



Y1_Mult_Model(T) $(ord(T)>2) .. x(T, '1') =e=  sum(M, x(T-1,M)*Coeff_ARX(M, '1', '1')) + sum(M, x(T-2,M)*Coeff_ARX(M, '1', '2'));

Y2_Mult_Model(T) $(ord(T)>2) .. x(T, '2') =e=  sum(M, x(T-1,M)*Coeff_ARX(M, '2', '1')) + sum(M, x(T-2,M)*Coeff_ARX(M, '2', '2'));

Y3_Mult_Model(T) $(ord(T)>2) .. x(T, '3') =e=  sum(M, x(T-1,M)*Coeff_ARX(M, '3', '1')) + sum(M, x(T-2,M)*Coeff_ARX(M, '3', '2'));

Y4_Mult_Model(T) $(ord(T)>2) .. x(T, '4') =e=  sum(M, x(T-1,M)*Coeff_ARX(M, '4', '1')) + sum(M, x(T-2,M)*Coeff_ARX(M, '4', '2'));

Y5_Mult_Model(T) $(ord(T)>2) .. x(T, '5') =e=  sum(M, x(T-1,M)*Coeff_ARX(M, '5', '1')) + sum(M, x(T-2,M)*Coeff_ARX(M, '5', '2'));

Y6_Mult_Model(T) $(ord(T)>2) .. x(T, '6') =e=  sum(M, x(T-1,M)*Coeff_ARX(M, '6', '1')) + sum(M, x(T-2,M)*Coeff_ARX(M, '6', '2'));

Quality_eqn(Q) .. Quality(Q) =e= sum((T ,M),  x(T, M)*Coeff_Quality(M, Q, T));

Quality_low(Q) .. Quality(Q) =g= 0;



Objective .. JMPC =e= sum(Q, ((quality(Q) - qual_set(Q))/q_pred_std(Q))* ((quality(Q) - qual_set(Q))/q_pred_std(Q)));
*+.001*sum(T$(ord(T)>3), (x(T,'7') - 0.04)*(x(T,'7')-0.04)+(x(T,'8') - 8)*(x(T,'8') - 8)+(x(T,'9') - 30)*(x(T,'9') - 30)+(x(T,'10') - 0.06)*(x(T,'10') - 0.06));



MODEL QBMPC /all/;

QBMPC.SCALEOPT = 1;

option nlp = ipopt;

solve QBMPC using nlp minimizing JMPC;

execute_unload %matout%;
