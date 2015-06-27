set i 'points' /p1*p10/;
parameters
x(i) 'x coordinates',
y(i) 'y coordinates';
* fill with random data
x(i) = uniform(1,10);
y(i) = uniform(1,10);
variables
a 'x coordinate of center of circle'
b 'y coordinate of center of circle'
r 'radius';
equations
e(i) 'points must be inside circle';
e(i).. sqr(x(i)-a) + sqr(y(i)-b) =l= sqr(r);
r.lo = 0;
model m /all/;
option nlp=minos;
solve m using nlp minimizing r;
