function theta=integ1(a,b)
fun1 = @(x,y) 0.64*((1./(((a-x).^2+(b-y).^2).^0.5)+1./(((a-2.4+x).^2+(b-y).^2).^0.5)));
q1 = integral2(fun1,0,1.2,-2.5-0.1*x,2.5-0.1*x);
fun2 = @(x,y) 0.64*((1./(((a-x).^2+(b-y).^2).^0.5)+1./(((a-2.4+x).^2+(b-y).^2).^0.5))).*(x/1.2).^0.26;
q2 = integral2(fun2,0,1.2,-2.5-0.1*x,2.5-0.1*x);
theta=2070*(q1-q2);
end