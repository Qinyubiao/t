

%　aはx軸の座標，ｂはy軸の座標，ｔ(i,j)はx-y平面つまり二次元切削の工具側温度数値 friction heat
a=0:0.1:2; b=-3:0.1:3;
m=1;t=[];
for i=1:length(a)
    s=a(i);
    for j=1:length(b)
      v=b(j);
      t(i,j)=integ1(s,v);
    end
end
disp(t);
contour(b,a,t,'ShowText','on')

%　single point of temperature by friction motion in tool side
function theta=integ1(a,b)
    fun1 = @(x) abs(-1.*log(-1.*b - 2.5 + 0.1.*x + sqrt(a.^2 - 2.*a.*x + b.^2 + 6.25 + 5.0.*b - 0.2*b.*x + 1.01*x.^2 - 0.50*x))) + abs(log(2.5 + 0.1*x - 1.*b + sqrt(a.^2 - 2.*a.*x + b.^2 + 6.25 - 5.0*b - 0.2*b.*x + 1.01*x.^2 + 0.50*x)));
    q1 = integral(fun1,0,1.2);
    fun2 = @(x) (abs(-1.*log(-1.*b - 2.5 + 0.1.*x + sqrt(a.^2 - 2.*a.*x + b.^2 + 6.25 + 5.0.*b - 0.2*b.*x + 1.01*x.^2 - 0.50*x))) + abs(log(2.5 + 0.1*x - 1.*b + sqrt(a.^2 - 2.*a.*x + b.^2 + 6.25 - 5.0*b - 0.2*b.*x + 1.01*x.^2 + 0.50*x)))).*(x/1.2).^0.26;
    q2 = integral(fun2,0,1.2);
    fun3 = @(x) abs(-1.*log(-25. - 10.*b + x + sqrt(100.*a.^2 + 200.*a.*x + 100.*b.^2 - 20.*b.*x + 101.*x.^2 - 480.*a + 500.*b - 530.*x + 1201.))) + abs(log(25. + x - 10.*b + sqrt(100.*a.^2 + 200.*a.*x + 100.*b.^2 - 20.*b.*x + 101.*x.^2 - 480.*a - 500.*b - 430.*x + 1201.)));
    q3 =integral(fun3,0,1.2);
    fun4 = @(x) (abs(-1.*log(-25. - 10.*b + x + sqrt(100.*a.^2 + 200.*a.*x + 100.*b.^2 - 20.*b.*x + 101.*x.^2 - 480.*a + 500.*b - 530.*x + 1201.))) + abs(log(25. + x - 10.*b + sqrt(100.*a.^2 + 200.*a.*x + 100.*b.^2 - 20.*b.*x + 101.*x.^2 - 480.*a - 500.*b - 430.*x + 1201.)))).*(x/1.2).^0.26;
    q4= integral(fun4,0,1.2);
    theta=331.6*0.64*(q1+q3-q2-q4);
end