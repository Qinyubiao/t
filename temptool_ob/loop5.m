%　aはx軸の座標，ｂはy軸の座標，ｔ(i,j)はx-z平面つまり二次元切削の工具側温度数値 friction heat
a=0:0.05:1.2; b=-0.5:0.05:0;
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
    fun1 = @(x,y) 1./((a-x).^2+y.^2+b.^2).^0.5+1./((a-2.4+x).^2+y.^2+b.^2).^0.5;
    q1 = integral2(fun1,0,1.2,-2.5,2.5);
    fun2 = @(x,y) (1./((a-x).^2+y.^2+b.^2).^0.5+1./((a-2.4+x).^2+y.^2+b.^2).^0.5).*(x/1.2).^0.26;
    q2 = integral2(fun2,0,1.2,-2.5,2.5);

    theta=331.6*0.64*(q1-q2);
end