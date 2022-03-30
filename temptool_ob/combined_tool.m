%Combine effect of fric+rub in tool
a=0:0.1:1.1; b=0:0.05:0.5;
m=1;t=[];
for i=1:length(a)
    s=a(i);
    for j=1:length(b)
      v=b(j);
      t(i,j)=integ1(s,v);
    end
end
disp(t);
[c,h]=contour(b,a,t,10,'ShowText','on')
h.LevelList=round(h.LevelList,0)  %rounds levels to 3rd decimal place
clabel(c,h)
colorbar
set(gca,'YDir','reverse')  
axis equal

function theta=integ1(a,b)
    fun1 = @(x,y) 754*x.*(1./((a-x).^2+y.^2+b.^2).^0.5+1./((a-2.2+x).^2+y.^2+b.^2).^0.5);
    q1 = integral2(fun1,0,0.98,-2.5,2.5);
    fun11 = @(x,y) 739*(1./((a-x).^2+y.^2+b.^2).^0.5+1./((a-2.2+x).^2+y.^2+b.^2).^0.5);
    q11 = integral2(fun11,0.98,1.1,-2.5,2.5);
    fun2 = @(x,y) 754*x.*(1./((a-x).^2+y.^2+b.^2).^0.5+1./((a-2.2+x).^2+y.^2+b.^2).^0.5).*(x/1.1).^0.26;
    q2 = integral2(fun2,0,0.98,-2.5,2.5);
    fun22 = @(x,y) 738*(1./((a-x).^2+y.^2+b.^2).^0.5+1./((a-2.2+x).^2+y.^2+b.^2).^0.5).*(x/1.1).^0.26;
    q22 = integral2(fun22,0.98,1.1,-2.5,2.5);
    fun3 = @(x,y) 753*x.*(1./((a-x).^2+y.^2+b.^2).^0.5+1./((a-2.2+x).^2+y.^2+b.^2).^0.5).*(x/1.1).^15;
    q3 = integral2(fun3,0,0.98,-2.5,2.5);
    fun33 = @(x,y) 738*(1./((a-x).^2+y.^2+b.^2).^0.5+1./((a-2.2+x).^2+y.^2+b.^2).^0.5).*(x/1.1).^15;
    q33 = integral2(fun33,0.98,1.1,-2.5,2.5);
    
    fun4 = @(y,l) 0.085.*(1./((a-0.1*l-1.1).^2+y.^2+(b+l).^2).^0.5+1./((a-0.1*l-1.1).^2+y.^2+(b-l).^2).^0.5);
    q4 = integral2(fun4,-2.5,2.5,0,0.08);

    theta=0.73*(q1+q11)-0.66*(q2+q22)-0.66*(q3+q33)+2287*q4;
end