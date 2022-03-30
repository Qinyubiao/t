%chip fric
a=0:0.1:1.2;b=0:0.05:0.5;
r=[];
for i=1:length(a)
    s=a(i);
    for j=1:length(b)
        v=b(j);
        r(i,j)=equ5(s,v);
    end
end
disp(r);
contour(b,a,r,'ShowText','on')


function y=equ5(a,b)
f1=@(l) (exp(-347.*(a-l))).*(bessely(0, 347.*((a-l).^2+b.^2).^0.5)+bessely(0, 347.*((a-l).^2+(1-b).^2).^0.5));
p1=integral(f1,0, 1.2);
f2=@(l) ((l./1.2).^0.24).*exp(-347*(a-l)).*(bessely(0, 347.*((a-l).^2+b.^2).^0.5)+bessely(0, 347.*((a-l).^2+(1-b).^2).^0.5));
p2=integral(f2,0, 1.2);

y=14483.7*(0.27*p1+0.7*p2);
end