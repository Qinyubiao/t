%　aはｚ軸の座標，ｂはｘ軸の座標，r(i,j)はx-z平面つまり二次元切削の切りくず側温度数値
a=0.01:0.01:0.12;b=0:0.005:0.06;
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
set(gca,'YDir','reverse')
set(gca,'XDir','reverse') 

function y=equ5(x,z)

f1=@(l) (exp(-347.*(x-l))).*(bessely(0, 347.*((x-l).^2+z.^2).^0.5)+bessely(0, 347.*((x-l).^2+(0.1327-z).^2).^0.5));
p1=integral(f1,0, 0.1209);

f2=@(l) ((l./0.1209).^0.24).*exp(-347*(x-l)).*(bessely(0, 347.*((x-l).^2+z.^2).^0.5)+bessely(0, 347.*((x-l).^2+(0.1327-z).^2).^0.5));
p2=integral(f2,0, 0.1209);

f3=@(l) ((l./0.1209).^16).*exp(-347*(x-l)).*(bessely(0, 347.*((x-l).^2+z.^2).^0.5)+bessely(0, 347.*((x-l).^2+(0.1327-z).^2).^0.5));
p3=integral(f3,0 ,0.1209);
%fp3=vpaintegral(f3,l,[0 0.1209]);
y=14483.7*(0.27*p1+0.7*p2+0.7*p3);
end