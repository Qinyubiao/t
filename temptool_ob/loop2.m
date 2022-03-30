%　aはｚ軸の座標，ｂはｘ軸の座標，ｔ(i,j)はx-z平面つまり二次元切削の工具側温度数値
a=0:0.05:1.2; b=0:0.02:0.48;
m=1;t=[];
for i=1:length(a)
    s=a(i);
    for j=1:length(b)
      v=b(j);
      t(i,j)=equ5(s,v);
    end
end
disp(t);
contour(a,b,t,'ShowText','on')

function y=equ5(x,z)
% f=37 degree, a=6 degress, sin(f-a)=0.51,cos(f-a)=0.86, xi=1.2-0.51li, zi=0.86li,
% Vx=779.7mm/s, e=14mm2/s.
f1=@(l) exp(-27.9.*(x-1.2+0.51*l)).*(bessely(0, 27.9.*((x-1.2+0.51*l).^2+(z-0.86*l).^2).^0.5)+bessely(0, 27.9.*((x-1.2+0.51*l).^2+(0.56-z-0.86*l).^2).^0.5));
y=integral(f1,0,0.5);
end