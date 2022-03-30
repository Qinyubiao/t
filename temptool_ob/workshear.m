% workpiece shear heat
% rubbing heat
% workpiece combined
%a=-0.5:0.01:0.1;b=-0.6:0.01:0;
a=-0.5:0.02:0.5;b=-0.2:0.02:0.7;
r=[];
for i=1:length(a)
    s=a(i);
    for j=1:length(b)
        v=b(j);
        r(i,j)=equ5(s,v);
    end
end
disp(r)
[c,h]=contour(a,b,r',10,"ShowText", "On");
h.LevelList=round(h.LevelList,0)  %rounds levels to 3rd decimal place
clabel(c,h)
set(gca,'XDir','reverse')
axis equal


function y=equ5(x,z)
syms l;
q1=besselk(0, 61.*((x-0.79*l).^2+(z-0.62*l).^2).^0.5);
q2=besselk(0, 61.*((x-0.79*l).^2+(0.8-z-0.62*l).^2).^0.5);
f1=@(l) (exp(-61.*(x-0.79*l))).*(q1+q2);
fp1=vpaintegral(f1,l,[0 0.65]);

q3=besselk(0, 61*((x+l).^2+(z)^2).^0.5 );
f2=@(l) (exp(-61*(x+l))).*q3;
fp2=vpaintegral(f2,l,[0,0.08]);

y=2457*fp1+0.915*3173*2*fp2;
end