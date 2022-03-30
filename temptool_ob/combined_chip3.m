% workpiece shear heat
% chip shear
% chip combined
a=0:0.1:1.1;b=0:0.115:0.46;
r=[];t=[];
for i=1:length(a)
    s=a(i);
    for j=1:length(b)
        v=b(j);
        r(i,j)=equ5(s,v);
    end
end
disp(r);

[c,h]=contour(b,a,r,20,"ShowText", "On");
h.LevelList=round(h.LevelList,0)  %rounds levels to 3rd decimal place
clabel(c,h)
set(gca,'XDir','reverse')
set(gca,'YDir','reverse')
axis equal


%chip fric
function y=equ5(x,z)
syms l;
%qfric=476.91,0~0.12;486.64(1.1-x),0.12~1.1 
q1=besselk(0, 44*0.238*((x-l).^2+z.^2).^0.5);
q2=besselk(0, 44*0.238*((x-l).^2+(1.1-z).^2).^0.5);
f1=@(l) 2051*0.238*(exp(-44*0.238*(x-l))).*(q1+q2);
f11=@(l) 2093*0.238*l.*(exp(-44*0.238*(x-l))).*(q1+q2);
%p1=int(f1,l,[0 0.1209]);
fp1=vpaintegral(f1,l,[0.98 1.1]);
fp11=vpaintegral(f11,l,[0 0.98]);
f2=@(l) 2051*0.238*((l./1.1).^0.24).*exp(-44*0.238*(x-l)).*(q1+q2);
f22=@(l) 2093*0.238*l.*((l./1.1).^0.24).*exp(-44*0.238*(x-l)).*(q1+q2);
%p2=int(f2,l,[0 0.1209]);
fp2=vpaintegral(f2,l,[0.98 1.1]);
fp22=vpaintegral(f22,l,[0 0.98]);
f3=@(l) 2051*0.238*((l./1.1).^15).*exp(-44*0.238*(x-l)).*(q1+q2);
f33=@(l) 2093*0.238*l.*((l./1.1).^15).*exp(-44*0.238*(x-l)).*(q1+q2);
%p3=int(f3,l,[0 0.1209]);
fp3=vpaintegral(f3,l,[0.98 1.1]);
fp33=vpaintegral(f33,l,[0 0.98]);

%chip_shear cos(38.1-(-6.5))=0.71,sin(38.1-(-6.5))=0.7
q3=besselk(0, 44.*0.238*((x-1.1+0.7*l).^2+(z-0.71*l).^2).^0.5);
q4=besselk(0, 44.*0.238*((x-1.1+0.7*l).^2+(0.92-z-0.71*l).^2).^0.5);
f4=@(l) (exp(-44.*0.238*(x-1.1+0.7*l))).*(q3+q4);
fp4=vpaintegral(f4,l,[0 0.65]);

y=585*fp4+0.13*(fp1+fp11)+0.58*(fp2+fp22)+1.5*(fp3+fp33);
end

