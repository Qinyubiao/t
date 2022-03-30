% workpiece shear heat
% chip shear
% chip combined
a=-0.26:0.065:0.26;b=0:0.06:0.12;
r=[];t=[];
for i=1:length(a)
    s=a(i);
    for j=1:length(b)
        v=b(j);
        r(i,j)=equ5(s,v);
    end
end
disp(r);

[c,h]=contour(b,a,r,10,"ShowText", "On");
h.LevelList=round(h.LevelList,0)  %rounds levels to 3rd decimal place
clabel(c,h)
set(gca,'XDir','reverse')
set(gca,'YDir','reverse')
axis equal


%chip fric
function y=equ5(x,z)
syms l;
%qfric=476.91,0~0.12;486.64(1.1-x),0.12~1.1 
q1=besselk(0, 44*((x-l).^2+z.^2).^0.5);
q2=besselk(0, 44*((x-l).^2+(0.26-z).^2).^0.5);
f1=@(l) 2051*(exp(-44*(x-l))).*(q1+q2);
f11=@(l) 2093*0.13*l.*(exp(-44*(x-l))).*(q1+q2);
%p1=int(f1,l,[0 0.1209]);
fp1=vpaintegral(f1,l,[0.13 0.26]);
fp11=vpaintegral(f11,l,[0 0.13]);
f2=@(l) 2051*((l./0.26).^0.24).*exp(-44*(x-l)).*(q1+q2);
f22=@(l) 2093*0.13*l.*((l./0.26).^0.24).*exp(-44*(x-l)).*(q1+q2);
%p2=int(f2,l,[0 0.1209]);
fp2=vpaintegral(f2,l,[0.13 0.26]);
fp22=vpaintegral(f22,l,[0 0.13]);
f3=@(l) 2051*((l./0.26).^15).*exp(-44*(x-l)).*(q1+q2);
f33=@(l) 2093*0.13*l.*((l./0.26).^15).*exp(-44*(x-l)).*(q1+q2);
%p3=int(f3,l,[0 0.1209]);
fp3=vpaintegral(f3,l,[0.13 0.26]);
fp33=vpaintegral(f33,l,[0 0.13]);

%chip_shear cos(38.1-(-6.5))=0.71,sin(38.1-(-6.5))=0.7
q3=besselk(0, 44.*((x-0.26+0.7*l).^2+(z-0.71*l).^2).^0.5);
q4=besselk(0, 44.*((x-0.26+0.7*l).^2+(0.24-z-0.71*l).^2).^0.5);
f4=@(l) (exp(-44.*(x-0.26+0.7*l))).*(q3+q4);
fp4=vpaintegral(f4,l,[0 0.16]);

y=2457*fp4+0.27*(fp1+fp11)+0.66*(fp2+fp22)+0.66*(fp3+fp33);
end
