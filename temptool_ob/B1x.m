% workpiece shear heat
% chip shear
% chip combined
a=0:0.1:1.2;b=zeros(1,13);
r=[];t=[];
for i=1:length(a)
    s=a(i);
    for j=1:length(b)
        v=b(j);
        r(i,j)=equ5(s,v);
        t(i,j)=integ1(s,v);
    end
end
disp(r);
disp(t);
plot(a,r,'-or')
hold on
plot(a,t,'-ob')
xlabel('li')
ylabel('temperature')
legend('chip','tool')


%chip
function y=equ5(x,z)
syms l;
%qfric=476.91,0~0.12;486.64(1.1-x),0.12~1.1 
q1=besselk(0, 44*((x-l).^2+z.^2).^0.5);
q2=besselk(0, 44*((x-l).^2+(1.1-z).^2).^0.5);
f1=@(l) 2051*(exp(-44*(x-l))).*(q1+q2);
f11=@(l) 2093*l.*(exp(-44*(x-l))).*(q1+q2);
%p1=int(f1,l,[0 0.1209]);
fp1=vpaintegral(f1,l,[0.98 1.1]);
fp11=vpaintegral(f11,l,[0 0.98]);
f2=@(l) 2051*((l./1.1).^0.24).*exp(-44*(x-l)).*(q1+q2);
f22=@(l) 2093*l.*((l./1.1).^0.24).*exp(-44*(x-l)).*(q1+q2);
%p2=int(f2,l,[0 0.1209]);
fp2=vpaintegral(f2,l,[0.98 1.1]);
fp22=vpaintegral(f22,l,[0 0.98]);
f3=@(l) 2051*((l./1.1).^15).*exp(-44*(x-l)).*(q1+q2);
f33=@(l) 2093*l.*((l./1.1).^15).*exp(-44*(x-l)).*(q1+q2);
%p3=int(f3,l,[0 0.1209]);
fp3=vpaintegral(f3,l,[0.98 1.1]);
fp33=vpaintegral(f33,l,[0 0.98]);

%chip_shear cos(38.1-(-6.5))=0.71,sin(38.1-(-6.5))=0.7
q3=besselk(0, 44.*((x-1.1+0.7*l).^2+(z-0.71*l).^2).^0.5);
q4=besselk(0, 44.*((x-1.1+0.7*l).^2+(0.92-z-0.71*l).^2).^0.5);
f4=@(l) (exp(-44.*(x-1.1+0.7*l))).*(q3+q4);
fp4=vpaintegral(f4,l,[0 0.65]);

y=0.27*(fp1+fp11)+0.66*(fp2+fp22)+0.66*(fp3+fp33)+2457*fp4;
end

%tool
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