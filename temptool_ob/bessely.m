syms x z l w
y=exp(-w-(347*((x-0.8*l)^2+z^2)^0.5)^2/(4*w))/w
F=int(y,w,0 ,inf)