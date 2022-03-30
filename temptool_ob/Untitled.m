a=-0.5:0.01:0.1;b=-0.6:0.01:0
for i=1:length(a)
    m=a(i)
    for j=1:length(b)
        n=b(j)
        z(i,j)=m+n
    end
end
contour(a,b,z)

