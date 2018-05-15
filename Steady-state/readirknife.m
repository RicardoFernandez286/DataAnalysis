data=[];
directory=pwd;
files = dir(directory);
files(1:2,:) = [];
names = transpose({files.name});
start=6.2;
finish=6.56;
deltax = 0.005;
n=(finish-start)/deltax;
for i = 1:length(names)
    name = char(names(i));
    d = csvread(name);
    data = [data d(:,1)];
end
y = sum(data,1);
x = 0:deltax:n*deltax;
x=transpose(x); y=transpose(y);
dy = diff(y)/deltax;
ysort = sort(y);
dysort = diff(ysort)/deltax;
dx=x(1:end-1);
plot(x,ysort); figure;
plot(dx,dysort);