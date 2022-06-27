f = 100;
d = 100;

refract = [1,0,;-1/f,0];
distance = [1,d;0,1];

osys = refract*distance;
[v,d] = eig(osys)
[U,S,V] = svd(distance)