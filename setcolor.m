h = gcf;
lines = h.Children(1).Children;

c = [         0         0  212.5000
         0   42.5000  255.0000
         0  127.5000  255.0000
         0  212.5000  255.0000
   42.5000  255.0000  212.5000
  127.5000  255.0000  127.5000
  212.5000  255.0000   42.5000
  255.0000  212.5000         0
  255.0000  127.5000         0
  255.0000   42.5000         0
  212.5000         0         0
  127.5000         0         0];

indices = [11,10,9,8,7,6,5,1];
%indices = fliplr(indices);
for idx = 1:8
    lines(idx).Color = c(indices(idx),:)/255;
end