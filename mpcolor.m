function h = mpcolor(x,y,C)
% function h = mpcolor(x,y,C)
% Pads a row of nan's at the bottom and a column on the right side of input
% C so that the resulting image matches the output of imagesc
% MS 2017-08-13
% Inputs:
% x, y are vectors - they will be expanded to the requirement of pcolor
% function
% C - 2d matrix

[r,c] = size(C);
cnew = [C nan(r,1); nan(1,c+1)];
pcolor(cnew);
set(gca,'xtick',x+0.5,'xticklabel',x,'ytick',y+0.5,'yticklabel',y)
h = gca;