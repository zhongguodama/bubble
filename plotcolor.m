function [ ] = plotcolor(x, data, tag_array )
%PLOTCOLOR Summary of this function goes here
%   Detailed explanation goes here
% @input {data}: plot each column as a line
% @input {tag_array}: cell array of tag string to show on the colormap
data = data';
a = tag_array;
tag = cellfun(@(x)sprintf(num2str(x)),a,'uniformoutput',false);
colors = jet;
[~, col] = size(data);
numberOfLines = col;
xcolor = [1:64]'/64;
xcolormap = (1:numberOfLines)'/numberOfLines;
cmap = interp1(xcolor,colors,xcolormap);
axes('ColorOrder',cmap,'NextPlot','replacechildren')
plot(x,data)
colormap(cmap)
delta = 1/(numberOfLines );
colorbar('Ticks',(delta:delta:1)-delta/2,...
         'TickLabels',tag)
end