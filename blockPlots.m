function [handle,tag] = blockPlots(X,Y,block,blockType,varargin)
% BLOCKPLOTS plot the Y in different line type according to Block
%
% Input: X: vector<double> [m,1]
%        Y: vector<double> [m,1]
%        block: cell<string> [m,1]
%        blockType: cell<string> [t,1] (t<m) (blockName <- unique(Block))
%        varargin: lineType: cell<string> [t,1]
%                  {xlabel,title}
% Output:
% 
% Example:
%        blockPlots(X,Y,Block,BlockName)
%        or 
%        blockPlots(X,Y,Block,BlockName,lineType);
%
% by Yuting CHEN
handle = setFigureProperty;
tag = 0;
% handle = figure;
if isempty(varargin)
    lineType = {'k:.','r-','b:.','r.:','rp--'};
else
    lineType = varargin{1};
    
end

for i = 1:length(blockType)
    ind = strcmp(block, blockType{i});
    plot(X(ind),Y(ind),lineType{i});
    hold on;
end
if ~isempty(varargin) && length(varargin)>1
    xlabel(varargin{2}{1});
    ylabel(varargin{2}{2});
end
tag = 1;
end

