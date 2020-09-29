function savePlot(filename,NameValueArgs)
% 
% 
% Example:
%    filename = [cd,filesep,'CellDensity_MAM_CPM'];
%    savePlot(filename,'XYWH',[150,0,600,700]);
%
% @Yuting Chen
% yuting.chen17@imperial.ac.uk


arguments
    filename (1,:) char = [cd,filesep,'test_temp'];
    NameValueArgs.XYWH (1,4) double = [150,150,250,225];
    NameValueArgs.needreply (1,1) char = 'Y'
    NameValueArgs.onlyPng (1,1) logical = false
    NameValueArgs.Units (1,:) char = 'points';%'centimeters'
    NameValueArgs.wholepage (1,1) logical = false
    NameValueArgs.targetSize (1,:) char = ''% '1c',
end
if NameValueArgs.wholepage
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
else
    set(gcf,'units',NameValueArgs.Units,'position',NameValueArgs.XYWH);
    % set(gcf,'units','centimeters','position',NameValueArgs.XYWH);
end


if ~isempty(NameValueArgs.targetSize)
    % ref: https://www.elsevier.com/authors/author-schemas/artwork-and-media-instructions/artwork-sizing
    % Number of pixels versus resolution and print size, for bitmap images
    switch(NameValueArgs.targetSize)
        case 'ms'% minimal size
            set(gcf,'units','centimeters','position',[5,5,3,5]);
        case '1c'% single column
            set(gcf,'units','centimeters','position',[5,5,9,10]);
        case '1.5c'% 1.5 column
            set(gcf,'units','centimeters','position',[5,5,14,15]);
        case 'dc'% double column
            set(gcf,'units','centimeters','position',[5,5,19,20]);
        case '2/3c'% 2/3 column
            set(gcf,'units','centimeters','position',[5,5,6,7]);
        case '1/2c'
            set(gcf,'units','centimeters','position',[5,5,4.5,5]);
    end
end

if strcmpi(NameValueArgs.needreply,'Y')
    reply = input('Are you sure to save the new file? Y/N [Y]:','s');
else
    reply = 'Y';
end

if reply == 'Y'
    saveas(gcf,filename,'png')
    
    if ~(NameValueArgs.onlyPng)
        % saveas(gcf,filename,'eps')
        % saveas(gcf,filename,'pdf')
        saveas(gcf,filename,'fig')
    end
    
else
    fprintf('Not saved.\n');
end

end