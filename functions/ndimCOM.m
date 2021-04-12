function [OUT] = ndimCOM(IN,type)    
    if strcmp(type,'manual')        
        imagesc(log10(IN));axis image
        h = impoly;
        hMask = createMask(h);
        IN = IN.*hMask;
    end
    C = cellfun(@(n) 1:n, num2cell(size(IN)),'uniformoutput',0);
    [C{:}] = ndgrid(C{:});
    C = cellfun(@(x) x(:), C,'uniformoutput',0);
    C = [C{:}];
    OUT(:,:) = IN(:).'*C/sum(IN(:),'double');
end