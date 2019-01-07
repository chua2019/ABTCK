function r=isPerfectSubset(A,B)
scanIt2 = @(A,B) bsxfun(@eq,A,B);

scanIt1 = @(A,B) arrayfun(@(i) ...
   scanIt2(A,B(i:length(A)+i-1)),...
   1:length(B)-length(A),'UniformOutput',false);

r = any(cellfun(@all,scanIt1(A,B)));
end