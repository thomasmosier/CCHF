function prec = precision(number)
%Determine number of non-zero decimal places in number.

if isempty(number) || number == round(number)
    prec = 0;
else
    precTemp = diff(round2(number,(0:20)));
    prec = find(precTemp ~= 0, 1, 'last') + 1;
end
    
    