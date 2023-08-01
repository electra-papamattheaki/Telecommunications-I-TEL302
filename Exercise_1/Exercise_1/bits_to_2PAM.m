function [X] = bits_to_2PAM(b)

len = length(b);

for i=1:len
    if (b(i) == 0)
        X(i) = 1 ;
    else
        X(i)=-1;
    end   
end