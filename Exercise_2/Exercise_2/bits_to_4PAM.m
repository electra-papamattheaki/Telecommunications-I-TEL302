function [X] = bits_to_4PAM(b)

len = length(b);

    for i=1:len
        if (i == len)
            if (b(i) == 0)
            X(i) = 3;
            break;
            else 
            X(i) = -3;
            break;
            end
        end
        if (b(i) == 0 && b(i+1) == 0)
            X(i) = 3 ;
        elseif (b(i) == 0 && b(i+1) == 1)
            X(i) = 1 ;
        elseif (b(i) == 1 && b(i+1) == 1)
            X(i) = -1 ;
        else 
            X(i) = -3 ;
        end
    end   
end