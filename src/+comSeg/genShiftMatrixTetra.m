function [xxshift, yyshift, zzshift] = genShiftMatrixTetra
xxshift = [];
yyshift = [];
zzshift = [];
for i = -1:1
    for j = -1:1
        for k = -1:1
            if((i==0 && j==0 && k~=0))
               for m = [-1,1]
                   for n = [-1,1]
                      xxshift = [xxshift, [0, 0,   m*1,   0, 0,   0, m*1, m*1]]; 
                      yyshift = [yyshift, [0, 0,   n*1, n*1, 0,   0, n*1, 0]];
                      zzshift = [zzshift, [0, k*2, 1*k, 1*k ,0, k*2, 1*k, 1*k]];
                   end
               end
            elseif(i==0 && j~=0 && k == 0)
               for m = [-1,1]
                   for n = [-1,1]
                      xxshift = [xxshift, [0,   0, m, 0, 0, 0,   m, m]]; 
                      yyshift = [yyshift, [0, j*2, j, j, 0, 2*j, j, j]];
                      zzshift = [zzshift, [0,   0, n, n ,0, 0,   n, 0]];
                   end
               end                
            elseif(i~=0 && j == 0 && k == 0)
               for m = [-1,1]
                   for n = [-1,1]
                      xxshift = [xxshift, [0, 2*i, i, i, 0, 2*i, i, i]]; 
                      yyshift = [yyshift, [0, 0,   m, 0, 0,   0, m, m]];
                      zzshift = [zzshift, [0, 0,   n, n, 0,   0, n, 0]];
                   end
               end            
            end
        end
    end
end


end