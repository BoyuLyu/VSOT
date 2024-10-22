%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%FUNCTION Indices = ISMEMBERROWS(A,B,Flag)
%%
%%      - When A and B are matrices with the same number of columns, then IsMemberRows
%%        returns a logical array containing true (1) for every row of A which occurs
%%        in B. False (0) otherwise.
%%      - A and B can be character matrices.
%%      - If a Flag is specified as a third argument, such as '', then all possible
%%        permutations of column indices are considered.
%%      - Gives same results as ismember(A,B,'rows'), except provides the answer much,
%%        much quicker. Approximately 20310 times faster (on my computer) when A and B  
%%        are both 3000x10 matrices.
%%
%%      Example:  >> A = [ 1   2   3 
%%                         4   5   6 
%%                         7   8   9 ]
%%
%%                >> B = [ 4   5   6 
%%                         10  11  12]
%%
%%                >> I = ismemberrows(A,B)
%%
%%                   I = [ 0
%%                         1
%%                         0 ]
%%
%%
%%
%%      Example:  >> A = [ 1   2   3 
%%                         4   5   6 
%%                         7   8   9 ]
%%
%%                >> B = [ 5   4   6 
%%                         10  11  12]
%%
%%                >> I = ismemberrows(A,B,'')
%%
%%                   I = [ 0
%%                         1
%%                         0 ]
%%
%%      - See also UNIQUE, INTERSECT, SETDIFF, SETXOR, UNION.
%%
%%   Tested under MATLAB 5.3.1
%%
%%   Steven Holden
%%   sholden@engr.mun.ca
%%   16/03/2000
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
