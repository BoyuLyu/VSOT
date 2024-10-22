terminalWeights=[
16,0;
13,0;
0,20;
0,4;
];
pInPtr=[
1,2,10.9,4.5;
1,3,12.1,-1.5;
2,3,-1.5,9.3;
2,3,14.2,0;
1,4,0,7.5];
ooutx = solve_min_cut(terminalWeights, pInPtr);