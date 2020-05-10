AAA=[-2,1,-3; -1,0,-2; 0,-1,-2]; 
% BBB=[2,2,-2; -1 5 -2; -1 1 2];
BBB=[2,2; -1 5; -1 1];
CCC=[5 -4; 1 0; 1 -1]; 
X=are(AAA,BBB,CCC), norm(AAA'*X+X*AAA-X*BBB*X+CCC)