func = @(x,y)(1);

% http://people.scs.fsu.edu/~burkardt/m_src/stroud/stroud.html

%  Parameters:
%
%    Input, integer RULE, the rule desired.
%      1, 1 point 1-st degree;
%      2, 4 point 3-rd degree, Stroud S2:3-1;
%      3, 4 point 3-rd degree, Lether #1;
%      4, 4 point 3-rd degree, Stroud S2:3-2;
%      5, 5 point 3-rd degree;
%      6, 7 point 5-th degree;
%      7, 9 point 5-th degree;
%      8, 9 point 5-th degree, Lether #2;
%      9, 12 point 7-th degree;
%     10, 16 point 7-th degree, Lether #3;
%     11, 21 point 9-th degree, Stroud S2:9-3;
%     12, 25 point 9-th degree, Lether #4 (after correcting error);
%     13, 64 point 15-th degree Gauss product rule.
rule = 10;

[ xtab, ytab, weight ] = circle_xy_set ( rule, 0 );
order = circle_xy_size ( rule);
result = circle_xy_sum ( func, 0, 0, 1, order, xtab, ytab, weight )