% TEST_MGS_REMARK_VEC 
% test mgs_remark_vec correctly simplifies triggers
%
% run like 
%  result = runtests('test_mgs_remark_vec') 
% see https://www.mathworks.com/help/matlab/matlab_prog/write-script-based-unit-tests.html
%
% 20230904WF  - used as example for matlab testing

%% test output
x = mgs_remark_vec(ones(5,1));
assert(size(x,1) == 5, 'size')

%% many
% cue=50,img=100,isi=150,mgs=200 + side(1-4) + sceen(10,20,30)
%NB. need 0+1 value b/c vec's non-zero min will be removed from values
ttls =     [0, 254 61 64 111 114 161 164 231 234]+1;
expected = [     1  2  2  -3   3   4   4  -5   5];
enums = mgs_remark_vec(ttls),
% NaN == NaN is never true in matlab. remove missing values
assert(all(rmmissing(enums) == expected), '1-5 left/right')

