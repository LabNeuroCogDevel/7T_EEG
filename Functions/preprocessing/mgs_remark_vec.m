function [simple] = mgs_remark_vec(mark)
% REMARK_VEC simplify  biosemi LPT TTL triggers from mgs into events 1 to 5
% use sign to indicate side: negative is left, positive is right (?)
% 1. iti, 2. cue, 3. img, 4. dly, 5. mgs
% 20230904WF - could be used by preprocessing/remark.m; currently just to demonstrate testing

     mark = mark - min(mark(mark>0));
     mark(mark>65000) = 0;
     % isi (150+x) and iti (254) are different
     %    event inc in 50: (50-200: cue=50,img=100,isi=150,mgs=200)
     %    category inc in 10 (10->30: None,Outdoor,Indoor)
     %    side inc in 1 (1->4: Left -> Right)
     %        61 == cue:None,Left
     %        234 == mgs:Indoor,Right
     %1 254 = ITI
     %2 50<cue<100 [50+(c 10,20,30)+(s,1-4)]
     %3 100<img.dot<150 [100+(c 10,20,30)+(s,1-4)] +/-
     %4 150<delay<200 [150+(c 10,20,30)]
     %5 200<mgs<250 [200+(c 10,20,30)+(s,1-4)]

     ending = mod(mark',10);


     simple = nan(size(mark));
     simple(mark == 254)= 1;
     simple(mark>=50 & mark<100)= 2;
     simple(mark>=100 & mark<150 & ((ending == 1) + (ending == 2))')= -3;
     simple(mark>=100 & mark<150 & ((ending == 3) + (ending == 4))')= 3;

     simple(mark>=150 & mark<200)= 4;
     simple(mark>=200 & mark<250 & ((ending == 1) + (ending == 2))')= -5;
     simple(mark>=200 & mark<250 & ((ending == 3) + (ending == 4))')= 5;
end
