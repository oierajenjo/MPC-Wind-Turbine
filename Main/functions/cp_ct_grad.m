function res = cp_ct_grad(la,be,dcl,lambdaVec,pitchVec)
[~, i_la] = min(abs(lambdaVec-abs(la)));
[~, i_be] = min(abs(pitchVec-be));
res = dcl(i_la,i_be);
end