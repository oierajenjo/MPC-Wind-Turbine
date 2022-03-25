function [la,res] = cp_max(be,cl,lambdaVec,pitchVec)
[~,i_be] = min(abs(pitchVec-be));
l_c = cl(:,i_be);
[res,i_la] = max(l_c);
la = lambdaVec(i_la);
end
