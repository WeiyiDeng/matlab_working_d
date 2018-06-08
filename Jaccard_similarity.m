function sim_vec1_vec2 = Jaccard_similarity(vec1,vec2)
% w: vec1 and vec2 should have the same length and should both be composed of binary values (eg. 0 1) only!

sim_vec1_vec2 = sum(vec1 & vec2)/sum(vec1 | vec2);

end