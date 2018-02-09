function sim_vec1_vec2 = cosine_similarity_by_random_TF_IDF(vec1,vec2,IDF_weighting_matrix)

temp_nume1 = sqrt(vec1*IDF_weighting_matrix*vec1');
temp_nume2 = sqrt(vec2*IDF_weighting_matrix*vec2');
temp_deno = sum(sum(IDF_weighting_matrix));
sim_vec1_vec2 = (temp_nume1*temp_nume2)/temp_deno;

end