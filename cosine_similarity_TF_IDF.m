function sim_vec1_vec2 = cosine_similarity_TF_IDF(vec1,vec2,IDF_weighting_matrix)

temp_nume = vec1*IDF_weighting_matrix*vec2';
temp_deno1 = sqrt(vec1*IDF_weighting_matrix*vec1');
temp_deno2 = sqrt(vec2*IDF_weighting_matrix*vec2');
sim_vec1_vec2 = temp_nume/(temp_deno1*temp_deno2);

end