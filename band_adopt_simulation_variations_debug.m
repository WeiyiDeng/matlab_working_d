% % simulate friend matrix
% NF = 100;                    % # of users in friend matrix
% PF = 0.05;                   % percentage of people who are friends in the friend matrix
% f_score = rand(NF,NF);
% f_friend = f_score>(1-PF);
% f_temp = triu(f_friend);
% f_temp2 = tril(f_friend');
% f_combi = f_temp+f_temp2;
% friend_matrix = f_combi.*(1-diag(ones(1,NF)));
% 
% mean(sum(friend_matrix))     % avg # of friends per user
% 
% % simulate similarity matrix
% PN = 0.05;                   % percentage of neighbours 
% s_score = 0.5+randn(NF,NF);
% 
% s_temp = triu(s_score);
% s_temp2 = tril(s_score');
% s_combi = s_temp+s_temp2;
% 
% neighbour_threshold = quantile(s_combi(:),1-PN);
% s_neighbour = s_combi>neighbour_threshold;
% neighbour_matrix = s_neighbour.*(1-diag(ones(1,NF)));
% 
% s_combi_std = (s_combi-min(s_combi(:)))/(max(s_combi(:))-min(s_combi(:)));      % standardize similarity scores to be in 0-1 range
% similarity_matrix = s_combi_std.*(1-diag(ones(1,NF)))+diag(ones(1,NF));
% similarity_matrix(similarity_matrix<0) = 0;
% 
% mean(sum(neighbour_matrix))     % avg # of neighbours per user
% mean(sum(neighbour_matrix.*friend_matrix))        % avg # of overlap between neighbour and friend per user
% 
% % simulate band attraction
clc
clear

tic
%% no overlaps among members and among friends
sth = [];
rep_times = 10;
% rep_times = 2;
% rep_times = 1;
NM_set = [];
NM_set = repmat(30,rep_times,1);
% NM_set = [NM_set; repmat(12,rep_times,1)];
% NM_set = [NM_set; repmat(20,rep_times,1)];
% NM_set = [NM_set; repmat(40,rep_times,1)];
% NM_set = [NM_set; repmat(60,rep_times,1)];
% NM_set = [NM_set; repmat(150,rep_times,1)];
% NM_set = [NM_set; repmat(200,rep_times,1)];

% CHUNKS = 50;
% CHUNKS = 5;
CHUNKS = 1;
% combined_set = NM_set*CHUNKS;

metrics_out_avg_set = zeros(length(NM_set),13);
average_member_adopts = zeros(length(NM_set),1);

beta_A_set = 1:10;
% beta_A_set = 1:2;
% beta_A_set = [1 5 10];
% beta_A_set = [1 1];

alpha_R = -12;
se_R = -1;
% beta_C = 1;
beta_C = 0;
se_C = 1;
beta_P1 = -0.3;
se_P1 = 2;
beta_P2 = -0.3;
se_P2 = 3;
% beta_A = 10;                               % social influence parameter

beta_A_set_disp = repmat(beta_A_set,rep_times,1);
beta_A_set_disp = beta_A_set_disp(:);

metrics_out_avg_set_SI_variant = zeros(length(beta_A_set_disp),13);
average_member_adopts_SI_variant = zeros(length(beta_A_set_disp),1);
a_ind = 1;
for a = 1:length(beta_A_set)
        beta_A = beta_A_set(a);
        rng(1)

        seed1 = 1;
        
    for d = 1:length(NM_set)
        rng(100+seed1)
        % generate members first 
        N_pool = 10000;                                            % # of users in the whole pool
        % NM = 100;                                              % number of members
        NM = NM_set(d);                                              % number of neighbours per member
        NN = 50;
        PF2 = 0.005;                                           % prob of becoming friends with another user in the pool
        % num_friends = binornd(N_pool,PF2,NM,1);               % simulate # of friends drawn from NM members

        r = 3;
        p = 0.95;
        num_friends = nbinrnd(r,1-p,1,NM);                      % simulate # of friends drawn from NM members (drawn from skewed negative binomial distribution)
        disp(sum(num_friends))
        
        % multinomial distribution
        % p = [0.2,0.3,0.5];
        % R = mnrnd(10000,p,1000);

        metrics_chunks_cell = cell(NM*CHUNKS,1);

        seed2 = 1;
        ind_c = 1:NM;

        for c = 1:CHUNKS
            rng(100000+seed2)
            %     debug = rand(1,1);                                       % debug
            %     disp(debug)

            % simulate friend matrix
            sizeF = sum(num_friends)+NM;
            % prep_friend_index = cumsum(num_friends)+100;
            temp = repelem(1:NM,num_friends);
            F_sparse = sparse(temp,NM+(1:sum(num_friends)),ones(sum(num_friends),1),sizeF+NM*NN,sizeF+NM*NN);      % F is the friend matrix for all users
            F_sparse_inv = F_sparse';
            F = full(F_sparse+F_sparse_inv);
            rowF = temp;
            colF = NM+(1:sum(num_friends));
            M_filter = ones(sizeF+NM*NN,sizeF+NM*NN);                      % create m_filter that filter out all members (set to be 0) 
            M_filter(1:NM,:) = 0;

            % simulate neighbour matrix
            rowN = repmat(1:NM,NN,1);
            rowN = rowN(:);
            colN = NM + sum(num_friends)+(1:NM*NN);
            N_sparse = sparse(rowN,colN,ones(size(rowN)),sizeF+NM*NN,sizeF+NM*NN);                               % N is the neighbour matrix for all users
            N_sparse_inv = N_sparse';
            N = full(N_sparse+N_sparse_inv);

            % simulate similarity matrix
            similarity_vec = zeros(sum(num_friends),1);
            ind = 0;
            for i = 1:NM
                similarity_kj = rand(num_friends(i),1);           % simulate similarity scores between each member and friends
                row_indice = (ind+1):(ind+num_friends(i));
                similarity_vec(row_indice) = similarity_kj;
                ind = ind+num_friends(i);
            end

            similarity_vec2 = zeros(NM*NN,1);
            ind = 0;
            for i = 1:NM
                similarity_kj = 0.9+0.1*rand(NN,1);               % simulate similarity scores between each member and neighbours (similarity within 0.9-1 range)
                row_indice = (ind+1):(ind+NN);
                similarity_vec2(row_indice) = similarity_kj;
                ind = ind+NN;
            end
%             disp(sum(sum(similarity_vec2)))

            S_sparse = sparse([rowF rowN'],[colF colN],[similarity_vec; similarity_vec2],sizeF+NM*NN,sizeF+NM*NN);
            S_sparse_inv = S_sparse';
            S = full(S_sparse+S_sparse_inv);
            
            % simulate extra similarity matrices to be used for temp shocks and band attractiveness 
            similarity_vec3 = zeros(sum(num_friends),1);
            ind = 0;
            for i = 1:NM
                similarity_kj = rand(num_friends(i),1);           % simulate similarity scores between each member and friends
                row_indice = (ind+1):(ind+num_friends(i));
                similarity_vec3(row_indice) = similarity_kj;
                ind = ind+num_friends(i);
            end

            similarity_vec4 = zeros(NM*NN,1);
            ind = 0;
            for i = 1:NM
                similarity_kj = 0.9+0.1*rand(NN,1);               % simulate similarity scores between each member and neighbours (similarity within 0.9-1 range)
                row_indice = (ind+1):(ind+NN);
                similarity_vec4(row_indice) = similarity_kj;
                ind = ind+NN;
            end
            
            similarity_vec5 = zeros(sum(num_friends),1);
            ind = 0;
            for i = 1:NM
                similarity_kj = rand(num_friends(i),1);           % simulate similarity scores between each member and friends
                row_indice = (ind+1):(ind+num_friends(i));
                similarity_vec5(row_indice) = similarity_kj;
                ind = ind+num_friends(i);
            end

            similarity_vec6 = zeros(NM*NN,1);
            ind = 0;
            for i = 1:NM
                similarity_kj = 0.9+0.1*rand(NN,1);               % simulate similarity scores between each member and neighbours (similarity within 0.9-1 range)
                row_indice = (ind+1):(ind+NN);
                similarity_vec6(row_indice) = similarity_kj;
                ind = ind+NN;
            end            
            
            % simulate temp shocks 
            T = 52;                                               % number of weeks
            NB = 100;
%             Cikt = beta_C+se_C*randn(NB,T,NM);                  % ??? w: not generating equal number of positive and negative numbers?
            % test with mean(mean(mean(Cikt)))
            Cikt = beta_C+se_C*(rand(1,1)*2-1).*randn(NB,T,NM);
%             disp(rand(1))
            disp(sum(sum(sum(Cikt))))                    % debug
            Cijt = zeros(NB,T,sum(num_friends));                  % simulate shocks for friends
            % Cijt = beta_C+se_C*randn(NB,T,sum(num_friends));
            ind = 0;
            for i = 1:NM
                row_indice = (ind+1):(ind+num_friends(i));
                Sjk = similarity_vec3(row_indice);
                temp = repmat(Sjk,1,NB,T);
                Sjk_ext = permute(temp,[2 3 1]);        
                Cikt_ext = repmat(Cikt(:,:,i),1,1,length(Sjk));
%                 Cijt_part = Sjk_ext.*Cikt_ext+(1-Sjk_ext).*(beta_C+se_C*randn(NB,T,length(Sjk)));
                Cijt_part = Sjk_ext.*Cikt_ext+(1-Sjk_ext).*(beta_C+se_C*(rand(1,1)*2-1).*randn(NB,T,length(Sjk)));
            %     row_indice = (ind+1):(ind+num_friends(i));
                Cijt(:,:,row_indice)=Cijt_part;
                ind = ind+num_friends(i);
            end
%             disp(sum(sum(sum(Cijt))))

            Cint = zeros(NB,T,NN*NM);                               % simulate shocks for neighbours
            % Cint = beta_C+se_C*randn(NB,T,NN*NM);
            ind = 0;
            for i = 1:NM
                row_indice = (ind+1):(ind+NN);
                Snk = similarity_vec4(row_indice);
                temp = repmat(Snk,1,NB,T);
                Snk_ext = permute(temp,[2 3 1]);        
                Cikt_ext = repmat(Cikt(:,:,i),1,1,length(Snk));
%                 Cint_part = Snk_ext.*Cikt_ext+(1-Snk_ext).*(beta_C+se_C*randn(NB,T,length(Snk)));
                Cint_part = Snk_ext.*Cikt_ext+(1-Snk_ext).*(beta_C+se_C*(rand(1,1)*2-1).*randn(NB,T,length(Snk)));
                Cint(:,:,row_indice)=Cint_part;
                ind = ind+NN;
            end
%             disp(sum(sum(sum(Cint))))

            C = cat(3,Cikt,Cijt,Cint);          % shock per band per t for all users
            C = permute(C,[3 1 2]);                                                                     % NU*NB*t

            % simulate band attraction/preference
            % band_attractiveness_mean = rand(1,NB);                        % can also add different means for different bands
            % Pik = -8*repmat(band_attractiveness_mean,NM,1)+3*randn(NM,NB);
            % Pik = -8+3*randn(NM,NB);
            bandPref = @(beta1,sigma1, beta2, sigma2, nrow, ncol) repmat(beta1+sigma1*randn(1,ncol),nrow,1)+ (beta2+sigma2*randn(nrow,ncol));            % this approach also gives different means for different band and different draws for different users as well
            Pik = bandPref(beta_P1,se_P1,beta_P2,se_P2,NM,NB);
%             disp(sum(sum(Pik)))                                      % debug
            ind = 0;
            Pij = zeros(sum(num_friends),NB);                  % simulate band preferences for friends
            for i = 1:NM
                row_indice = (ind+1):(ind+num_friends(i));
                Sjk = similarity_vec5(row_indice);
                Sjk_ext = repmat(Sjk,1,NB);
                Pik_ext = repmat(Pik(i,:),length(Sjk),1);
            %     Pij_part = Sjk_ext.*Pik_ext+(1-Sjk_ext).*(-8+3*randn(length(Sjk),NB));
            %     Pij_part = Sjk_ext.*Pik_ext+(1-Sjk_ext).*(-8*repmat(band_attractiveness_mean,length(Sjk),1)+3*randn(length(Sjk),NB));
                Pij_part = Sjk_ext.*Pik_ext+(1-Sjk_ext).*bandPref(beta_P1,se_P1,beta_P2,se_P2,length(Sjk),NB);
                Pij(row_indice,:) = Pij_part;
                ind = ind+num_friends(i);
            end
            % Pij = bandPref(beta_P1,se_P1,beta_P2,se_P2,sum(num_friends),NB);                % band preference without similarity between users

            ind = 0;
            Pin= zeros(NN*NM,NB);                  % simulate band preferences for neighbours
            for i = 1:NM
                row_indice = (ind+1):(ind+NN);
                Snk = similarity_vec6(row_indice);
                Snk_ext = repmat(Snk,1,NB);
                Pik_ext = repmat(Pik(i,:),length(Snk),1);
            %     Pin_part = Snk_ext.*Pik_ext+(1-Snk_ext).*(-8+3*randn(length(Snk),NB));
            %     Pin_part = Snk_ext.*Pik_ext+(1-Snk_ext).*(-8*repmat(band_attractiveness_mean,length(Snk),1)+3*randn(length(Snk),NB));
                Pin_part = Snk_ext.*Pik_ext+(1-Snk_ext).*bandPref(beta_P1,se_P1,beta_P2,se_P2,length(Snk),NB);
                Pin(row_indice,:) = Pin_part;
                ind = ind+NN;
            end
            % Pin = bandPref(beta_P1,se_P1,beta_P2,se_P2,NN*NM,NB);                         % band preference without similarity between users

            P = cat(1,Pik,Pij,Pin);          % band preference of all bands for all users                % NU*NB
%             disp(sum(sum(P)))
            
            % simulate baseline adoption tendency
            R = alpha_R+se_R*randn(NM+sum(num_friends)+NM*NN,NB,T);                                                   % NU*NB*t
%             disp(sum(sum(sum(R))))


            % simulate band adoption t=1
            % betaP = 0.5;
            Draw_threshold = rand(NM+sum(num_friends)+NM*NN,NB,T);                                       % NU*NB*t
            % Draw_threshold = rand(NM+sum(num_friends)+NM*NN,NB);                                           % NU*NB
            % Draw_threshold = repmat(Draw_threshold,1,1,T);                      % keep the same draws across different T                % NU*NB*t
%             disp(sum(sum(sum(Draw_threshold))))                                         % debug
            Prob1 = exp(R(:,:,1)+C(:,:,1)+P)./(1+exp(R(:,:,1)+C(:,:,1)+P));
            A1 = Prob1>Draw_threshold(:,:,1);
            % sth = R(:,:,1)+C(:,:,1)+P;
            % hist(sth(:))
            NU = NM+sum(num_friends)+NM*NN;
            sum(sum(A1))/(NU*NB)

            % simulate band adoption t=2
%             Prob2 = exp(R(:,:,2)+beta_A*F.*S*A1+C(:,:,2)+P)./...
%                 (1+exp(R(:,:,2)+beta_A*F.*S*A1+C(:,:,2)+P));
            Prob2 = exp(R(:,:,2)+beta_A*F.*M_filter.*S*A1+C(:,:,2)+P)./...         % use member filter to filter out all social influence to members' adoption decisions
                (1+exp(R(:,:,2)+beta_A*F.*M_filter.*S*A1+C(:,:,2)+P));
            A2 = Prob2>Draw_threshold(:,:,2);
            % A2 = A2+A1;
            % A2(A2>0)=1;
            A2(A1>0)=0;                                % A2 is the adoptions occured in t=2 only, not including those who have adopted in t=1
            sum(sum(A2))/(NU*NB)

            % simulate adoptions for all t
            TI = 4;                                    % friend adoption is influential within this time interval
            A = zeros(NU,NB,T);
            A(:,:,1) = A1;
            for t = 2:T
                if t <TI+1
%                     Prob_t = exp(R(:,:,t)+beta_A*F.*S*sum(A(:,:,1:t-1),3)+C(:,:,t)+P)./...
%                         (1+exp(R(:,:,t)+beta_A*F.*S*sum(A(:,:,1:t-1),3)+C(:,:,t)+P));
                    Prob_t = exp(R(:,:,t)+beta_A*F.*M_filter.*S*sum(A(:,:,1:t-1),3)+C(:,:,t)+P)./...                       % use member filter to filter out all social influence to members' adoption decisions
                        (1+exp(R(:,:,t)+beta_A*F.*M_filter.*S*sum(A(:,:,1:t-1),3)+C(:,:,t)+P));
                    A_t = Prob_t>Draw_threshold(:,:,t);
                    A_t(sum(A(:,:,1:t-1),3)>0)=0;
                    A(:,:,t) = A_t;
                else
%                     Prob_t = exp(R(:,:,t)+beta_A*F.*S*sum(A(:,:,t-TI:t-1),3)+C(:,:,t)+P)./...
%                         (1+exp(R(:,:,t)+beta_A*F.*S*sum(A(:,:,t-TI:t-1),3)+C(:,:,t)+P));
                    Prob_t = exp(R(:,:,t)+beta_A*F.*M_filter.*S*sum(A(:,:,t-TI:t-1),3)+C(:,:,t)+P)./...                    % use member filter to filter out all social influence to members' adoption decisions
                        (1+exp(R(:,:,t)+beta_A*F.*M_filter.*S*sum(A(:,:,t-TI:t-1),3)+C(:,:,t)+P)); 
                    A_t = Prob_t>Draw_threshold(:,:,t);
                    A_t(sum(A(:,:,1:t-1),3)>0)=0;
                    A(:,:,t) = A_t;
                end
            end

            %% compute measures
            adoption = cell(T,1);                           % adoption observations in cells
            adoption_mat = zeros(sum(sum(sum(A))),3);       % adoption observations in matrix
            ind = 0;
            for t = 1:T
                A_temp = A(:,:,t);
                [i,j] = find(A_temp);
                adoption{t} = [i j];                % USER_ID*BAND_ID
                row_range = (ind+1):(ind+length(i));
                adoption_mat(row_range,1)=i;        % USER_ID
                adoption_mat(row_range,2)=j;        % BAND_ID
                adoption_mat(row_range,3)=t;        % WEEK_ID
                ind = ind+length(i);
            end

            member_adoption = adoption_mat(adoption_mat(:,1) <=NM,:);
            member_adoption = sortrows(member_adoption);
            disp(['average member adopts ' num2str(size(member_adoption,1)/NM)])
            average_member_adopts(d) = size(member_adoption,1)/NM;

            member_cells = cell(NM,1);
            for i = 1:NM
                member_cells{i}=member_adoption(find(member_adoption(:,1)==i),:);
            end

            friend_end = cumsum(num_friends);
            friend_start = [1 friend_end(1:end-1)+1];
            friend_end = friend_end+NM;
            friend_start = friend_start+NM;
            neighbour_end = (1:NM).*NN+sum(num_friends)+NM;
            neighbour_start = [sum(num_friends)+NM+1 neighbour_end(1:end-1)+1];
            % sum(num_friends)+NM

            adoption_mat_sparse = sparse(adoption_mat(:,1),adoption_mat(:,2),adoption_mat(:,3),NU,NB);

            metrics_store_cell = cell(NM,1);
            for m = 1:  NM
                m_adoption = member_cells{m};
                m_metrics_store = zeros(size(m_adoption,1),6);
                for i = 1:size(m_adoption,1)
                    band_id = m_adoption(i,2);
                    week_id = m_adoption(i,3);
                    f_adoption = adoption_mat_sparse(friend_start(m):friend_end(m),band_id);
                    within4Weeks_f = sum(f_adoption>week_id & f_adoption-week_id<=TI);
                    random_week_id = datasample(1:T,1);
            %         random4Weeks_f = sum(f_adoption>random_week_id & f_adoption-random_week_id<=TI);
                    fixed_4Weeks_f = sum(f_adoption(:)>0)/52*4;
                    n_adoption = adoption_mat_sparse(neighbour_start(m):neighbour_end(m),band_id);
                    within4Weeks_n = sum(n_adoption>week_id & n_adoption-week_id<=TI);
            %         random4Weeks_n = sum(n_adoption>random_week_id & n_adoption-random_week_id<=TI);
                    fixed_4Weeks_n = sum(n_adoption(:)>0)/52*4;
                    m_metrics_store(i,1) = within4Weeks_f;
            %         m_metrics_store(i,2) = random4Weeks_f;
                    m_metrics_store(i,2) = fixed_4Weeks_f;
                    m_metrics_store(i,3) = within4Weeks_n;
            %         m_metrics_store(i,4) = random4Weeks_n;
                    m_metrics_store(i,4) = fixed_4Weeks_n;
                end
                m_metrics_store(:,5) = num_friends(m);
                m_metrics_store(:,6) = NN;
                metrics_store_cell{m} = m_metrics_store;
            end

            metrics_chunks_cell(ind_c) = metrics_store_cell;

            ind_c = ind_c+NM;

            disp(a)
            disp(d)
            disp(c)

            seed2 = seed2+1;
            seed1 = seed1 +1;

        end

        metrics_out_cell = cell(NM*CHUNKS,1);
        metrics_out_avg = zeros(NM*CHUNKS,13);
        for m = 1:(NM*CHUNKS)
            m_metrics = metrics_chunks_cell{m};
            m_out = zeros(size(m_metrics,1),6);
            if isempty(m_metrics)
                m_out = 0;
                metrics_out_avg(m,:) = 0;
            else
            m_out(:,1) = m_metrics(:,1)./m_metrics(:,5);        % Adoption by friends within 4 weeks after member adopts/# of friends
            m_out(:,2) = m_metrics(:,2)./m_metrics(:,5);        % Adoption by friends in random 4 weeks/# of friends
            m_out(:,3) = m_metrics(:,3)./m_metrics(:,6);        % Adoption by neighbours within 4 weeks after member adopts/# of friends
            m_out(:,4) = m_metrics(:,4)./m_metrics(:,6);        % Adoption by neighbours in random 4 weeks/# of friends
            m_out(:,5) = m_out(:,1)-m_out(:,2);
            m_out(:,6) = m_out(:,3)-m_out(:,4);
            m_out(:,7) = m_metrics(:,1)./m_metrics(:,3)./m_metrics(:,5).*m_metrics(:,6);
            metrics_out_avg(m,1:6) = mean(m_metrics,1);         % average across bands
            m_out(find(m_out==Inf))=NaN;
            metrics_out_avg(m,7:13) = nanmean(m_out,1);            % average across bands (remove NaNs)
            metrics_out_avg(isnan(metrics_out_avg))=0;
            metrics_out_cell{m} = m_out;
            end
        end

        % run Wilcoxon rank sum test (equivalent to Mann-Whitney U-test)
        [p,h,stats] = ranksum(metrics_out_avg(:,11),metrics_out_avg(:,12))

        metrics_out_avg_set(d,:) = mean(metrics_out_avg,1);                 % average across members

        metrics_out_avg_set_SI_variant(a_ind,:) = metrics_out_avg_set(d,:);
        average_member_adopts_SI_variant(a_ind) = average_member_adopts(d);
        a_ind = a_ind+1;

        sth = cat(3,sth,metrics_out_avg);

        clearvars -EXCEPT NM_set metrics_out_avg_set average_member_adopts d alpha_R se_R beta_C se_C beta_P1 se_P1 beta_P2 se_P2 beta_A CHUNKS metrics_out_avg_set_SI_variant a_ind beta_A_set_disp beta_A_set beta_A a average_member_adopts_SI_variant seed1 sth R1 A1 C1 P1 S1 F1 D1 MSC1
        
%         seed1 = seed1 +1;

    end

end

% figure
% subplot(3,2,1)
% scatter(NM_set*CHUNKS, metrics_out_avg_set(:,1))
% title('within4Weeks_f')
% subplot(3,2,2)
% scatter(NM_set*CHUNKS, metrics_out_avg_set(:,2))
% % title('random4Weeks_f')
% title('fixed4Weeks_f')
% subplot(3,2,3)
% scatter(NM_set*CHUNKS, metrics_out_avg_set(:,3))
% title('within4Weeks_n')
% subplot(3,2,4)
% scatter(NM_set*CHUNKS, metrics_out_avg_set(:,4))
% % title('random4Weeks_n')
% title('fixed4Weeks_n')
% subplot(3,2,5)
% scatter(NM_set*CHUNKS, metrics_out_avg_set(:,5))
% title('# of friends')
% subplot(3,2,6)
% scatter(NM_set*CHUNKS, metrics_out_avg_set(:,6))
% title('# of neighbours')
% 
% figure
% subplot(3,2,1)
% scatter(NM_set*CHUNKS, metrics_out_avg_set(:,7))
% title('A_f(TI)/N_f')
% subplot(3,2,2)
% scatter(NM_set*CHUNKS, metrics_out_avg_set(:,8))
% % title('A_f(R)/N_f')
% title('A_f(Fix)/N_f')
% subplot(3,2,3)
% scatter(NM_set*CHUNKS, metrics_out_avg_set(:,9))
% title('A_n(TI)/N_n')
% subplot(3,2,4)
% scatter(NM_set*CHUNKS, metrics_out_avg_set(:,10))
% % title('A_n(R)/N_n')
% title('A_n(Fix)/N_n')
% subplot(3,2,5)
% scatter(NM_set*CHUNKS, metrics_out_avg_set(:,11))
% % title('A_f(TI)/N_f-A_f(R)/N_f')
% title('A_f(TI)/N_f-A_f(Fix)/N_f')
% subplot(3,2,6)
% scatter(NM_set*CHUNKS, metrics_out_avg_set(:,12))
% % title('A_n(TI)/N_n-A_n(R)/N_n')
% title('A_n(TI)/N_n-A_n(Fix)/N_n')

% figure
% scatter(NM_set*CHUNKS, average_member_adopts)
% title('average member adoptions')

% figure
% scatter(NM_set*CHUNKS, metrics_out_avg_set(:,8))
% title('Fixed 4 weeks friend adoptions/# of friends')

figure
scatter(beta_A_set_disp, average_member_adopts_SI_variant)
title('average member adoptions')

% figure
% scatter(beta_A_set_disp, metrics_out_avg_set_SI_variant(:,8))
% title('Fixed 4 weeks friend adoptions/# of friends')

figure
scatter(beta_A_set_disp, metrics_out_avg_set_SI_variant(:,13))
title('after 4 weeks friend adoptions/after 4 weeks neighbour adoptions')

toc

% save test.mat                             % NANs for NM>100 ??? 

% clear
% load('test.mat')

% w: change se_C values [0.1 1 5]
% w: adoption always increases
% w: when utility is already very low, only (large) positive shocks will have an impact on P (see oneNote)

figure
y = metrics_out_avg_set_SI_variant(:,13);
plot(reshape(y,10,10)')

% % test
% w = mean(Pin,2);
% ww = reshape(w,50,30);
% mean(ww)
% ttest(ww(:,find(mean(ww)==min(mean(ww)))),ww(:,find(mean(ww)==max(mean(ww)))))

% % first find the max and min populations of column 13 of metrics_out_avg_set (eg. row 4 and row 7)
% sth_max = sth(:,:,[4 14 24 34 44 54 64 74 84 94]);
% sth_min = sth(:,:,[7 17 27 37 47 57 67 77 87 97]);
% f_temp = squeeze(sth_max(:,1,:));
% mean(f_temp)
% n_temp = squeeze(sth_max(:,3,:));
% mean(n_temp)
% f_temp2 = squeeze(sth_min(:,1,:));
% mean(f_temp2)
% n_temp2 = squeeze(sth_min(:,3,:));
% mean(n_temp2)
% plot(1:10,mean(f_temp),1:10,mean(f_temp2))
% plot(1:10,mean(f_temp)./mean(n_temp),1:10,mean(f_temp2)./mean(n_temp2))
% plot(1:10,mean(n_temp),1:10,mean(n_temp2))