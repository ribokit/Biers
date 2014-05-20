function [ D_comb ] = combine_diagonal( D_diag, DEL, thres, shift)

% INPUT
%       
% Clarence Cheng, May 2014
%
% % Examples
% D_comb = combine_diag(D_filt, 1, 0.025, 2);
% D_comb_comb = combine_diag(D_comb, 0, 0.025, 2);
%

D_comb = D_diag*0;


%% Diagonal comparison

for i = 1:size(D_diag, 1)-10
    for j = 1:size(D_diag, 2)-10
        if D_diag(i,j) > 0
            tot = D_diag(i,j);
            for n = 3:6
                if D_diag(i+n,j+n) > D_diag(i,j)*thres
                    tot = [tot D_diag(i+n,j+n)];
                    size(tot);
                end
            end
            if length(tot) > 1
                D_comb(i+shift,j+shift) = mean(tot);
            else
                if DEL == 1
                else
                    D_comb(i,j) = D_diag(i,j);
                end
            end
        end
    end
end

    %mohcaplot( D_comb/3, seqpos, ligpos, 'Combine diagonals', '', '', secstr );
end


%% Vertical comparison 
% for h = 1:length(thres)
%     for i = 1:size(D_diag, 1)-10
%         for j = 1:size(D_diag, 2)-10
%             if j > i
%                 if D_diag(i,j) > 0
%                     tot = D_diag(i,j);
%                     for n = 8:10
%                         if D_diag(i,j+n) > D_diag(i,j)*thres(h)
%                             tot = [tot D_diag(i,j+n)];
%                             size(tot);
%                         end
%                     end
%                     if length(tot) > 1
%                         D_comb(i,j+shift) = mean(tot);
%                     else
%                         if DEL == 1
%                             D_comb(i,j+shift) = 0;
%                         elseif DEL == 0
%                             D_comb(i,j+shift) = D_diag(i,j+shift);
%                         end
%                     end
%                 end
%             elseif i > j
%                 if D_diag(i,j) > 0
%                     tot = D_diag(i,j);
%                     for n = 8:10
%                         if D_diag(i+n,j) > D_diag(i,j)*thres(h)
%                             tot = [tot D_diag(i+n,j)];
%                             size(tot);
%                         end
%                     end
%                     if length(tot) > 1
%                         D_comb(i+shift,j) = mean(tot);
%                     else
%                         if DEL == 1
%                             D_comb(i+shift,j) = 0;
%                         elseif DEL == 0
%                             D_comb(i+shift,j) = D_diag(i+shift,j);
%                         end
%                     end
%                 end
%             end
%         end
%     end
% 
%     mohcaplot( D_comb/3, seqpos, ligpos, 'Combine diagonals', '', '', secstr );
% end
