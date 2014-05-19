function structure = convert_bps_to_structure( bps, nres )

% structure = repmat('.',1,nres);

% for pseudoknot (pk) annotation
% pk pairs are in the end of bps
% bps_right_previous = 0;
% pk_flag = 0;
% 
% for i = 1:size( bps, 1 )
%     structure( bps(i,1) ) = '(';
%     structure( bps(i,2) ) = ')';
%     if bps_right_previous>bps(i,2) || bps_left_previous>bps(i,1)|| pk_flag;
%         structure( bps(i,1) ) = '[';
%         structure( bps(i,2) ) = ']';
%         pk_flag = 1;
%     end;
%     
%     bps_right_previous = bps(i,2);
% end

bps_all = zeros(nres,2);
for i = 1:nres
    bps_all(i,2) = i;
    if find(bps(:,2) == i) ~= 0;
        bps_all(i,1) = bps(find(bps(:,2) == i),1);
    end;
    if find(bps(:,1) == i) ~= 0;
        bps_all(i,1) = bps(find(bps(:,1) == i),2);
    end;
end;

fid = fopen('bps_temp.txt','w');
fprintf(fid,'  %d  ENERGY = 0  temp\n',nres);
for i = 1:nres
    fprintf(fid,'    %d A       %d    %d    %d    %d\n',i,i-1,i+1,bps_all(i,1),bps_all(i,2));
end;
fclose(fid);

structure = load_ct_file( 'bps_temp.txt' );
system('rm bps_temp.txt');
    