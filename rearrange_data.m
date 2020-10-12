
copyfile *RVKC437* RVKC437 


A=dir('*RVKC437*')
i=1;
B=[];
for i=1:length(A)
    date1 = A(i).name(3:10);
    datenew=strrep(date1,'-','');
    
    B(i) = string(datenew);
end
B=B';
BB=string(B);

for i=1:length(A)
% mkdir(BB(i))
cd(A(i).name)
C = dir('*rhd');

E = append('I:\BackUp_Data\OPTO_data\RVKC437\',BB(i));
if length(C)~=0; 
copyfile('*rhd', E)
end
cd ..

end

A=dir('*19*');
for i=1:length(A)
    date1(i) = string(A(i).name);
end


date1=date1';
date1=str2double(date1)
B=B';
BB=string(B);
cd log
a = dir('*.xlsx');
b=readtable(a(1).name);
for i=1:length(date1)
b(i,1)=date1(i);
end
writetable(b,a(1).name) %write back to excel in log folder

