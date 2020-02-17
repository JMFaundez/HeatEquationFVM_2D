m = [21, 51, 101];
Nt = [50,100,200];

count = 1;
for t=1:length(Nt)    
    for i=1:length(m)
    fname = [num2str(count),'S1.mat'];
    [X,Y,Q] = solver(m(i),m(i),Nt(t),1);
    save(fname,'X','Y','Q');
    count = count+1;
    end
end


