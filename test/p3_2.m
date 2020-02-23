m = [21, 41, 81, 161, 321];
Nt = [50,100,200,400];

count = 1;
for t=1:length(Nt)    
    for i=1:length(m)
    fname = [num2str(count),'S2.mat'];
    [X,Y,Q] = solver(m(i),m(i),Nt(t),2);
    save(fname,'X','Y','Q');
    count = count+1;
    end
end


