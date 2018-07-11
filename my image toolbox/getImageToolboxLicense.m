sp = ones(2000, 2000);
nres = 1000;

for j = 1:1000
    try
        disp(num2str(j))
        imresize(sp, [nres,nres]);
        break
    catch err
        disp(err.identifier)
        pause(10)
        continue
    end
end
