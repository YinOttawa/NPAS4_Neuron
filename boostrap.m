%boostrap
d = bootstrp(1000,@mean,group1) - bootstrp(1000,@mean,group2);
figure;hist(d);

prctile(d,5)
prctile(d,1)
prctile(d,0.1)


