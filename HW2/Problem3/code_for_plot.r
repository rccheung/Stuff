tw = read.csv('C:/Users/user/School Work/Fall 2013/STA 250/HW2/Problem3/means/000000_0.csv',header = F)

means = numeric(1000)
for(i in 1:1000){
	a = as.character(tw[i,])
	means[i] = as.numeric(unlist(strsplit(a,"\001"))[2])
}

tw = read.csv('C:/Users/user/School Work/Fall 2013/STA 250/HW2/Problem3/vars/000000_0.csv',header = F)

vars = numeric(1000)
for(i in 1:1000){
	a = as.character(tw[i,])
	vars[i] = as.numeric(unlist(strsplit(a,"\001"))[2])
}

plot(means,vars, main = 'mean variance plot for problem 3')