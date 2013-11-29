#Problem 1

#The functions
#Bisection Method
bisection = function(fn, interval, tol, max.itr, option)
{
  for(i in 1:max.itr){
    c = mean(interval)    #(l+u)/2
    g.c = eval({l = c; fn})           #calculate g(c)
    if(abs(g.c) < tol) break    #check convergence
    if( eval({l = interval[1]; fn})*g.c < 0 ){
      interval[2] = c
    }else{interval[1] = c}
    if(option == T & i %% 10000 == 0) print(paste('Finishing the ',i,'th iteration.'))
  }
  print(paste('Number of iterations =',i))
  return(c)
}

#Newton Raphson
NR = function(fn, fn.der, start, tol, max.itr, option)
{
  for(i in 1:max.itr){
    new = start - eval({l = start; fn})/eval({l = start; fn.der})
    if(abs(eval({l = new; fn})) < tol){ 
      break
    }else{start = new}   #update
       if(option == T & i %% 10000 == 0) print(paste('Finishing the ',i,'th iteration.'))
  }
  return(new)
}

f = expression(log( (2+l)^125*(1-l)^38*l^34) )
f = expression( 125*log(2+l) + 38*log(1-l) + 34*log(l) )
D.f = D(f,'l')
D.f2 = D(D.f, 'l')



#Analysis
#Bisection
bisection(D.f, c(0.01,0.99), 1e-6, 100000, T)

#NR
tmp = seq(0.01, 0.99, by = 0.01)
NR.res = sapply(tmp, function(i) NR(D.f, D.f2, i, 1e-6, 10000,T))
plot(tmp, NR.res, xlab = 'Initial Values', ylab = 'NR solutions')