## Comparando resultados obtidos com a dsad que usa o cubature e a dsad original
dsad(c(200,20000),0.1,0.00001)
dpoix(c(200,20000),0.1,0.00001)
dsad2(c(200,20000),0.1,exp,rate=0.00001)

debug(dsad2)
