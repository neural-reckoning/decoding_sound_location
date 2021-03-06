'''
From Pecka et al. 2008, Suppl Fig 2A
'''
from pylab import *

x=array([1.95611E+2,8.09537E+2,
2.44568E+2,8.45101E+2,
2.34168E+2,1.68302E+3,
5.10560E+2,2.25036E+2,
5.15338E+2,4.78386E+2,
5.34191E+2,4.15720E+2,
6.15433E+2,3.86553E+2,
6.45106E+2,3.87379E+2,
6.46594E+2,1.50056E+2,
6.86005E+2,1.19445E+2,
7.08922E+2,3.73208E+2,
7.08160E+2,4.83958E+2,
7.44130E+2,2.31640E+2,
7.45388E+2,5.76039E+1,
7.71709E+2,1.21509E+2,
6.93383E+2,2.30402E+2,
7.77975E+2,5.01431E+2,
7.78692E+2,4.06502E+2,
7.88760E+2,2.95958E+2,
7.89245E+2,2.32672E+2,
8.87300E+2,2.98022E+2,
8.77455E+2,2.34529E+2,
9.86167E+2,3.31522E+2,
9.41384E+2,2.67410E+2,
9.75674E+2,2.20565E+2,
9.31654E+2,1.24811E+2,
9.77323E+2,4.65289E+1,
1.19294E+3,1.29145E+2,
1.05915E+3,2.53653E+2,
1.05834E+3,3.32760E+2,
1.03451E+3,2.53240E+2
]).reshape((31,2))

BF=x[:,0]
BD=x[:,1]
#BD=BD[BF>400]
#BF=BF[BF>400]

#BF=hstack((BF,BF))
#BD=hstack((BD,-BD))

subplot(211)
plot(BF,BD*BF*1e-6,'.')
#plot(BF,BD,'.')
ylabel('Best phase')
xlabel('Best frequency (Hz)')
subplot(212)
hist(BD*BF*1e-6)
ylabel('Number of cells')
xlabel('Best phase')
xlim(0,0.5)

#figure()
#hist(BD,20)

figure()
x=rand(1000)
y=sin(.5*pi*x)*150
hist(y)

show()
