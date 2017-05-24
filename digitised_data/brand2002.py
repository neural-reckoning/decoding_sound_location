'''
BDs from Brand et al. (2002)
'''

from pylab import *

x=array([5.05833E+2,3.87099E+2,
4.99615E+2,3.52481E+2,
5.95787E+2,2.77908E+2,
6.29316E+2,2.17063E+2,
6.19982E+2,1.65828E+2,
6.98524E+2,2.04732E+2,
7.01943E+2,1.92280E+2,
8.00870E+2,1.68931E+2,
8.51077E+2,8.59680E+1,
8.80999E+2,5.55696E+1,
9.81730E+2,1.74802E+2,
9.98754E+2,1.19463E+2,
1.10726E+3,1.23816E+2,
1.20600E+3,1.18461E+2,
1.20010E+3,5.33900E+1,
1.31522E+3,5.49867E+1]).reshape((16,2))

BF=x[:,0]
BD=x[:,1]
#BD=BD[BF>400]
#BF=BF[BF>400]

#BF=hstack((BF,BF))
#BD=hstack((BD,-BD))

subplot(211)
plot(BF,BD*BF*1e-6,'.')
ylabel('Best phase')
xlabel('Best frequency (Hz)')
subplot(212)
hist(BD*BF*1e-6)
ylabel('Number of cells')
xlabel('Best phase')
show()
