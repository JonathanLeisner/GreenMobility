    
# write an array
a = ['first','second','third']

print('with index i the loop prints in the correct order')
for i in range(3):
    print(i,a[i])
    
print('with index i-1 the loop starts from the last element')
for i in range(3):
    print(i,a[i-1])
    
print('')


# now build a shock vector of size 10*4
# reshaping in a 10*4 two dimensional array
# the first index will be from 0 to 9
# the second will be from 0 to 3
import numpy as np
ush = np.random.normal(0,1,size=10*4).reshape(10,4)

try: 
    print(ush[9,4])
except:
    print('without try statement gives an out of bounds error')

print('')
# no error in what follows
print ('no error ',ush[9,3])