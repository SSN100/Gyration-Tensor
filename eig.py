import math
import numpy as np
from numpy import linalg as LA 
gyr_ten=np.genfromtxt("Gyration_tensor_2.dat", usecols=(0, 1, 2)) #Gyration_tensor_2.dat is provided by the fortran program. Reading into an numpy array.
gyr_ten=np.reshape(gyr_ten, (1250, 3, 3))
for i in range(1250):
    a, b=LA.eig(gyr_ten[i, :, :]) # a is an array of eigen values, b is an array of eigen vetors(from origin) corresponding to each eigenvalues.
    print("%8.5f" % math.sqrt(sum(a))) # Print Radius of Gyration
    #Rgx=math.sqrt(a(1)+a(2))
    #Rgy=math.sqrt(a(2)+a(0))
    #Rgz=math.sqrt(a(0)+a(1))
    
    #print(b)
