import timeit
from multiprocessing import Pool
import time
import math

fpool = False
N = 5000000

t1 = timeit.default_timer() 
def cube(x):
    return math.sqrt(x)

if fpool:
  if __name__ == "__main__":
      with Pool() as pool:
        result = pool.map(cube, range(10,N))
      print("Pool Program finished!")
else:
  for x in range(10,N):
    result = cube(x)
  print("Serial Program finished!")

t2 = timeit.default_timer()
print(f"Process time: {t2-t1}")

