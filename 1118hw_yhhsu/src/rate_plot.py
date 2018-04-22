import numpy as np
import matplotlib.pyplot as plt
L=np.array((3,4,5,6,7,8))
N=2**L+1
H=1./(N-1)
v_time=np.array((0.106,0.339,2.1354,12.374,65.050,373.788))
v_rate=np.array((0.803,0.777,0.384,0.265,0.197,0.154))

w_time=np.array((0.208,0.752,3.173,16.723,86.205,427.954))
w_rate=np.array((0.718,0.675,0.454,0.310,0.231,0.182))

#Plot
v, = plt.semilogy(L,v_time,"r+")
w, = plt.semilogy(L,w_time,"b^")
plt.legend([v, w],['V-circle $\gamma=1$','W-circle $\gamma=2$'],loc=4)
plt.axis([2,9,0,1000])
plt.title("Computing time versus lattice size L, (m1,m2)=(3,3)")
plt.xlabel("Lattice size")
plt.ylabel("Time (s)")
plt.grid(True)
plt.savefig("time_size.png")
#plt.text(5,-20,"time={:.3f}\nLinear Fitting\ny={:.5f} x + {:.5f}".format(elapsed,slope,const))
plt.show()

#Plot
# v, = plt.plot(L,v_rate,"r+")
# w, = plt.plot(L,w_rate,"b^")
# plt.legend([v, w],['V-circle $\gamma=1$','W-circle $\gamma=2$'],loc=1)
# plt.axis([2,9,0,1])
# plt.title("Convergence rate versus lattice size L, (m1,m2)=(3,3)")
# plt.xlabel("Lattice size")
# plt.ylabel("Convergence Rate")
# plt.grid(True)
# plt.savefig("rate_size.png")
# #plt.text(5,-20,"time={:.3f}\nLinear Fitting\ny={:.5f} x + {:.5f}".format(elapsed,slope,const))
# plt.show()