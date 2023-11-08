#
# Reference: Maziar Raissi et.al
# "Physics Informed Deep Learning (Part I): Data-driven Solutions of Nonlinear Partial Differential Equations" 2017
#
import sys
sys.path.insert(0, 'Utilities/')
import os
import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.io
from scipy.interpolate import griddata
from pyDOE import lhs
from mpl_toolkits.mplot3d import Axes3D
import time
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from tensorflow.python.client import device_lib

np.random.seed(1234)
tf.set_random_seed(1234)

#---------------------
def figsize(scale, nplots = 1):
    fig_width_pt = 390.0                          # Get this from LaTeX using \the\textwidth
    inches_per_pt = 1.0/72.27                       # Convert pt to inch
    golden_mean = (np.sqrt(5.0)-1.0)/2.0            # Aesthetic ratio (you could change this)
    fig_width = fig_width_pt*inches_per_pt*scale    # width in inches
    fig_height = nplots*fig_width*golden_mean              # height in inches
    fig_size = [fig_width,fig_height]
    return fig_size

pgf_with_latex = {                      # setup matplotlib to use latex for output
    "pgf.texsystem": "pdflatex",        # change this if using xetex or lautex
    "text.usetex": True,                # use LaTeX to write all text
    "font.family": "serif",
    "font.serif": [],                   # blank entries should cause plots to inherit fonts from the document
    "font.sans-serif": [],
    "font.monospace": [],
    "axes.labelsize": 10,               # LaTeX default is 10pt font.
    "font.size": 10,
    "legend.fontsize": 8,               # Make the legend/label fonts a little smaller
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "figure.figsize": figsize(1.0),     # default fig size of 0.9 textwidth
    "pgf.preamble": [
        r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts becasue your computer can handle it :)
        r"\usepackage[T1]{fontenc}",        # plots will be generated using this preamble
        ]
    }

def newfig(width, nplots = 1):
    fig = plt.figure(figsize=figsize(width, nplots))
    ax = fig.add_subplot(111)
    return fig, ax

def savefig(filename, crop = True):
    if crop == True:
        plt.savefig('{}.pdf'.format(filename), bbox_inches='tight', pad_inches=0)
        plt.savefig('{}.eps'.format(filename), bbox_inches='tight', pad_inches=0)
    else:
        plt.savefig('{}.pdf'.format(filename))
        plt.savefig('{}.eps'.format(filename))

class PhysicsInformedNN:
    def __init__(self, X_u, u, X_IC0, IC0, X_BC0, X_BC, BC, X_f, layers, lb, ub, rhoA, EI):

        self.layers = layers
        self.rhoA = rhoA
        self.EI = EI
	
        self.lb = lb
        self.ub = ub
		
        self.x_u = X_u[:,0:1]# variable x 
        self.t_u = X_u[:,1:2]# varialbe t
        
        self.u = u
    
        self.x_IC0 = X_IC0[:,0:1]
        self.t_IC0 = X_IC0[:,1:2]

        self.IC0 = IC0
        
        self.x_BC0 = X_BC0[:,0:1]
        self.t_BC0 = X_BC0[:,1:2]

        self.x_BC = X_BC[:,0:1]
        self.t_BC = X_BC[:,1:2]
		
        self.BC = BC/rhoA
		
        self.x_f = X_f[:,0:1]
        self.t_f = X_f[:,1:2]
		     

        
        # Initialize NNs
        self.weights, self.biases = self.initialize_NN(layers)
        
        # tf placeholders and graph
        self.sess = tf.Session(config=tf.ConfigProto(allow_soft_placement=True,
                                                     log_device_placement=True))
													 
        self.x_u_tf = tf.placeholder(tf.float32, shape=[None, self.x_u.shape[1]])
        self.t_u_tf = tf.placeholder(tf.float32, shape=[None, self.t_u.shape[1]])        
        self.u_tf = tf.placeholder(tf.float32, shape=[None, self.u.shape[1]])		
        
        self.x_IC0_tf = tf.placeholder(tf.float32, shape=[None, self.x_IC0.shape[1]])
        self.t_IC0_tf = tf.placeholder(tf.float32, shape=[None, self.t_IC0.shape[1]])        
        self.IC0_tf = tf.placeholder(tf.float32, shape=[None, self.IC0.shape[1]]) 

        self.x_BC0_tf = tf.placeholder(tf.float32, shape=[None, self.x_BC0.shape[1]])
        self.t_BC0_tf = tf.placeholder(tf.float32, shape=[None, self.t_BC0.shape[1]]) 
		
        self.x_BC_tf = tf.placeholder(tf.float32, shape=[None, self.x_BC.shape[1]])
        self.t_BC_tf = tf.placeholder(tf.float32, shape=[None, self.t_BC.shape[1]]) 
        self.BC_tf = tf.placeholder(tf.float32, shape=[None, self.BC.shape[1]]) 

        self.x_f_tf = tf.placeholder(tf.float32, shape=[None, self.x_f.shape[1]])
        self.t_f_tf = tf.placeholder(tf.float32, shape=[None, self.t_f.shape[1]])  		

        self.u_pred = self.net_u(self.x_u_tf, self.t_u_tf) 		
        self.IC0_pred1, self.IC0_pred2 = self.net_IC0(self.x_IC0_tf, self.t_IC0_tf) 
        self.BC0_pred1, self.BC0_pred2 = self.net_BC0(self.x_BC0_tf, self.t_BC0_tf)      
        self.BC1_pred1, self.BC1_pred2 = self.net_BC(self.x_BC_tf, self.t_BC_tf) 
        self.f_pred = self.net_f(self.x_f_tf, self.t_f_tf) 		
        
        self.loss = tf.reduce_mean(tf.square(self.u_tf - self.u_pred)) + \
		            tf.reduce_mean(tf.square(self.f_pred)) + \
                    tf.reduce_mean(tf.square(self.BC_tf-self.BC1_pred1)) + \
                    tf.reduce_mean(tf.square(self.BC1_pred2)) + \
                    tf.reduce_mean(tf.square(self.BC0_pred1)) + \
                    tf.reduce_mean(tf.square(self.BC0_pred2)) + \
                    tf.reduce_mean(tf.square(self.IC0_pred1-self.IC0_tf)) + \
                    tf.reduce_mean(tf.square(self.IC0_pred2)) 				
					
        self.optimizer = tf.contrib.opt.ScipyOptimizerInterface(self.loss, 
                                                                method = 'L-BFGS-B', 
                                                                options = {'maxiter': 100000,
                                                                           'maxfun': 100000,
                                                                           'maxcor': 500,
                                                                           'maxls': 500,
                                                                           'ftol' : 0.2 * np.finfo(float).eps})		
        
        init = tf.global_variables_initializer()
        self.sess.run(init)
          
    def initialize_NN(self, layers):        
        weights = []
        biases = []
        num_layers = len(layers) 
        for l in range(0,num_layers-1):
            W = self.xavier_init(size=[layers[l], layers[l+1]])
            b = tf.Variable(tf.zeros([1,layers[l+1]], dtype=tf.float32), dtype=tf.float32)
            weights.append(W)
            biases.append(b)        
        return weights, biases
        
    def xavier_init(self, size):
        in_dim = size[0]
        out_dim = size[1]        
        xavier_stddev = np.sqrt(2/(in_dim + out_dim))
        return tf.Variable(tf.truncated_normal([in_dim, out_dim], stddev=xavier_stddev), dtype=tf.float32)
    
    def neural_net(self, X, weights, biases):
        num_layers = len(weights) + 1
        
        H = 2.0*(X - self.lb)/(self.ub - self.lb) - 1.0
        for l in range(0,num_layers-2):
            W = weights[l]
            b = biases[l]
            H = tf.tanh(tf.add(tf.matmul(H, W), b))
        W = weights[-1]
        b = biases[-1]
        Y = tf.add(tf.matmul(H, W), b)
        return Y

    def net_u(self, x, t):
        u = self.neural_net(tf.concat([x,t],1), self.weights, self.biases)	
        return u	
		
    def net_IC0(self, x, t):
        u = self.net_u(x,t)
        u_t = tf.gradients(u, t)[0]		
        return u,u_t
    
    def net_BC0(self, x,t):
        u = self.net_u(x,t)
        u_x = tf.gradients(u, x)[0]
        return u, u_x	
		
    def net_BC(self, x,t):
        u = self.net_u(x,t)
        u_x = tf.gradients(u, x)[0]
        u_xx = tf.gradients(u_x, x)[0]
        u_xxx = tf.gradients(u_xx, x)[0]
        return u_xx, u_xxx	

    def net_f(self, x,t):
        u = self.net_u(x,t)
        u_t = tf.gradients(u, t)[0]
        u_tt = tf.gradients(u_t, t)[0]
        u_x = tf.gradients(u, x)[0]
        u_xx = tf.gradients(u_x, x)[0]
        u_xxx = tf.gradients(u_xx, x)[0]
        u_xxxx = tf.gradients(u_xxx, x)[0]		
        my_delta = lambda xx: 1 if xx==1 else 0
        my_time = lambda tt: 5.*tf.exp(-0.2*tt)*tf.sin(8.*tt)
        force = my_delta(x)*my_time(t)	
        f = self.rhoA*u_tt + self.EI*u_xxxx-force
        return f		
    
    def callback(self, loss):
        print('Loss:', loss)
        
    def train(self):
        
        tf_dict = {self.x_IC0_tf: self.x_IC0, self.t_IC0_tf: self.t_IC0, self.IC0_tf: self.IC0,
                   self.x_f_tf: self.x_f, self.t_f_tf: self.t_f, self.x_BC_tf: self.x_BC,
				   self.t_BC_tf: self.t_BC,self.x_BC0_tf: self.x_BC0, self.t_BC0_tf: self.t_BC0,
				   self.x_u_tf: self.x_u, self.t_u_tf: self.t_u,self.u_tf: self.u,self.BC_tf: self.BC}
                                                                                                                          
        self.optimizer.minimize(self.sess, 
                                feed_dict = tf_dict,         
                                fetches = [self.loss], 
                                loss_callback = self.callback)                                    
    
    def predict(self, X_star):
                
        u_star = self.sess.run(self.u_pred, {self.x_u_tf: X_star[:,0:1], self.t_u_tf: X_star[:,1:2]})  
        f_star = self.sess.run(self.f_pred, {self.x_f_tf: X_star[:,0:1], self.t_f_tf: X_star[:,1:2]})
               
        return u_star, f_star
    
if __name__ == "__main__": 
     
    rho = 2.375*(10**-9)
    E = 72000
    I = 416.7
    A = 200
    rhoA = rho*A
    EI = E*I
    noise = 0.0        

    N_u = 100 #number of boundray data used for training
    N_f = 10000#number of interior data used for training
    layers = [2, 45, 45, 45, 45, 45, 45, 45, 45, 1] #Node for each layer
    
    data = scipy.io.loadmat('data.mat')

    t = data['t'].flatten()[:,None]# 100 
    x = data['x'].flatten()[:,None]# 256
    Exact = np.real(data['usol'])# 100x256 analytical solution
    Force = np.real(data['ufoc'])#  force history term

    X, T = np.meshgrid(x,t)# domain meshgrid 100x256 

    X_star = np.hstack((X.flatten()[:,None], T.flatten()[:,None]))# independent variables as input X<-(x,t)
    u_star = Exact.flatten()[:,None]# dependent variable u<-             

    # Doman bounds
    lb = X_star.min(0) # 0-axis lower bound for each independent variable 
    ub = X_star.max(0) # upper bound for each independent variable
        
    xx1 = np.hstack((X[0:1,:].T, T[0:1,:].T))
    uu1 = Exact[0:1,:].T # training data @t=0
    xx2 = np.hstack((X[:,0:1], T[:,0:1]))
    uu2 = Exact[:,0:1] # training data @x=0
    xx3 = np.hstack((X[:,-1:], T[:,-1:]))
    uu3 = Exact[:,-1:] # training data @x=1
    xx4 = np.hstack((Force[:,0:1],Force[:,1:2]))
    uu4 = Force[:,-1:]
    
    X_IC0_train = np.vstack([xx1]) # all train data
    IC0_train = np.vstack([uu1])
    X_BC0_train = xx2	
    X_BC_train = xx3 # all train data	
    BC_train = uu4
    X_f_train = lb + (ub-lb)*lhs(2, N_f) # Latin hybercube design 
    X_f_train = np.vstack([X_f_train, xx1, xx2, xx3])#put the boundary point into training?

    X_u_train = np.vstack([xx1, xx2, xx3]) # for plotting
    u_train = np.vstack([uu1, uu2, uu3])
        
    model = PhysicsInformedNN(X_u_train, u_train, X_IC0_train, IC0_train, X_BC0_train, 
                              X_BC_train, BC_train, X_f_train, layers, lb, ub, rhoA, EI)
    
    start_time = time.time()                
    model.train()
    elapsed = time.time() - start_time                
    print('Training time: %.4f' % (elapsed))
    
    u_pred, f_pred = model.predict(X_star)
            
    error_u = np.linalg.norm(u_star-u_pred,2)/np.linalg.norm(u_star,2)
    print('Error u: %e' % (error_u))                     

    
    U_pred = griddata(X_star, u_pred.flatten(), (X, T), method='cubic')
    Error = np.abs(Exact - U_pred)   
    
    ######################################################################
    ############################# Plotting ###############################
    ######################################################################    
	
    fig, ax = newfig(1.0, 1.1)
    ax.axis('off')
    
    ####### Row 0: u(t,x) ##################    
    gs0 = gridspec.GridSpec(1, 2)
    gs0.update(top=1-0.06, bottom=1-1/3, left=0.15, right=0.85, wspace=0)
    ax = plt.subplot(gs0[:, :])
    
    h = ax.imshow(U_pred.T, interpolation='nearest', cmap='rainbow', 
                  extent=[t.min(), t.max(), x.min(), x.max()], 
                  origin='lower', aspect='auto')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(h, cax=cax)
    
    ax.plot(X_u_train[:,1], X_u_train[:,0], 'kx', label = 'Data (%d points)' % (u_train.shape[0]), markersize = 4, clip_on = False)
    
    line = np.linspace(x.min(), x.max(), 2)[:,None]
    ax.plot(t[125]*np.ones((2,1)), line, 'w-', linewidth = 1)
    ax.plot(t[250]*np.ones((2,1)), line, 'w-', linewidth = 1)
    ax.plot(t[375]*np.ones((2,1)), line, 'w-', linewidth = 1)    
    
    ax.set_xlabel('$t$')
    ax.set_ylabel('$x$')
    ax.legend(frameon=False, loc = 'best')
    ax.set_title('$u(t,x)$', fontsize = 10)
    
    ####### Row 1: u(t,x) slices ##################    
    gs1 = gridspec.GridSpec(1, 3)
    gs1.update(top=1-1/3, bottom=0, left=0.1, right=0.9, wspace=0.5)
    
    ax = plt.subplot(gs1[0, 0])
    ax.plot(x,Exact[125,:], 'b-', linewidth = 2, label = 'Exact')       
    ax.plot(x,U_pred[125,:], 'r--', linewidth = 2, label = 'Prediction')
    ax.set_xlabel('$x$')
    ax.set_ylabel('$u(t,x)$')    
    ax.set_title('$t = 0.5$', fontsize = 10)
    ax.axis('square')
    ax.set_xlim([0,1.1])
    ax.set_ylim([-0.4,0.4])
    
    ax = plt.subplot(gs1[0, 1])
    ax.plot(x,Exact[250,:], 'b-', linewidth = 2, label = 'Exact')       
    ax.plot(x,U_pred[250,:], 'r--', linewidth = 2, label = 'Prediction')
    ax.set_xlabel('$x$')
    ax.set_ylabel('$u(t,x)$')
    ax.axis('square')
    ax.set_xlim([0,1.1])
    ax.set_ylim([-0.4,0.4])
    ax.set_title('$t = 1.0$', fontsize = 10)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.35), ncol=5, frameon=False)
    
    ax = plt.subplot(gs1[0, 2])
    ax.plot(x,Exact[375,:], 'b-', linewidth = 2, label = 'Exact')       
    ax.plot(x,U_pred[375,:], 'r--', linewidth = 2, label = 'Prediction')
    ax.set_xlabel('$x$')
    ax.set_ylabel('$u(t,x)$')
    ax.axis('square')
    ax.set_xlim([0,1.1])
    ax.set_ylim([-0.4,0.4]) 
    ax.set_title('$t = 1.5$', fontsize = 10)
    
    savefig('./figures/beam')
    np.savetxt('./tables/pred.csv',U_pred)
    



