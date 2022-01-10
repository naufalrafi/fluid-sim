import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from FDMFuncs import *
from myClasses import *

#### SPATIAL AND TEMPORAL INPUTS
length=4 #Length of computational domain in the x-direction
breadth=4 #Breadth of computational domain in the y-direction
colpts=257 #Number of grid points in the x-direction #KEEP ODD
rowpts=257 #Number of grid points in the y-direction #KEEP ODD

#Create an object of the class Space called cavity
cavity=Space()
cavity.CreateMesh(rowpts,colpts)
cavity.SetDeltas(breadth,length)

#Create mesh for X and Y inputs to the figure
x=np.linspace(0,length,colpts)
y=np.linspace(0,breadth,rowpts)
[X,Y]=np.meshgrid(x,y)
#Determine indexing for stream plot (10 points only)
index_cut_x=int(colpts/10)
index_cut_y=int(rowpts/10)

#Create blank figure
fig=plt.figure(figsize=(16,8))
ax=plt.axes(xlim=(0,length),ylim=(0,breadth))

#### FLUID PROPERTIES
rho=1 #Density of fluid
mu=0.01 #Viscosity of fluid
#Create an object of the class Fluid called water
water=Fluid(rho,mu)

#### BOUNDARY SPECIFICATIONS
u_in=1 #Lid velocity
v_wall=0 #Velocity of fluid at the walls
p_out=0 #Gauge pressure at the boundaries
#Create objects of the class Boundary having either Dirichlet ("D") or Neumann ("N") type boundaries
flow=Boundary("D",u_in)
noslip=Boundary("D",v_wall)
zeroflux=Boundary("N",0)
pressureatm=Boundary("D",p_out)

#### Simulation inputs
rowpts=257
colpts=257
length=4
breadth=4

#### SIMULATION PARAMETERS
time=150 #Simulation time
CFL_number=0.8 #Reduce this if solution diverges
file_flag=1 #Keep 1 to print results to file
interval=100 #Record values in file per interval number of iterations

#Go to the Result directory
cwdir=os.getcwd()
dir_path=os.path.join(cwdir,"Result")
os.chdir(dir_path)

#Go through files in the directory and store filenames
filenames=[]
iterations=[]
for root,dirs,files in os.walk(dir_path):
    for datafile in files:
        if "PUV" in datafile:
            filenames.append(datafile)
            no_ext_file=datafile.replace(".txt","").strip()
            iter_no=int(no_ext_file.split("V")[-1])
            iterations.append(iter_no)

#Discern the final iteration and interval
initial_iter=np.amin(iterations)            
final_iter=np.amax(iterations)
inter=(final_iter - initial_iter)/(len(iterations)-1)
number_of_frames=len(iterations)
sorted_iterations=np.sort(iterations)


def MakeResultDirectory(wipe = False):
    cwdir=os.getcwd()
    dir_path=os.path.join(cwdir,"Result")
    if not os.path.isdir(dir_path):
        os.makedirs(dir_path,exist_ok=True)
    else:
        if wipe:
            os.chdir(dir_path)
            filelist=os.listdir()
            for file in filelist:
                os.remove(file)
    
    os.chdir(cwdir)


def read_datafile(iteration):
    #Set filename and path according to given iteration
    filename="PUV{0}.txt".format(iteration)
    filepath=os.path.join(dir_path,filename)
    #Load text file as numpy array
    arr=np.loadtxt(filepath,delimiter="\t")
    rows,cols=arr.shape
    #Define empty arrays for pressure and velocities
    p_p=np.zeros((rowpts,colpts))
    u_p=np.zeros((rowpts,colpts))
    v_p=np.zeros((rowpts,colpts))
    #Organize imported array into variables
    p_arr=arr[:,0]
    u_arr=arr[:,1]
    v_arr=arr[:,2]
    
    #Reshape 1D data into 2D
    p_p=p_arr.reshape((rowpts,colpts))
    u_p=u_arr.reshape((rowpts,colpts))
    v_p=v_arr.reshape((rowpts,colpts))
    
    return p_p,u_p,v_p
            
def animate(i):
    #Print frames left to be added to the animation
    sys.stdout.write("\rFrames remaining: {0:03d}".format(len(sorted_iterations)-i))
    sys.stdout.flush()
    #Get iterations in a sequential manner through sorted_iterations
    iteration=sorted_iterations[i]
    #Use the read_datafile function to get pressure and velocities
    p_p,u_p,v_p=read_datafile(iteration)
    #Clear previous plot and make contour and stream plots for current iteration
    ax.clear()
    ax.set_xlim([0,length])
    ax.set_ylim([0,breadth])
    ax.set_xlabel("$x$",fontsize=12)
    ax.set_ylabel("$y$",fontsize=12)
    ax.set_title("Frame No: {0}".format(i))
    cont=ax.contourf(X,Y,p_p)
    stream=ax.streamplot(X[::index_cut_y,::index_cut_x],\
                         Y[::index_cut_y,::index_cut_x],\
                         u_p[::index_cut_y,::index_cut_x],\
                         v_p[::index_cut_y,::index_cut_x],\
                         color="k")
    return cont,stream


def WriteToFile(space: mc.Space, iteration, interval):
    if(iteration%interval==0):
        dir_path=os.path.join(os.getcwd(),"Result")
        filename="PUV{0}.txt".format(iteration)
        path=os.path.join(dir_path,filename)
        with open(path,"w") as f:
            for i in range(space.rowpts):
                for j in range(space.colpts):
                    f.write("{}\t{}\t{}\n".format(space.p_c[i,j],space.u_c[i,j],space.v_c[i,j]))


#### RUN SIMULATION
# Print general simulation information
print("######## Beginning FlowPy Simulation ########")
print("#############################################")
print("# Simulation time: {0:.2f}".format(time))
print("# Mesh: {0} x {1}".format(colpts,rowpts))
print("# Re/u: {0:.2f}\tRe/v:{1:.2f}".format(rho*length/mu,rho*breadth/mu))
print("# Save outputs to text file: {0}".format(bool(file_flag)))
## Initialization
# Make directory to store results
MakeResultDirectory(wipe=True)
# Initialize counters
t=0
i=0
## Run
while(t<time):
    #Print time left
    sys.stdout.write("\rSimulation time left: {0:.2f}".format(time-t))
    sys.stdout.flush()
    #Set the time-step
    SetTimeStep(CFL_number,cavity,water)
    timestep=cavity.dt
    
    #Set boundary conditions
    SetUBoundary(cavity,noslip,noslip,flow,noslip)
    SetVBoundary(cavity,noslip,noslip,noslip,noslip)
    SetPBoundary(cavity,zeroflux,zeroflux,pressureatm,zeroflux)
    
    #Calculate starred velocities
    GetStarredVelocities(cavity,water)
    
    #Solve the pressure Poisson equation
    SolvePressurePoisson(cavity,water,zeroflux,zeroflux,\
pressureatm,zeroflux)
    #Solve the momentum equation
    SolveMomentumEquation(cavity,water)
    #Save variables and write to file
    SetCentrePUV(cavity)
    if(file_flag==1):
        WriteToFile(cavity,i,interval)
    #Advance time-step and counter
    t+=timestep
    i+=1

print("######## Making FlowPy Animation ########")
print("#########################################")
anim=animation.FuncAnimation(fig,animate,frames=number_of_frames,interval=50,blit=False)
movie_path=os.path.join(dir_path,"FluidFlowAnimation.mp4")
anim.save(r"{0}".format(movie_path))
print("\nAnimation saved as FluidFlowAnimation.mp4 in Result")